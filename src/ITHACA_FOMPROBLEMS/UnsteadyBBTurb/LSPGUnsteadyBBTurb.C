/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
-------------------------------------------------------------------------------

  License
  This file is part of ITHACA-FV

  ITHACA-FV is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ITHACA-FV is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with ITHACA-FV. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the LSPGUnsteadyBBTurb class.

#include "LSPGUnsteadyBBTurb.H"
#include "viscosityModel.H"
#include "alphatJayatillekeWallFunctionFvPatchScalarField.H" // Used to implement BCs
#include "calculatedFvPatchField.H" // Used to implement BCs
#include "fvCFD.H"
#include <functional>
#include <cmath>
#include "ITHACAPOD.H"
#include "pisoControl.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //
LSPGUnsteadyBBTurb::LSPGUnsteadyBBTurb(std::shared_ptr<ITHACAcontext> context)
: LSPGProblem(std::move(context))
{
  fvMesh& mesh = this->mesh();
  Time& runTime = this->runTime();
  pimpleControl& pimple = this->pimple();
  fv::options& fvOptions = this->fvOptions();
  #include "createFields.H"
  offline = ITHACAutilities::check_off();
  podex = ITHACAutilities::check_pod();
  readITHACAdict();
}

void LSPGUnsteadyBBTurb::readITHACAdict()
{
    Time& runTime = this->runTime();
    fvMesh& mesh = this->mesh();
    IOMRFZoneList& MRF = this->MRF();

    ITHACAdict = new IOdictionary(
        IOobject(
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE));
    // Are these calls necessary? Some of these also present in the constructor of unsteadyNS. Check please.
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "None");
    timeDependentBC = ITHACAdict->lookupOrDefault<bool>("timeDependentBC", false);
    derivativeInRBF = ITHACAdict->lookupOrDefault<bool>("derivativeInRBF", false);
    centerSnapshots = ITHACAdict->lookupOrDefault<bool>("centerSnapshots", false);
        
    M_Assert(bcMethod == "lift" || bcMethod == "penalty"
             || bcMethod == "Gunzburger",
                         "The BC method must be set to lift, penalty or Gunzburger in ITHACAdict");
    viscDict = ITHACAdict->subDict("viscDict");
    
    NUmodes = ITHACAdict->lookupOrDefault<label>("NmodesUproj", 10);
    NTmodes = ITHACAdict->lookupOrDefault<label>("NmodesTproj", 5);
    NPrghmodes = ITHACAdict->lookupOrDefault<label>("NmodesPrghproj", 5);
    NNutModes = ITHACAdict->lookupOrDefault<label>("NmodesNutproj", 5);

    dimInputRBF = ITHACAdict->lookupOrDefault<label>("dimInputRBF", 0);
    M_Assert(dimInputRBF <= NUmodes ,
             "The dimension of the input to the RBF must be less than or equal to the number of velocity modes.");
    M_Assert(dimInputRBF >= 0,
             "The dimension of the input to the RBF must be greater than or equal to zero.");
    Info << "### INFO ### " << nl
         << "BC method: " << bcMethod << nl
         << "Time dependent BCs: " << timeDependentBC << nl
         << "Derivative in RBF for nut interpolation: " << derivativeInRBF << nl
         << "Number of velocity modes for projection: " << NUmodes << nl
         << "Number of temperature modes for projection: " << NTmodes << nl
         << "Number of pressure modes for projection: " << NPrghmodes << nl
         << "Number of eddy viscosity modes for projection: " << NNutModes << endl;
}

// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //
void LSPGUnsteadyBBTurb::truthSolve(const List<scalar> mu_now, label nSample)
{
    fvMesh& mesh = this->mesh();
    Time& runTime = this->runTime();
    pimpleControl& pimple = this->pimple();
    #include "initContinuityErrs.H"
    fv::options& fvOptions = this->fvOptions();
    IOMRFZoneList& MRF = this->MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& p_rgh = _p_rgh();
    volScalarField& T = _T();
    volScalarField& nut = _nut();
    volScalarField& alphat = _alphat();
    volScalarField& rhok = _rhok();
    volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();
    surfaceScalarField& phi = _phi();
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();

    instantList Times = runTime.times();

    runTime.setEndTime(finalTime);
    runTime.setTime(startTime, 0);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime + writeEvery;
    label nSavedTimesteps = (finalTime - startTime) / writeEvery;
    M_Assert(timeSnapshots.size() > nSample,
             "The timeSnapshots list does not have enough space for the current sample index.");
    timeSnapshots[nSample].resize(nSavedTimesteps);
    label stepCounter = 0;
    Info << "Starting time loop." << nl << endl;

    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime++;
        // runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
#include "UEqn.H"
#include "TEqn.H"

            while (pimple.correct())
            {
#include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            nut = turbulence->nut();
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p_rgh, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(nut, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                             runTime.timeName());
            timeSnapshots[nSample](stepCounter) = runTime.value();
            stepCounter++;
            Ufield.append(U.clone());
            Prghfield.append(p_rgh.clone());
            Tfield.append(T.clone());
            Nutfield.append(nut.clone());
            nextWrite += writeEvery;
            writeMu(mu_now);
            counter++;
        }
    }
}

void LSPGUnsteadyBBTurb::truthSolve(fileName folder)
{
    fvMesh& mesh = this->mesh();
    Time& runTime = this->runTime();
    pimpleControl& pimple = this->pimple();
    fv::options& fvOptions = this->fvOptions();
    IOMRFZoneList& MRF = this->MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    
    scalar cumulativeContErr = 0.0;
    #include "initContinuityErrs.H"
    
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& p_rgh = _p_rgh();
    volScalarField& T = _T();
    volScalarField& nut = _nut();
    volScalarField& alphat = _alphat();
    volScalarField& rhok = _rhok();
    volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();
    surfaceScalarField& phi = _phi();
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();
    
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(startTime, 0);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime + writeEvery;

    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime++;
        // runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;

        while (pimple.loop())
        {
#include "UEqn.H"
#include "TEqn.H"

            while (pimple.correct())
            {
#include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;

        if (checkWrite(runTime))
        {
            nut = turbulence->nut();
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(p_rgh, name(counter), folder);
            ITHACAstream::exportSolution(T, name(counter), folder);
            ITHACAstream::exportSolution(nut, name(counter), folder);
            std::ofstream of(folder + "/" + name(counter) + "/" + runTime.timeName());
            Ufield.append(U.clone());
            Prghfield.append(p_rgh.clone());
            Tfield.append(T.clone());
            Nutfield.append(nut.clone());
            nextWrite += writeEvery;
            counter++;
        }
    }
}

void LSPGUnsteadyBBTurb::removeMean()
{
    Info << "Not implemented yet" << endl;

}

// * * * * * * * * * * * * * * RBF Prep Methods * * * * * * * * * * * * * * //
void LSPGUnsteadyBBTurb::offlineRBFInterpolation()
{
    Eigen::MatrixXd weights;
    Eigen::MatrixXd coeffL2nut = ITHACAutilities::getCoeffs(fluctNutfield, nutmodes,
        NNutModes);
    Eigen::MatrixXd coeffL2vel;
    coeffL2vel.resize(0, 0);
    label inputModes = dimInputRBF;

    if (inputModes == 0)
    {
        inputModes = NUmodes;
    }

    if (bcMethod == "lift")
    {
        coeffL2vel = ITHACAutilities::getCoeffs(Uomfield, Umodes,
                                                inputModes); // Returns a [modes x snapshots]
    }
    else
    {
        coeffL2vel = ITHACAutilities::getCoeffs(Ufield, Umodes,
                                                inputModes); // Returns a [modes x snapshots]
    }

    Info << "Shape of the L2 velocity coeff matrix: " << coeffL2vel.rows() << " x "
         << coeffL2vel.cols() << endl;
    Info << "Shape of the L2 eddy viscosity coeff matrix: " << coeffL2nut.rows() <<
            " x " << coeffL2nut.cols() << endl;
    List<Eigen::MatrixXd> velDerCoeff(2);

    // Returns a list of two matrices: [0] = velocity derivative coeffs, [1] = eddy viscosity coeffs. Each matrix is [snapshots x coeffs]
    if (derivativeInRBF)
    {
        velDerCoeff = velDerivativeCoeff(coeffL2vel.transpose(), coeffL2nut.transpose(),
                                         timeSnapshots);
        dimA = velDerCoeff[0].cols();
    }
    else
    {
        velDerCoeff[0] = coeffL2vel.transpose();
        velDerCoeff[1] = coeffL2nut.transpose();
        dimA = velDerCoeff[0].cols();
    }

    if (Pstream::master())
    {
        ITHACAutilities::createSymLink("./ITHACAoutput/Debug");
        ITHACAstream::exportMatrix(velDerCoeff[0], "A_RBF", "numpy",
                                   "./ITHACAoutput/Debug/");
        ITHACAstream::exportMatrix(velDerCoeff[1], "G_RBF", "numpy",
                                   "./ITHACAoutput/Debug/");
    }

    if (ITHACAutilities::check_file("./ITHACAoutput/shapeParameters"))
    {
        Info << "### MESSAGE - Reading RBF shape parameters from file is not yet implemented for mathtoolbox usage. Using the one provided in the dictionary."
             << endl;
    }

    rbfSplines.resize(NNutModes);
    Eigen::MatrixXd x = velDerCoeff[0].transpose();

    for (label i = 0; i < NNutModes; i++)
    {
        // Create a RBF interpolator instance
        rbfSplines[i] = std::make_shared<ithacaInterpolator>(viscDict);
        Eigen::VectorXd y = velDerCoeff[1].col(i);
        // rbfSplines[i]->optimizeShapeParameter(x, y, 5);
        rbfSplines[i]->fit(x, y);
        Info << "### INTERPOLATION - Fitting ithacaInterpolator for mode " << i + 1 <<
             " completed. Here's some information:" << endl;
        rbfSplines[i]->printInfo();
    }
}

List<Eigen::MatrixXd> LSPGUnsteadyBBTurb::velDerivativeCoeff(
    const Eigen::MatrixXd& A, const Eigen::MatrixXd& G,
    const List<Eigen::VectorXd>& snapshotTimes)
{
    const label velCoeffsNum = A.cols();
    const label parsSamplesNum = snapshotTimes.size();
    // 1. Calculate total rows without storing per-sample counts in a vector
    label newRowsNum = 0;

    for (label j = 0; j < parsSamplesNum; ++j)
    {
        newRowsNum += (snapshotTimes[j].size() - 1);
    }

    List<Eigen::MatrixXd> newCoeffs;
    newCoeffs.setSize(2);
    const label newColsNum = 2 * velCoeffsNum;
    newCoeffs[0].resize(newRowsNum, newColsNum);
    newCoeffs[1].resize(newRowsNum, G.cols());
    label outOffset = 0;
    label blockStart = 0;
    // Pre-alias sub-matrices for readability and safety
    Eigen::MatrixXd& C0 = newCoeffs[0];
    Eigen::MatrixXd& C1 = newCoeffs[1];

    for (label j = 0; j < parsSamplesNum; ++j)
    {
        const Eigen::VectorXd& timeSnap = snapshotTimes[j];
        const label nSnaps = timeSnap.size();
        const label rowsPerBlock = nSnaps - 1;

        if (rowsPerBlock <= 0)
        {
            continue;
        }

        // Direct reference blocks into matrix A (zero allocation)
        auto b0 = A.middleRows(blockStart, rowsPerBlock);
        auto b2 = A.middleRows(blockStart + 1, rowsPerBlock);
        // Compute time steps: timeSnap[k+1] - timeSnap[k]
        // Uses Eigen block expression to avoid tail/head array copies
        auto deltaT = timeSnap.tail(rowsPerBlock) - timeSnap.head(rowsPerBlock);
        // Direct block assignment for newCoeffs[0] left half (b2)
        C0.block(outOffset, 0, rowsPerBlock, velCoeffsNum) = b2;
        // Direct vectorized calculation into newCoeffs[0] right half (derivative)
        // Eliminates intermediate matrices b0, b2, derivative, and bNew completely
        C0.block(outOffset, velCoeffsNum, rowsPerBlock, velCoeffsNum) =
              (b2 - b0).array().colwise() / deltaT.array();
        // Direct block assignment for newCoeffs[1]
        C1.middleRows(outOffset, rowsPerBlock) = G.middleRows(blockStart + 1,
            rowsPerBlock);
        outOffset += rowsPerBlock;
        blockStart += nSnaps;
    }

    return newCoeffs;
}

void LSPGUnsteadyBBTurb::splitEddyViscositySnapshots()
{
    const label nSamples = timeSnapshots.size();
    avgNutfield.setSize(nSamples);
    const label totalSnapshots = Nutfield.size();
    fluctNutfield.setSize(totalSnapshots);
    label globalIndex = 0;
    label flatIndex = 0;

    for (label i = 0; i < nSamples; i++)
    {
        label nSnap = timeSnapshots[i].size();
        M_Assert(nSnap > 0, "Each parameter sample must have at least one snapshot");
        // 1. Calculate Average
        // Initialize with the first snapshot of the current sample
        avgNutfield.set(i, new volScalarField(
                            IOobject(
                                "avgNut",
                                Nutfield[globalIndex].time().timeName(),
                                Nutfield[globalIndex].mesh(),
                                IOobject::NO_READ,
                                IOobject::NO_WRITE),
                            Nutfield[globalIndex]));
        volScalarField& avg = avgNutfield[i];

        for (label j = 1; j < nSnap; j++)
        {
            avg += Nutfield[globalIndex + j];
        }

        avg /= scalar(nSnap);

        // 2. Calculate Fluctuations directly
        for (label j = 0; j < nSnap; j++)
        {
            fluctNutfield.set(flatIndex,
                              new volScalarField(Nutfield[globalIndex + j] - avg));
            fluctNutfield[flatIndex].rename("fluctNut");
            flatIndex++;
        }

        globalIndex += nSnap;
    }
}

// * * * * * * * * * * * * * * Preparation of modes * * * * * * * * * *  //
void LSPGUnsteadyBBTurb::prepareModes()
{
    if (centerSnapshots)
    {
        Umodes.append(UavgPtr->clone());
        Tmodes.append(TavgPtr->clone());
        // Add the Prgh mean field to the top of the volScalarModes (derivated from autoPtr<volScalarField>), without substitutuing
        int oldPsize = Prghmodes.size();
        Prghmodes.resize(oldPsize + 1);

        for (label i = oldPsize; i > 0; i--)
        {
            Prghmodes.set(i, Prghmodes.release(i - 1));
        }

        Prghmodes.set(0, p_rghavgPtr->clone());
    }
}

// * * * * * * * * * * * * * * Restart Methods * * * * * * * * * * * * * //

void LSPGUnsteadyBBTurb::restart()
{
    Time& runTime = this->runTime();
    fvMesh& mesh = this->mesh();
    pimpleControl& pimple = this->pimple();

    runTime.setTime(0, 0);
    // Read transportProperties dictionary
    IOdictionary transportProperties(
        IOobject(
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE));
    // Update constants
    _nu() = dimensionedScalar("nu", dimViscosity, transportProperties.lookup("nu"));
    _Pr() = dimensionedScalar("Pr", dimless, transportProperties.lookup("Pr"));
    _Prt() = dimensionedScalar("Prt", dimless, transportProperties.lookup("Prt"));
    _beta() = dimensionedScalar("beta", dimless / dimTemperature,
                                transportProperties.lookup("beta"));
    _TRef() = dimensionedScalar("TRef", dimTemperature,
                                transportProperties.lookup("TRef"));
    // Now we reset the fields, reading from disk (OpenFOAM folders)
    volVectorField U_new(IOobject("U", runTime.timeName(), mesh,
                                  IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _U() = U_new;
    // _U().correctBoundaryConditions();
    volScalarField T_new(IOobject("T", runTime.timeName(), mesh,
                                  IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _T() = T_new;
    // _T().correctBoundaryConditions();
    volScalarField p_rgh_new(IOobject("p_rgh", runTime.timeName(), mesh,
                                      IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _p_rgh() = p_rgh_new;
    // _p_rgh().correctBoundaryConditions();
    volScalarField nut_new(IOobject("nut", runTime.timeName(), mesh,
                                    IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _nut() = nut_new;
    // _nut().correctBoundaryConditions();
    volScalarField alphat_new(IOobject("alphat", runTime.timeName(), mesh,
                                       IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _alphat() = alphat_new;
    // _alphat().correctBoundaryConditions();
    volScalarField UliftBC_new(IOobject("UliftBC", runTime.timeName(), mesh,
                                        IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _UliftBC() = UliftBC_new;
    _phi() = linearInterpolate(_U()) & mesh.Sf();
    _rhok() = 1.0 - _beta() * (_T() - _TRef());
    _laminarTransport.clear();
    _laminarTransport = autoPtr<singlePhaseTransportModel>(new
        singlePhaseTransportModel(_U(), _phi()));
    turbulence.clear();
    turbulence = autoPtr<incompressible::turbulenceModel>
                 (incompressible::turbulenceModel::New(_U(), _phi(), _laminarTransport()));
    turbulence->validate();
    _p() = _p_rgh() + _rhok() * _gh();
    setRefCell(_p(), _p_rgh(), pimple.dict(), pRefCell, pRefValue);

    if (_p_rgh().needReference())
    {
        _p() += dimensionedScalar("p", _p().dimensions(),
                                  pRefValue - getRefCellValue(_p(), pRefCell));
    }

    mesh.setFluxRequired(_p_rgh().name());
    Info << "Restart complete." << endl;
}

void LSPGUnsteadyBBTurb::resizeModes()
{
    Umodes.resize(NUmodes);
    Tmodes.resize(NTmodes);
    nutmodes.resize(NNutModes);
    Prghmodes.resize(NPrghmodes);
}

void LSPGUnsteadyBBTurb::computePOD(label nModesU, label nModesPrgh, label nModesT, label nModesNut)
{
  if (nModesU == 0)
  {
    nModesU = NUmodes;
  }
  if (nModesPrgh == 0)
  {
    nModesPrgh = NPrghmodes;
  }
  if (nModesT == 0)
  {
    nModesT = NTmodes;
  }
  if (nModesNut == 0)
  {
    nModesNut = NNutModes;
  }
  Info << "### MESSAGE - Computing POD modes for velocity, pressure, temperature and eddy viscosity." << endl;
  if (bcMethod == "lift")
  {
    ITHACAPOD::getModes(Uomfield, Umodes, _U().name(), podex, 0, 0, nModesU, true);
    ITHACAPOD::getModes(Prghfield, Prghmodes, _p_rgh().name(), podex, 0, 0, nModesPrgh, false);
    ITHACAPOD::getModes(Tomfield, Tmodes, _T().name(), podex, 0, 0, nModesT, true);
    ITHACAPOD::getModes(fluctNutfield, nutmodes, fluctNutfield[0].name(), podex, 0, 0, nModesNut, true);
  }
  else
  {
    ITHACAPOD::getModes(Ufield, Umodes, _U().name(), podex, 0, 0, nModesU, true);
    ITHACAPOD::getModes(Prghfield, Prghmodes, _p_rgh().name(), podex, 0, 0, nModesPrgh, false);
    ITHACAPOD::getModes(Tfield, Tmodes, _T().name(), podex, 0, 0, nModesT, true);
    ITHACAPOD::getModes(fluctNutfield, nutmodes, fluctNutfield[0].name(), podex, 0, 0, nModesNut, true);
  }
}

void LSPGUnsteadyBBTurb::setupLift()
{
  liftSolve();
  liftSolveT();
  computeLift(Ufield, liftfield, Uomfield);
  computeLiftT(Tfield, liftfieldT, Tomfield);
  ITHACAstream::exportFields(liftfield, "./ITHACAoutput/Lift", "ULift");
  ITHACAstream::exportFields(liftfieldT, "./ITHACAoutput/Lift", "TLift");
  // ITHACAstream::exportFields(Uomfield, "./ITHACAoutput/Lift", "U_om");
  // ITHACAstream::exportFields(Tomfield, "./ITHACAoutput/Lift", "T_om");
}

void LSPGUnsteadyBBTurb::liftSolve()
{
    for (label k = 0; k < inletIndex.rows(); k++)
    {
        Time& runTime = this->runTime();
        fvMesh& mesh = this->mesh();
        IOMRFZoneList& MRF = this->MRF();
        
        surfaceScalarField& phi = _phi();
        pisoControl potentialFlow(mesh, "potentialFlow");
        volVectorField& U = _U();
        volScalarField& UliftBC = _UliftBC();
        label BCind = inletIndex(k, 0);
        volVectorField Ulift("Ulift" + name(k), U);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        Info << "Solving a lifting Problem" << endl;
        Vector<double> v1(0, 0, 0);
        v1[inletIndex(k, 1)] = 1;
        Vector<double> v0(0, 0, 0);

        for (label j = 0; j < U.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Ulift, j, v1);
            }
            else if (U.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Ulift, j, v0);
            }
            else
            {
            }

            assignIF(Ulift, v0);
            phi = linearInterpolate(Ulift) & mesh.Sf();
        }

        Info << "Constructing velocity potential field Phi\n" << endl;
        volScalarField Phi
        (
            IOobject
            (
                "Phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Phi", dimLength * dimVelocity, 0),
            UliftBC.boundaryField().types()
        );
        label PhiRefCell = 0;
        scalar PhiRefValue = 0.0;
        setRefCell
        (
            Phi,
            potentialFlow.dict(),
            PhiRefCell,
            PhiRefValue
        );
        mesh.setFluxRequired(Phi.name());
        runTime.functionObjects().start();
        MRF.makeRelative(phi);
        adjustPhi(phi, Ulift, UliftBC);

        while (potentialFlow.correctNonOrthogonal())
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(dimensionedScalar("1", dimless, 1), Phi)
                ==
                fvc::div(phi)
            );
            PhiEqn.setReference(PhiRefCell, PhiRefValue);
            PhiEqn.solve();

            if (potentialFlow.finalNonOrthogonalIter())
            {
                phi -= PhiEqn.flux();
            }
        }

        MRF.makeAbsolute(phi);
        Info << "Continuity error = "
             << mag(fvc::div(phi))().weightedAverage(mesh.V()).value()
             << endl;
        Ulift = fvc::reconstruct(phi);
        Ulift.correctBoundaryConditions();
        Info << "Interpolated velocity error = "
             << (sqrt(sum(sqr((fvc::interpolate(U) & mesh.Sf()) - phi)))
                 / sum(mesh.magSf())).value()
             << endl;
        Ulift.write();
        liftfield.append(Ulift.clone());
    }
}

void LSPGUnsteadyBBTurb::liftSolveT()
{
    for (label k = 0; k < inletIndexT.rows(); k++)
    {
        Time& runTime = this->runTime();
        fvMesh& mesh = this->mesh();

        volScalarField& T = _T();
        volVectorField& U = _U();
        surfaceScalarField& phi = _phi();
        phi = linearInterpolate(U) & mesh.Sf();
        simpleControl simple(mesh);
        volScalarField& alphat = _alphat();

        dimensionedScalar& Pr = _Pr();
        dimensionedScalar& Prt = _Prt();
        label BCind = inletIndexT(k, 0);
        volScalarField Tlift("Tlift" + name(k), T);
        instantList Times = runTime.times();
        runTime.setTime(Times[1], 1);
        Info << "Solving a lifting Problem" << endl;
        scalar t1 = 1;
        scalar t0 = 0;
        alphat = turbulence->nut() / Prt;
        alphat.correctBoundaryConditions();
        volScalarField alphaEff("alphaEff", turbulence->nu() / Pr + alphat);

        for (label j = 0; j < T.boundaryField().size(); j++)
        {
            if (j == BCind)
            {
                assignBC(Tlift, j, t1);
                assignIF(Tlift, t0);
            }
            else if (T.boundaryField()[BCind].type() == "fixedValue")
            {
                assignBC(Tlift, j, t0);
                assignIF(Tlift, t0);
            }
            else
            {
            }
        }

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::div(phi, Tlift)
                - fvm::laplacian(alphaEff, Tlift)
            );
            TEqn.solve();
            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
        }

        Tlift.write();
        liftfieldT.append(Tlift.clone());
    }
}
