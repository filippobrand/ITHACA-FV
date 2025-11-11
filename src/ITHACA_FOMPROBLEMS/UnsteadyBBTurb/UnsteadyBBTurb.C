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
/// Source file of the UnsteadyBBTurb class.

#include "UnsteadyBBTurb.H"
#include "viscosityModel.H"
#include "alphatJayatillekeWallFunctionFvPatchScalarField.H" // Used to implement BCs
#include "calculatedFvPatchField.H" // Used to implement BCs
#include <cmath>

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
UnsteadyBBTurb::UnsteadyBBTurb() { }

UnsteadyBBTurb::UnsteadyBBTurb(int argc, char* argv[])
{
    _args = autoPtr<argList>(
        new argList(argc, argv));

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _pimple = autoPtr<pimpleControl>(
        new pimpleControl(
            mesh));
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "createFields.H"
#pragma GCC diagnostic pop
#include "createFvOptions.H"
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //
#include "fvCFD.H"

// Method to performa a truthSolve
void UnsteadyBBTurb::truthSolve(List<scalar> mu_now)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
#include "initContinuityErrs.H"
    fv::options& fvOptions = _fvOptions();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    volScalarField& p = _p();
    volVectorField& U = _U();
    volScalarField& p_rgh = _p_rgh();
    volScalarField& T = _T();
    volScalarField _nut(turbulence->nut());
    volScalarField& alphat = _alphat();
    volScalarField& rhok = _rhok();
    volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();
    surfaceScalarField& phi = _phi();
    pimpleControl& pimple = _pimple();
    IOMRFZoneList& MRF = _MRF();
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    turbulence->validate();
    // save initial condition in folder 0
    ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(p_rgh, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
    ITHACAstream::exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
    std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
        runTime.timeName());
    Ufield.append(U.clone());
    Pfield.append(p.clone());
    Prghfield.append(p_rgh.clone());
    Tfield.append(T.clone());
    Nutfield.append(_nut.clone());
    counter++;
    nextWrite += writeEvery;
    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
#include "UEqn.H"
#include "TEqn.H"
            // --- Pressure corrector loop
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
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p_rgh, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(_nut, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" +
                runTime.timeName());
            Ufield.append(U.clone());
            Pfield.append(p.clone());
            Prghfield.append(p_rgh.clone());
            Tfield.append(T.clone());
            Nutfield.append(_nut.clone());
            counter++;
            nextWrite += writeEvery;
            writeMu(mu_now);
        }

        runTime++;
    }
}

void UnsteadyBBTurb::truthSolve(fileName folder)
{
#include "initContinuityErrs.H"
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
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
    fv::options& fvOptions = _fvOptions();
    pimpleControl& pimple = _pimple();
    IOMRFZoneList& MRF = _MRF();
    singlePhaseTransportModel& laminarTransport = _laminarTransport();
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(Times[1], 1);
    runTime.setDeltaT(timeStep);
    nextWrite = startTime;
    // Save initial condition
    ITHACAstream::exportSolution(U, name(counter), folder);
    ITHACAstream::exportSolution(p, name(counter), folder);
    ITHACAstream::exportSolution(T, name(counter), folder);
    ITHACAstream::exportSolution(p_rgh, name(counter), folder);
    ITHACAstream::exportSolution(nut, name(counter), folder);
    std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
    counter++;
    nextWrite += writeEvery;

    // Start the time loop
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;
        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
#include "UEqn.H"
#include "TEqn.H"
            // --- Pressure corrector loop
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
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(T, name(counter), folder);
            ITHACAstream::exportSolution(p_rgh, name(counter), folder);
            ITHACAstream::exportSolution(nut, name(counter), folder);
            std::ofstream of(folder + name(counter) + "/" + runTime.timeName());
            counter++;
            nextWrite += writeEvery;
        }

        runTime++;
    }
}

// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void UnsteadyBBTurb::projectSUP(fileName folder, label NUmodes, label NPrghmodes, label NTmodes,
    label NSUPmodes, label Nnutmodes)
{

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word M_str = "M_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
        {
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
        } else
        {
            M_matrix = mass_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word B_str = "B_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        } else
        {
            B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word C_str = "C_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_str))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_str);
        } else
        {
            C_tensor = convective_term_tens(NUmodes, NPrghmodes, NSUPmodes);
        }

        word K_str = "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPrghmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        } else
        {
            K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word H_str = "H_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + H_str))
        {
            ITHACAstream::ReadDenseMatrix(H_matrix, "./ITHACAoutput/Matrices/", H_str);
        } else
        {
            H_matrix = buoyant_term(NUmodes, NPrghmodes, NSUPmodes);
        }

        word W_str = "W_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + W_str))
        {
            ITHACAstream::ReadDenseMatrix(W_matrix, "./ITHACAoutput/Matrices/", W_str);
        } else
        {
            W_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
        }

        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        word Y_str = "Y_" + name(liftfieldT.size()) + "_" + name(NTmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Y_str))
        {
            ITHACAstream::ReadDenseMatrix(Y_matrix, "./ITHACAoutput/Matrices/", Y_str);
        } else
        {
            Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        }

        word P_str = "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPrghmodes);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
        } else
        {
            P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);
        }
    } else
    {
        Info << "No matrices found, assembling..." << endl; // DEBUG
        L_U_SUPmodes.resize(0);

        if (liftfield.size() != 0)
        {
            Info << "Adding lift modes..." << endl; // DEBUG
            for (label k = 0; k < liftfield.size(); k++)
            {
                L_U_SUPmodes.append(liftfield[k].clone());
            }
        }

        if (NUmodes != 0)
        {
            for (label k = 0; k < NUmodes; k++)
            {
                L_U_SUPmodes.append(Umodes[k].clone());
            }
        }

        if (NSUPmodes != 0)
        {
            for (label k = 0; k < NSUPmodes; k++)
            {
                L_U_SUPmodes.append(supmodes[k].clone());
            }
        }

        L_T_modes.resize(0);

        if (liftfieldT.size() != 0)
        {
            for (label k = 0; k < liftfieldT.size(); k++)
            {
                L_T_modes.append(liftfieldT[k].clone());
            }
        }

        if (NTmodes != 0)
        {
            for (label k = 0; k < NTmodes; k++)
            {
                L_T_modes.append(Tmodes[k].clone());
            }
        }
        M_matrix = mass_term(NUmodes, NPrghmodes, NSUPmodes);
        B_matrix = diffusive_term(NUmodes, NPrghmodes, NSUPmodes);
        C_tensor = convective_term_tens(NUmodes, NPmodes, NSUPmodes);
        K_matrix = pressure_gradient_term(NUmodes, NPrghmodes, NSUPmodes);
        H_matrix = buoyant_term(NUmodes, NTmodes, NSUPmodes);
        W_matrix = mass_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Q_matrix = convective_term_temperature(NUmodes, NTmodes, NSUPmodes);
        Y_matrix = diffusive_term_temperature(NUmodes, NTmodes, NSUPmodes);
        P_matrix = divergence_term(NUmodes, NPrghmodes, NSUPmodes);
        BT_matrix = BTturbulence(NUmodes, NSUPmodes);
        CT1_matrix = turbulenceTerm1(NUmodes, NSUPmodes, Nnutmodes);
        CT2_matrix = turbulenceTerm2(NUmodes, NSUPmodes, Nnutmodes);
        S_matrix = temperatureTurbulenceTerm(NTmodes, Nnutmodes);
    }
    B_total_matrix = B_matrix + BT_matrix;
    C_total_matrix.setSize(CT1_matrix.size());
    Info << "Summing the C tensors" << endl;
    for (label i = 0; i < C_total_matrix.size(); i++)
    {
        C_total_matrix[i] = CT1_matrix[i] + CT2_matrix[i];
    }
    // Eddy viscosity interpolation
    Info << "Starting the RBF interpolation process" << endl;
    Info << "The eddy viscosity field has Nnutmodes =" << Nnutmodes << " modes." << endl;
    Info << "The Nutfield variable has size: " << Nutfield.size() << endl;
    Info << "The nuTmodes variable has size: " << nuTmodes.size() << endl;
    Eigen::MatrixXd Ncoeff = ITHACAutilities::getCoeffs(Nutfield, nuTmodes);
    ITHACAstream::exportMatrix(Ncoeff, "Ncoeff", "python",
        "./ITHACAoutput/Matrices/");
    SAMPLES.resize(Nnutmodes);
    rbfsplines.resize(Nnutmodes);
    for (label i = 0; i < Nnutmodes; i++)
    {
        SAMPLES[i] = new SPLINTER::DataTable(1, 1);
        for (label j = 0; j < Ncoeff.cols(); j++)
        {
            SAMPLES[i]->addSample(mu.row(j), Ncoeff(i, j));
        }
        rbfsplines[i] = new SPLINTER::RBFSpline(*SAMPLES[i],
            SPLINTER::RadialBasisFunctionType::GAUSSIAN);
        Info << "Constructed RBF for mode " << i << endl;
    }
}

// * * * * * * * * * * * * * * Matrices Methods * * * * * * * * * * * * * * //

Eigen::MatrixXd UnsteadyBBTurb::BTturbulence(label NUmodes, label NSUPmodes)
{
    label BTsize = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd BT_matrix(BTsize, BTsize);
    BT_matrix = BT_matrix * 0;
    // Create PTRLIST with lift, velocities and supremizers
    PtrList<volVectorField> Together(0);

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k].clone());
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k].clone());
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k].clone());
        }
    }

    // Project everything
    for (label i = 0; i < BTsize; i++)
    {
        for (label j = 0; j < BTsize; j++)
        {
            BT_matrix(i, j) = fvc::domainIntegrate(Together[i] & (fvc::div(dev((T(fvc::grad(Together[j]))))))).value();
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "python",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "matlab",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(BT_matrix, "BT_matrix", "eigen",
        "./ITHACAoutput/Matrices/");
    return BT_matrix;
}

List<Eigen::MatrixXd> UnsteadyBBTurb::turbulenceTerm1(label NUmodes,
    label NSUPmodes, label Nnutmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    Info << "Creating CT matrix with: " << NUmodes << " velocity modes, "
         << NSUPmodes << " supremizer modes and " << Nnutmodes
         << " nut modes." << endl;
    Info << "Total size: " << Csize << "x" << Csize << endl;
    List<Eigen::MatrixXd> CT1_matrix;
    CT1_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        CT1_matrix[j].resize(Nnutmodes, Csize);
        CT1_matrix[j] = CT1_matrix[j] * 0;
    }

    PtrList<volVectorField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k].clone());
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k].clone());
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k].clone());
        }
    }
    for (label i = 0; i < Csize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix CT1_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                CT1_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & fvc::laplacian(nuTmodes[j], Together[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "python",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "matlab",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT1_matrix, "CT1_matrix", "eigen",
        "./ITHACAoutput/Matrices/CT1");
    return CT1_matrix;
}

List<Eigen::MatrixXd> UnsteadyBBTurb::turbulenceTerm2(label NUmodes,
    label NSUPmodes, label Nnutmodes)
{
    label Csize = NUmodes + NSUPmodes + liftfield.size();
    List<Eigen::MatrixXd> CT2_matrix;
    CT2_matrix.setSize(Csize);

    for (label j = 0; j < Csize; j++)
    {
        CT2_matrix[j].resize(Nnutmodes, Csize);
        CT2_matrix[j] = CT2_matrix[j] * 0;
    }

    PtrList<volVectorField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            Together.append(liftfield[k].clone());
        }
    }

    if (NUmodes != 0)
    {
        for (label k = 0; k < NUmodes; k++)
        {
            Together.append(Umodes[k].clone());
        }
    }

    if (NSUPmodes != 0)
    {
        for (label k = 0; k < NSUPmodes; k++)
        {
            Together.append(supmodes[k].clone());
        }
    }

    for (label i = 0; i < Csize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix CT2_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Csize; k++)
            {
                CT2_matrix[i](j, k) = fvc::domainIntegrate(Together[i] & (fvc::div(nuTmodes[j] * dev((fvc::grad(Together[k]))().T())))).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "python",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "matlab",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(CT2_matrix, "CT2_matrix", "eigen",
        "./ITHACAoutput/Matrices/CT2");
    return CT2_matrix;
}

List<Eigen::MatrixXd> UnsteadyBBTurb::temperatureTurbulenceTerm(
    label NTmodes, label Nnutmodes)
{
    label Stsize = NTmodes + liftfieldT.size();
    List<Eigen::MatrixXd> S_matrix;
    S_matrix.setSize(Stsize);
    Info << "Starting the creation of the S matrix with the following settings:" << endl;
    Info << "NTmodes = " << NTmodes << endl;
    Info << "Nnutmodes = " << Nnutmodes << endl;
    Info << "Stsize = " << Stsize << endl;

    for (label j = 0; j < Stsize; j++)
    {
        S_matrix[j].resize(Nnutmodes, Stsize);
    }

    PtrList<volScalarField> Together(0);

    // Create PTRLIST with lift, velocities and supremizers

    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            Together.append(liftfieldT[k].clone());
        }
    }

    if (NTmodes != 0)
    {
        for (label k = 0; k < NTmodes; k++)
        {
            Together.append(Tmodes[k].clone());
        }
    }
    Info << "Total size of the Together PTRLIST: " << Together.size() << endl;
    Info << "Total size of the nuTmodes PTRLIST: " << nuTmodes.size() << endl;
    for (label i = 0; i < Stsize; i++)
    {
        Info << "Filling layer number " << i + 1 << " in the matrix S_matrix" << endl;

        for (label j = 0; j < Nnutmodes; j++)
        {
            for (label k = 0; k < Stsize; k++)
            {
                S_matrix[i](j, k) = fvc::domainIntegrate(Together[i] * fvc::laplacian(nuTmodes[j], Together[k])).value();
            }
        }
    }

    // Export the matrix
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "python",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "matlab",
        "./ITHACAoutput/Matrices/");
    ITHACAstream::exportMatrix(S_matrix, "S_matrix", "eigen",
        "./ITHACAoutput/Matrices/S");
    return S_matrix;
}
