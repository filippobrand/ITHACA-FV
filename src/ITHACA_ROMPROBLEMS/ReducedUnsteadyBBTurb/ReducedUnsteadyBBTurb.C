/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2026 by the ITHACA-FV authors
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
/// Source file of the ReducedUnsteadyBBTurb class

#include "ReducedUnsteadyBBTurb.H"
#include "ReducedUnsteadyBBTurbSystems.H"
#include <atomic>
#include <fstream>
#include <chrono>

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedUnsteadyBBTurb::ReducedUnsteadyBBTurb()
{
}

ReducedUnsteadyBBTurb::ReducedUnsteadyBBTurb(UnsteadyBBTurb& fom):
    problem(&fom)
{
    // Construct the settings objects for the ROM
    romSettings_ = ROMSettings
    {
        problem->ITHACAdict->lookupOrDefault<word>("method", "PPE"),
        problem->ITHACAdict->lookupOrDefault<word>("bcMethod", "Gunzburger"),
        problem->ITHACAdict->lookupOrDefault<bool>("hasMonitors", false),
        problem->ITHACAdict->lookupOrDefault<bool>("timeDependentBC", false),
        (int)problem->inletIndex.rows(),
        (int)problem->inletIndexT.rows(),
        problem->ITHACAdict->lookupOrDefault<word>("solverODE", "BDF2")
    };
    M_Assert(romSettings_.bcMethod == "Gunzburger",
             "Only Gunzburger BC method is implemented for now. Please set bcMethod to Gunzburger in the ITHACA dictionary.");
    collectMatrices();
    initializeDimensions();

    if (romSettings_.bcMethod == "Gunzburger")
    {
        odeSystem_ = std::make_unique<SystemPPEGunzburger>(*this);
    }
    else if (romSettings_.bcMethod == "penalty")
    {
        penaltySettings_ = PenaltySettings
        {
            problem->ITHACAdict->lookupOrDefault<bool>("optimizeTau", false),
            problem->ITHACAdict->lookupOrDefault<int>("nTimestepsPenaltyAdapt", 100),
            problem->ITHACAdict->lookupOrDefault<int>("penaltyMaxIter", 50),
            problem->ITHACAdict->lookupOrDefault<float>("penaltyTolU", 1e-3),
            problem->ITHACAdict->lookupOrDefault<float>("penaltyTolT", 1e-2)
        };
        M_Assert(penaltySettings_.nTimestepsPenaltyAdapt > 0,
                 "Number of timesteps for penalty factor adaptation must be positive.");
        M_Assert(penaltySettings_.penaltyMaxIter > 0,
                 "Maximum number of iterations for penalty factor optimization must be positive.");
        M_Assert(penaltySettings_.penaltyTolU > 0 and penaltySettings_.penaltyTolT > 0,
                 "Penalty factor optimization tolerances must be positive.");
    }

    M_Assert(odeSystem_ != nullptr,
             "ODE system is not initialized. Please check the bcMethod in the ITHACA dictionary.");
    y = Eigen::VectorXd::Zero(expressionModes_.velocity + expressionModes_.pressure
                              + expressionModes_.temperature);
    interpolationSettings_ = InterpolationSettings
    {
        problem->ITHACAdict->lookupOrDefault<int>("firstRBFIndex", 0),
        problem->dimA,
        problem->ITHACAdict->lookupOrDefault<bool>("derivativeInRBF", false)
    };

    if (romSettings_.solverODE == "BDF2")
    {
        // The solver takes in a child of ReducedODESystem, which it whathever odeSystem_ points to
        odeSolver_ = std::make_unique<BDF2Solver>(*odeSystem_, y.size());
    }
    else if (romSettings_.solverODE == "BDF1"
             or romSettings_.solverODE == "BackwardEuler")
    {
        odeSolver_ = std::make_unique<ImplicitEulerSolver>(*odeSystem_, y.size());
    }
    else
    {
        Info << "ODE solver not recognized. Using BDF2 as default." << endl;
        odeSolver_ = std::make_unique<BDF2Solver>(*odeSystem_, y.size());
    }

    // Note: some other reduced class in ITHACA-FV in the constructor create a local copy of the
    // modes from the FOM problem. Maybe in the future the same could be done here.
}

ReducedUnsteadyBBTurb::~ReducedUnsteadyBBTurb() = default;

void ReducedUnsteadyBBTurb::collectMatrices()
{
    pCommonMatrices = problem->pCommonMatrices;

    if (romSettings_.method == "supremizer")
    {
        pSupremizerMatrices = problem->pSupremizerMatrices;
        pPPEMatrices = nullptr;
    }
    else if (romSettings_.method == "PPE")
    {
        pPPEMatrices = problem->pPPEMatrices;
        pSupremizerMatrices = nullptr;
    }

    if (romSettings_.bcMethod == "penalty")
    {
        pPenaltyMatrices = problem->pPenaltyMatrices;
        pGunzburgerMatrices = nullptr;
    }
    else if (romSettings_.bcMethod == "Gunzburger")
    {
        pGunzburgerMatrices = problem->pGunzburgerMatrices;
        pPenaltyMatrices = nullptr;
    }
    else
    {
        pPenaltyMatrices = nullptr;
        pGunzburgerMatrices = nullptr;
    }
}

void ReducedUnsteadyBBTurb::initializeDimensions()
{
    if (romSettings_.bcMethod == "Gunzburger")
    {
        // Here the test basis and the expression basis are different, so we need to differentiate
        testModes_ = NumberOfModes
        {
            pCommonMatrices->M.rows(), // Mass matrix <phi_i, phi_j> for velocity
            pCommonMatrices->K.cols(), // Should maybe implement a PPE based computation
            pCommonMatrices->W.rows(), // Temperature Mass matrix <xi_i, xi_j>
            pCommonMatrices->CTotal.dimension(1)
        };
        expressionModes_ = NumberOfModes
        {
            pCommonMatrices->M.cols(), // Mass matrix <phi_i, phi_j> for velocity
            pCommonMatrices->K.cols(),
            pCommonMatrices->W.cols(), // Temperature mass matrix <xi_i, xi_j>
            pCommonMatrices->CTotal.dimension(1)
        };
        // Check whether any of the number of modes is invalid (negative or zero)
        M_Assert(testModes_.velocity > 0,
                 "Number of test modes for velocity is invalid.");
        M_Assert(testModes_.pressure > 0,
                 "Number of test modes for pressure is invalid.");
        M_Assert(testModes_.temperature > 0,
                 "Number of test modes for temperature is invalid.");
        M_Assert(testModes_.nut > 0, "Number of test modes for nut is invalid.");
        M_Assert(expressionModes_.velocity > 0,
                 "Number of expression modes for velocity is invalid.");
        M_Assert(expressionModes_.pressure > 0,
                 "Number of expression modes for pressure is invalid.");
        M_Assert(expressionModes_.temperature > 0,
                 "Number of expression modes for temperature is invalid.");
        M_Assert(expressionModes_.nut > 0,
                 "Number of expression modes for nut is invalid.");
        Info << "### DEBUG --- Number of test modes: " << testModes_.velocity <<
             " velocity, "
             << testModes_.pressure << " pressure, "
             << testModes_.temperature << " temperature, "
             << testModes_.nut << " nut." << endl;
        Info << "### DEBUG --- Number of expression modes: " <<
             expressionModes_.velocity << " velocity, "
             << expressionModes_.pressure << " pressure, "
             << expressionModes_.temperature << " temperature, "
             << expressionModes_.nut << " nut." << endl;
    }
    else
    {
        Info << "Not implemented anything other than Gunzburger BC method for now." <<
             endl;
    }
}

void ReducedUnsteadyBBTurb::interpolateNutCoeffs()
{
    Eigen::VectorXd velocity_coeffs = y.head(expressionModes_.velocity);

    if (interpolationSettings_.derivative)
    {
        // We use the OdeSolver own method to compute the derivative
        Eigen::VectorXd velocity_derivative;
        velocity_derivative = odeSolver_->getLatestDerivative(velocity_coeffs);
        Eigen::VectorXd inputRBF(velocity_coeffs.size() + velocity_derivative.size());
        inputRBF << velocity_coeffs, velocity_derivative;

        for (int j = 0; j < expressionModes_.nut; j++)
        {
            nutCurrentCoeffs_(j) = problem->rbfSplines[j]->predict(inputRBF);
        }
    }
    else
    {
        for (int j = 0; j < expressionModes_.nut; j++)
        {
            nutCurrentCoeffs_(j) = problem->rbfSplines[j]->predict(velocity_coeffs);
        }
    }
}

void ReducedUnsteadyBBTurb::readEigenvalues()
{
    // We read the eigenvalues from the files in the ./ITHACAoutput/POD folder
    // These files follow this form:
    // %%MatrixMarket matrix array real general
    // 2720 1
    // 0.34428033233490540344
    // 0.12360804880946965612
    // 0.08097235508005383442
    // ...
    //
    // We ignore the first two lines and read the first n lines, where n is the number of modes for each variable.
    // Velocity - file named: Eigenvalues_U
    std::ifstream uFile("ITHACAoutput/POD/Eigenvalues_U");
    M_Assert(uFile.is_open(),
             "Could not open file ITHACAoutput/POD/Eigenvalues_U. Please make sure the file exists and is readable.");
    std::string line;
    std::getline(uFile, line); // Ignore first line
    std::getline(uFile, line); // Ignore second line
    uEigenvalues_.resize(expressionModes_.velocity);

    for (int i = 0; i < expressionModes_.velocity; i++)
    {
        std::getline(uFile, line);
        uEigenvalues_(i) = std::stod(line);
    }

    // Pressure - file named: Eigenvalues_p_rgh
    std::ifstream pFile("ITHACAoutput/POD/Eigenvalues_p_rgh");
    M_Assert(pFile.is_open(),
             "Could not open file ITHACAoutput/POD/Eigenvalues_p_rgh. Please make sure the file exists and is readable.");
    std::getline(pFile, line); // Ignore first line
    std::getline(pFile, line); // Ignore second line
    pEigenvalues_.resize(expressionModes_.pressure);

    for (int i = 0; i < expressionModes_.pressure; i++)
    {
        std::getline(pFile, line);
        pEigenvalues_(i) = std::stod(line);
    }

    // Temperature - file named: Eigenvalues_T
    std::ifstream tFile("ITHACAoutput/POD/Eigenvalues_T");
    M_Assert(tFile.is_open(),
             "Could not open file ITHACAoutput/POD/Eigenvalues_T. Please make sure the file exists and is readable.");
    std::getline(tFile, line); // Ignore first line
    std::getline(tFile, line); // Ignore second line
    tEigenvalues_.resize(expressionModes_.temperature);

    for (int i = 0; i < expressionModes_.temperature; i++)
    {
        std::getline(tFile, line);
        tEigenvalues_(i) = std::stod(line);
    }

    // Fluctnut - file named: Eigenvalues_fluctNut
    std::ifstream nutFile("ITHACAoutput/POD/Eigenvalues_fluctNut");
    M_Assert(nutFile.is_open(),
             "Could not open file ITHACAoutput/POD/Eigenvalues_fluctNut. Please make sure the file exists and is readable.");
    std::getline(nutFile, line); // Ignore first line
    std::getline(nutFile, line); // Ignore second line
    nutEigenvalues_.resize(expressionModes_.nut);

    for (int i = 0; i < expressionModes_.nut; i++)
    {
        std::getline(nutFile, line);
        nutEigenvalues_(i) = std::stod(line);
    }

    Info << "### EIGS - Velocity eigenvalues: " << uEigenvalues_.transpose() <<
         endl;
    Info << "### EIGS - Pressure eigenvalues: " << pEigenvalues_.transpose() <<
         endl;
    Info << "### EIGS - Temperature eigenvalues: " << tEigenvalues_.transpose() <<
         endl;
    Info << "### EIGS - FluctNut eigenvalues: " << nutEigenvalues_.transpose() <<
         endl;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyBBTurb::solveOnline(
    const Eigen::MatrixXd& vel_now_BC,
    const Eigen::MatrixXd& temp_now_BC,
    int startSnap, TimeManager& timeManager)
{
    boundaryConditions_ = BoundaryConditions(vel_now_BC, temp_now_BC, "linear");
    boundaryConditions_.initializeReducedCoeffs(
        startSnap, y, problem,
        expressionModes_.velocity, expressionModes_.pressure,
        expressionModes_.temperature
    );
    nutCurrentCoeffs_ = ITHACAutilities::getCoeffs(
                            problem->fluctNutfield[startSnap], problem->nutmodes);
    Info << "### DEBUG --- Initialized nut coefficients" << endl;
    // Temporary: for our test case, only the first two values in boundaryConditions_.getCurrentBCs() are used
    // in the interpolation.
    nutAvgCurrentCoeffs_ = interpolateIDW(boundaryConditions_.getCurrentBCs().head(
            2)); // Interpolate current BCs (parameters)
    onlineSolution.resize(timeManager.getNumberOfStepsToSave() + 1, y.size() + 1);
    nutCoeffMat.resize(timeManager.getNumberOfStepsToSave() + 1,
                       expressionModes_.nut + 1);
    onlineSolution.setZero();
    nutCoeffMat.setZero();
    onlineSolution(0, 0) = timeManager.getCurrentTime();
    onlineSolution.block(0, 1, 1, y.size()) = y.transpose();
    nutCoeffMat(0, 0) = timeManager.getCurrentTime();
    nutCoeffMat.block(0, 1, 1,
                      expressionModes_.nut) = nutCurrentCoeffs_.transpose();
    int save_index = 1;

    while (!timeManager.isFinished())
    {
        timeManager.advanceTime();
        boundaryConditions_.updateTimeDependentBC(timeManager.getCurrentTime());
        odeSolver_->solveStep(y, timeManager.getCurrentTime(),
                              timeManager.getActualTimeStep(), false);
        // ADD THE INTERPOLATION OF THE NUT COEFFICIENTS AVG HERE
        interpolateNutCoeffs();

        if (timeManager.shouldSaveCoefficients())
        {
            Info << "### DEBUG --- Saving coefficients at time " <<
                 timeManager.getCurrentTime() << endl;
            onlineSolution(save_index, 0) = timeManager.getCurrentTime();
            onlineSolution.block(save_index, 1, 1, y.size()) = y.transpose();
            nutCoeffMat(save_index, 0) = timeManager.getCurrentTime();
            nutCoeffMat.block(save_index, 1, 1,
                              expressionModes_.nut) = nutCurrentCoeffs_.transpose();
            save_index++;
            Info << "### DEBUG --- Saved coefficients at time " <<
                 timeManager.getCurrentTime() << endl;
            timeManager.updateAfterSave();
        }
    }

    if (save_index < onlineSolution.rows())
    {
        onlineSolution.conservativeResize(save_index, y.size() + 1);
        nutCoeffMat.conservativeResize(save_index, expressionModes_.nut + 1);
    }
}

void ReducedUnsteadyBBTurb::reconstructSolution(TimeManager& time_manager,
        bool exportFields, fileName folder)
{
    if (exportFields)
    {
        if (Pstream::master())
        {
            mkDir(folder);
            ITHACAutilities::createSymLink(folder);
        }
    }

    int exportEveryIndex = time_manager.getExportToEverySaved();

    if (exportEveryIndex < 1)
    {
        Info << "Export every saved index is less than 1. Setting it to 1." << endl;
        exportEveryIndex = 1;
    }

    List<Eigen::MatrixXd> CoeffU;
    List<Eigen::MatrixXd> CoeffPrgh;
    List<Eigen::MatrixXd> CoeffT;
    List<Eigen::MatrixXd> CoeffNut;
    // Since we know the size, preallocate the lists:
    CoeffU.resize(time_manager.getNumberOfStepsToExport() + 1);
    CoeffPrgh.resize(time_manager.getNumberOfStepsToExport() + 1);
    CoeffT.resize(time_manager.getNumberOfStepsToExport() + 1);
    CoeffNut.resize(time_manager.getNumberOfStepsToExport() + 1);

    for (int i = 0; i < onlineSolution.rows(); i++)
    {
        if (i % exportEveryIndex == 0)
        {
            Eigen::MatrixXd currentUCoeff = onlineSolution.block(i, 1, 1,
                expressionModes_.velocity);
            Eigen::MatrixXd currentPrghCoeff = onlineSolution.block(i,
                1 + expressionModes_.velocity, 1, expressionModes_.pressure);
            Eigen::MatrixXd currentTCoeff = onlineSolution.block(i,
                1 + expressionModes_.velocity + expressionModes_.pressure, 1,
                expressionModes_.temperature);
            Eigen::MatrixXd currentNutCoeff = nutCoeffMat.block(i, 1, 1,
                expressionModes_.nut);
            CoeffU[i / exportEveryIndex] = currentUCoeff;
            CoeffPrgh[i / exportEveryIndex] = currentPrghCoeff;
            CoeffT[i / exportEveryIndex] = currentTCoeff;
            CoeffNut[i / exportEveryIndex] = currentNutCoeff;
        }
    }

    volVectorField uRec("uRec", problem->L_U_SUPmodes[0]);
    volScalarField TRec("TRec", problem->L_Tmodes[0]);
    volScalarField prghRec("prghRec", problem->P_rghmodes[0]);
    volScalarField nutFluctRec("nutFluctRec", problem->nutmodes[0]);
    uRecFields = problem->L_U_SUPmodes.reconstruct(uRec, CoeffU, "uRec");
    TRecFields = problem->L_Tmodes.reconstruct(TRec, CoeffT, "TRec");
    nutFluctRecFields = problem->nutmodes.reconstruct(nutFluctRec, CoeffNut,
        "nutFluctRec");
    prghRecFields = problem->P_rghmodes.reconstruct(prghRec, CoeffPrgh, "prghRec");
    // Reconstruct the averaged eddy viscosity field as a linear combination
    volScalarField nutAvg(
        IOobject(
            "nutAvgRec",
            problem->nutmodes[0].time().timeName(),
            problem->nutmodes[0].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        problem->nutmodes[0].mesh(),
        dimensionedScalar("zero", problem->nutmodes[0].dimensions(), 0.0));

    for (int k = 0; k < nutAvgCurrentCoeffs_.size(); k++)
    {
        nutAvg += nutAvgCurrentCoeffs_(k) * problem->avgNutfield[k];
    }

    nutRecFields.resize(0);
    forAll(nutFluctRecFields, i)
    {
        nutRecFields.append(nutFluctRecFields[i] + nutAvg);
    }

    if (exportFields)
    {
        ITHACAstream::exportFields(uRecFields, folder, "uRec");
        ITHACAstream::exportFields(TRecFields, folder, "TRec");
        ITHACAstream::exportFields(nutFluctRecFields, folder, "nutFluctRec");
        ITHACAstream::exportFields(prghRecFields, folder, "prghRec");
        ITHACAstream::exportFields(nutRecFields, folder, "nutRec");
    }
}

void ReducedUnsteadyBBTurb::saveCoefficients(word folder, word modeIdentifier)
{
    if (Pstream::master())
    {
        mkDir("./ITHACAoutput/ReducedCoefficients/");
        word true_folder = "./ITHACAoutput/ReducedCoefficients/" + folder + "_" +
                           modeIdentifier + "/";
        ITHACAstream::exportMatrix(online_solution, "reducedCoefficients", "python",
                                   true_folder);
        ITHACAstream::exportMatrix(nutCoeffMat, "RBFCoefficients", "python",
                                   true_folder);
    }
}

// * * * * * * * * *  Inverse Distance Weighting Functions  * * * * * * * * //
Eigen::VectorXd ReducedUnsteadyBBTurb::interpolateIDW(const Eigen::VectorXd&
        input_parameters)
{
    // The MatrixXd of offline parameters is stored in problem->mu.
    label nOfflineSamples = problem->mu.cols();
    Eigen::VectorXd interpolatedNutCoeffs(nOfflineSamples);
    Eigen::VectorXd weights(nOfflineSamples);

    for (label i = 0; i < nOfflineSamples; i++)
    {
        weights(i) = 1.0 / ((input_parameters - problem->mu.col(
                                 i)).norm() + 1e-10); // Add a small value to avoid division by zero
    }

    double weightSum = weights.sum();

    for (label j = 0; j < nOfflineSamples; j++)
    {
        Eigen::VectorXd nutCoeffsAtJ(nOfflineSamples);
        nutCoeffsAtJ.setZero();
        nutCoeffsAtJ[j] = 1.0;
        interpolatedNutCoeffs(j) = weights.dot(nutCoeffsAtJ) / weightSum;
    }

    Info << "The interpolated eddy viscosity coefficients are: " <<
         interpolatedNutCoeffs << endl;
    return interpolatedNutCoeffs;
}

// * * * * * * * * *  Validation and setup helpers  * * * * * * * * //

void ReducedUnsteadyBBTurb::inf_sup_constant()
{
    double a;
    Eigen::VectorXd sup(expressionModes_.velocity);
    Eigen::VectorXd inf(expressionModes_.pressure);

    for (int i = 0; i < expressionModes_.pressure; i++)
    {
        for (int j = 0; j < expressionModes_.velocity; j++)
        {
            sup(j) = fvc::domainIntegrate(fvc::div(problem->L_U_SUPmodes[j]) *
                                          problem->P_rghmodes[i]).value() /
                     ITHACAutilities::H1Seminorm(problem->L_U_SUPmodes[j]) / ITHACAutilities::L2Norm(
                         problem->P_rghmodes[i]);
        }

        inf(i) = sup.maxCoeff();
    }

    a = inf.minCoeff();
    Info << "### STABILITY: The inf-sup constant is: " << a << endl;
}

void ReducedUnsteadyBBTurb::setFluidProperties(scalar nu, scalar Pr,
        scalar Pr_t)
{
    fluidProperties.nu = nu;
    fluidProperties.Pr = Pr;
    fluidProperties.Pr_t = Pr_t;
}

void ReducedUnsteadyBBTurb::reset()
{
    online_solution.clear();
    nutCoeffMat.setZero();
    uRecFields.clear();
    TRecFields.clear();
    nutFluctRecFields.clear();
    nutRecFields.clear();
}
// ************************************************************************ //
