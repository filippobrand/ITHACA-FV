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
#include <atomic>
#include <fstream>
#include <chrono>


// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedUnsteadyBBTurb::ReducedUnsteadyBBTurb()
{
}

ReducedUnsteadyBBTurb::ReducedUnsteadyBBTurb(UnsteadyBBTurb& FOMproblem):
    problem(&FOMproblem)
{
    pCommonMatrices = FOMproblem.pCommonMatrices;
    pSupremizerMatrices = FOMproblem.pSupremizerMatrices;
    pPPEMatrices = FOMproblem.pPPEMatrices;
    pPenaltyMatrices = FOMproblem.pPenaltyMatrices;

    N_BC = problem->inletIndex.rows();
    N_BC_t = problem->inletIndexT.rows();
    Nphi_u = pCommonMatrices->B.cols();
    Ntest_u = pCommonMatrices->B.rows();
    Nphi_prgh = pCommonMatrices->K.cols();
    Ntest_prgh = pCommonMatrices->K.cols();
    Nphi_t = pCommonMatrices->Y.cols();
    Ntest_t = pCommonMatrices->Y.rows();
    Nphi_nut = pCommonMatrices->CTotal.dimension(1);

    Info << "Reading the matrices, the following number of modes will be used: " << endl;
    Info << "Velocity modes and test functions: " << Nphi_u << " and " << Ntest_u << endl;
    Info << "Pressure modes and test functions: " << Nphi_prgh << " and " << Ntest_prgh << endl;
    Info << "Temperature modes and test functions: " << Nphi_t << " and " << Ntest_t << endl;
    Info << "Eddy viscosity modes: " << Nphi_nut << endl;

    dimA = problem->dimA;
    method = problem->ITHACAdict->lookupOrDefault<word>("method", "supremizer");

    y = Eigen::VectorXd::Zero(Nphi_u + Nphi_prgh + Nphi_t); // Solution vector

    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFIndex = problem->skipRBFIndex;
        Info << "Skipping lifting modes in the RBF evaluation. This means that the first RBF index is set to " << firstRBFIndex << endl;
    } else
    {
        firstRBFIndex = 0;
        Info << "## COMM - Not skipping lifting modes in the RBF evaluation. This means that the first RBF index is set to " << firstRBFIndex << endl;
    }

    // Penalty factor optimization settings
    optimizePenaltyFactor = problem->ITHACAdict->lookupOrDefault<bool>("optimizeTau", false);
    nTimestepsPenaltyAdapt = problem->ITHACAdict->lookupOrDefault<int>("nTimestepsPenaltyAdapt", 100);
    penaltyMaxIter = problem->ITHACAdict->lookupOrDefault<int>("penaltyMaxIter", 50);
    penaltyTolU = problem->ITHACAdict->lookupOrDefault<float>("penaltyTolU", 1e-3);
    penaltyTolT = problem->ITHACAdict->lookupOrDefault<float>("penaltyTolT", 1e-2);

    M_Assert(penaltyTolU > 0 && penaltyTolT > 0, "Penalty factor optimization tolerances must be positive.");
    inf_sup_constant();
    // Note: some other reduced class in ITHACA-FV in the constructor create a local copy of the
    // modes from the FOM problem. Maybe in the future the same could be done here
}


// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyBBTurb::solveOnline_PPE()
{
    // Set number of online solutions - Time related
    int numberOfStores = round((storeEvery) / dt); // Number of time steps between two stored solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores); // Total stored solutions, excluding initial condition
    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(Nphi_nut + 1, onlineSize);
    time = tstart;

    ODEStructurePPETurb ode_structure(*this);
    BDF2Solver solver(ode_structure, y.size());

    onlineTimeLoop(solver, firstRBFIndex, numberOfStores);
}

template <typename ODESolver>
void ReducedUnsteadyBBTurb::interpolateEddyViscosity(ODESolver& ode_solver)
{
  Eigen::VectorXd RBFInput = Eigen::VectorXd::Zero(dimA);
  if (problem->derivativeInRBF == true)
  {
    Eigen::VectorXd aDer = Eigen::VectorXd::Zero(Nphi_u);
    aDer = (y.head(Nphi_u) - ode_solver.y_old.head(Nphi_u)) / dt;
    RBFInput << y.segment(firstRBFIndex, dimA / 2), aDer.segment(firstRBFIndex, dimA / 2);
  } else
  {
      RBFInput << y.segment(firstRBFIndex, dimA);
  }
  for (int j = 0; j < Nphi_nut; j++)
  {
      nut_fluct(j) = problem->rbfSplines[j]->predict(RBFInput);
  }
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyBBTurb::reconstructSolution(bool exportFields, fileName folder)
{
    if (exportFields)
    {
        if (Pstream::master())
        {
            mkDir(folder);
            ITHACAutilities::createSymLink(folder);
        }
    }

    int timeStepCounter = 0;
    int nextWrite = 0;
    int exportEveryIndex = round(exportEvery / storeEvery);
    List<Eigen::MatrixXd> CoeffU;
    List<Eigen::MatrixXd> CoeffPrgh;
    List<Eigen::MatrixXd> CoeffT;
    List<Eigen::MatrixXd> CoeffNut;
    CoeffU.resize(0);
    CoeffPrgh.resize(0);
    CoeffT.resize(0);
    CoeffNut.resize(0);

    int reconstructionSizeU = Nphi_u;
    int reconstructionSizeT = Nphi_t;

    for (int i = 0; i < online_solution.size(); i++)
    {
        if (timeStepCounter == nextWrite)
        {
            Eigen::MatrixXd currentUCoeff;
            Eigen::MatrixXd currentPrghCoeff;
            Eigen::MatrixXd currentTCoeff;
            Eigen::MatrixXd currentNutCoeff;

            currentUCoeff = online_solution[i].block(1, 0, reconstructionSizeU, 1);
            if (method != "VMB")
            {
                currentPrghCoeff = online_solution[i].block(reconstructionSizeU + 1, 0, Nphi_prgh, 1);
            } else
            {
                currentPrghCoeff = online_solution[i].block(1, 0, reconstructionSizeU, 1);
            }
            currentTCoeff = online_solution[i].bottomRows(reconstructionSizeT);
            currentNutCoeff = rbfCoeffMat.block(1, i, Nphi_nut, 1);

            CoeffPrgh.append(currentPrghCoeff);
            CoeffU.append(currentUCoeff);
            CoeffT.append(currentTCoeff);
            CoeffNut.append(currentNutCoeff);

            nextWrite += exportEveryIndex;
        }
        timeStepCounter++;
    }
    volVectorField uRec("uRec", problem->L_U_SUPmodes[0]);
    volScalarField TRec("TRec", problem->L_Tmodes[0]);
    volScalarField prghRec("prghRec", problem->P_rghmodes[0]);
    volScalarField nutFluctRec("nutFluctRec", problem->nutmodes[0]);

    uRecFields = problem->L_U_SUPmodes.reconstruct(uRec, CoeffU, "uRec");
    TRecFields = problem->L_Tmodes.reconstruct(TRec, CoeffT, "TRec");
    nutFluctRecFields = problem->nutmodes.reconstruct(nutFluctRec, CoeffNut, "nutFluctRec");
    prghRecFields = problem->P_rghmodes.reconstruct(prghRec, CoeffPrgh, "prghRec");

    // Reconstruct the averaged eddy viscosity field as a linear combination of the nut_param coefficients and the PtrList<volScalarField> avgNutfield;
    volScalarField nutAvg(
        IOobject(
            "nutAvgRec",
            problem->nutmodes[0].time().timeName(),
            problem->nutmodes[0].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        problem->nutmodes[0].mesh(),
        dimensionedScalar("zero", problem->nutmodes[0].dimensions(), 0.0));

    for (int k = 0; k < nut_param.size(); k++)
    {
        nutAvg += nut_param(k) * problem->avgNutfield[k];
    }

    volScalarField nutRecField(
        IOobject(
            "nutRec",
            problem->nutmodes[0].time().timeName(),
            problem->nutmodes[0].mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        problem->nutmodes[0].mesh(),
        dimensionedScalar("zero", problem->nutmodes[0].dimensions(), 0.0));

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

void ReducedUnsteadyBBTurb::saveCoefficients(word folder)
{
    if (Pstream::master())
    {
        mkDir("./ITHACAoutput/ReducedCoefficients/");
        ITHACAstream::exportMatrix(online_solution, "reducedCoefficients", "python", "./ITHACAoutput/ReducedCoefficients/" + folder + "/");
        ITHACAstream::exportMatrix(rbfCoeffMat, "RBFCoefficients", "python", "./ITHACAoutput/ReducedCoefficients/" + folder + "/");
    }
}
// * * * * * * * * *  Inverse Distance Weighting Functions  * * * * * * * * //
Eigen::VectorXd ReducedUnsteadyBBTurb::interpolateIDW()
{
    // The MatrixXd of offline parameters is stored in problem->mu.
    label nOfflineSamples = problem->mu.cols();
    Eigen::VectorXd interpolatedNutCoeffs(nOfflineSamples);
    Eigen::VectorXd weights(nOfflineSamples);
    Info << "The shape of problem->mu is: " << problem->mu.rows() << " rows x " << problem->mu.cols() << " columns" << endl;
    Info << "mu_now has shape: " << mu_now.rows() << " rows x " << mu_now.cols() << " columns" << endl;
    for (label i = 0; i < nOfflineSamples; i++)
    {
        weights(i) = 1.0 / ((mu_now - problem->mu.col(i)).norm() + 1e-10); // Add a small value to avoid division by zero
    }
    double weightSum = weights.sum();
    for (label j = 0; j < nOfflineSamples; j++)
    {
        Eigen::VectorXd nutCoeffsAtJ(nOfflineSamples);
        nutCoeffsAtJ.setZero();
        nutCoeffsAtJ[j] = 1.0;
        interpolatedNutCoeffs(j) = weights.dot(nutCoeffsAtJ) / weightSum;
    }
    Info << "The interpolated eddy viscosity coefficients are: " << interpolatedNutCoeffs << endl;
    return interpolatedNutCoeffs;
}

void ReducedUnsteadyBBTurb::setTimeSettings(const Eigen::MatrixXd& timeMatrix, int index)
{
    tstart = timeMatrix(index, 0);
    finalTime = timeMatrix(index, 1);
    dt = timeMatrix(index, 2);
    storeEvery = timeMatrix(index, 3);
    exportEvery = timeMatrix(index, 4);
}

// * * * * * * * * *  Validation and setup helpers  * * * * * * * * //
void ReducedUnsteadyBBTurb::validateSettings()
{
    M_Assert(exportEvery >= dt,
        "The time step dt must be smaller than exportEvery.");
    M_Assert(storeEvery >= dt,
        "The time step dt must be smaller than storeEvery.");
    M_Assert(ITHACAutilities::isInteger(storeEvery / dt) == true,
        "The variable storeEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / dt) == true,
        "The variable exportEvery must be an integer multiple of the time step dt.");
    M_Assert(ITHACAutilities::isInteger(exportEvery / storeEvery) == true,
        "The variable exportEvery must be an integer multiple of the variable storeEvery.");
}

template <typename ODESolver>
void ReducedUnsteadyBBTurb::onlineTimeLoop(
    ODESolver& ode_solver, int firstRBFIndex, int numberOfStores)
{
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    Eigen::VectorXd res(y); // Residual vector

    int timeStepCounter = 0;
    int nextStore = 0;
    int storedSnapshotsCounter = 0;

    Eigen::MatrixXd tmp_sol(y.rows() + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    while (time < finalTime - dt)
    {
        time += dt;
        boundaryConditions.updateTimeDependentBC(time);
        res.setZero();
        ode_solver.solveStep(y, time, dt);
        interpolateEddyViscosity(ode_solver);

        if (problem->bcMethod == "lift")
        {
            boundaryConditions.correctLiftingCoeffs(y, N_BC, N_BC_t, Nphi_u, Nphi_prgh);
        }

        tmp_sol(0) = time;
        tmp_sol.col(0).tail(y.rows()) = y;

        if (timeStepCounter == nextStore)
        {
            if (storedSnapshotsCounter >= online_solution.size())
            {
                online_solution.append(tmp_sol);
            } else
            {
                online_solution[storedSnapshotsCounter] = tmp_sol;
            }
            rbfCoeffMat(0, storedSnapshotsCounter) = time;
            rbfCoeffMat.block(1, storedSnapshotsCounter, Nphi_nut, 1) = nut_fluct;
            nextStore += numberOfStores;
            storedSnapshotsCounter++;
        }
        timeStepCounter++;
    }
    Info << "Online solution computed, with total time steps solved: " << timeStepCounter << endl;
}

void ReducedUnsteadyBBTurb::inf_sup_constant()
{
    double a;
    Eigen::VectorXd sup(Nphi_u);
    Eigen::VectorXd inf(Nphi_prgh);

    for (int i = 0; i < Nphi_prgh; i++)
    {
        for (int j = 0; j < Nphi_u; j++)
        {
            sup(j) = fvc::domainIntegrate(fvc::div(problem->L_U_SUPmodes[j]) * problem->P_rghmodes[i]).value() /
                ITHACAutilities::H1Seminorm(problem->L_U_SUPmodes[j]) / ITHACAutilities::L2Norm(problem->P_rghmodes[i]);
        }
        inf(i) = sup.maxCoeff();
    }
    a = inf.minCoeff();
    Info << "### STABILITY: The inf-sup constant is: " << a << endl;
}

void ReducedUnsteadyBBTurb::solveOnline(const Eigen::MatrixXd& vel_now_BC, const Eigen::MatrixXd& temp_now_BC, int startSnap)
{
    Info << "########### Online solve N° " << count_online_solve << " ###########" << endl;
    validateSettings();
    boundaryConditions = BoundaryConditions(vel_now_BC, temp_now_BC, "linear");
    boundaryConditions.initializeReducedCoeffs(startSnap, y, problem, Nphi_u, Nphi_prgh, Nphi_t, N_BC, N_BC_t);
    if (problem->bcMethod == "lift")
    {
        boundaryConditions.correctLiftingCoeffs(y, N_BC, N_BC_t, Nphi_u, Nphi_prgh);
    }

    nut_fluct = ITHACAutilities::getCoeffs(problem->fluctNutfield[startSnap], problem->nutmodes);
    nut_param = interpolateIDW();

    M_Assert(method == "PPE", "Currently, only the PPE stabilization method is implemented for the online solve. Please select 'PPE' as method in the ITHACAdict or implement the selected method.");

    if (problem->bcMethod == "penalty" && optimizePenaltyFactor == true)
    {
        Info << "Penalty factor estimation for PPE is not available yet." << endl;
    }
    boundaryConditions.initializeReducedCoeffs(startSnap, y, problem, Nphi_u, Nphi_prgh, Nphi_t, N_BC, N_BC_t);
    auto clockStart = std::chrono::steady_clock::now();
    solveOnline_PPE();
    auto duration = std::chrono::steady_clock::now() - clockStart;
    Info << "The online solution took, in milliseconds: " << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " ms" << endl;

    count_online_solve++;
}

void ReducedUnsteadyBBTurb::reset()
{
    online_solution.clear();
    rbfCoeffMat.setZero();

    uRecFields.clear();
    TRecFields.clear();
    nutFluctRecFields.clear();
    nutRecFields.clear();
}

// * * * * * * * * *  BBTurbODESystemPPE class implementation  * * * * * * * * //

void ODEStructurePPETurb::evaluateResidual(const Eigen::VectorXd& state, const Eigen::VectorXd& state_dot, Eigen::VectorXd& residual, double t) const
{
    // Here we evaluate the residual for the PPE formulation
    Eigen::VectorXd a = state.head(rom.Nphi_u);
    Eigen::VectorXd a_dot = state_dot.head(rom.Nphi_u);
    Eigen::VectorXd b = state.segment(rom.Nphi_u, rom.Nphi_prgh);
    Eigen::VectorXd c = state.tail(rom.Nphi_t);
    Eigen::VectorXd c_dot = state_dot.tail(rom.Nphi_t);

    // Convective term momentum equation
    Eigen::MatrixXd cc(1, 1);
    // Non linear term turbulence in momentum equation - average nut
    Eigen::MatrixXd caveraged(1, 1);
    // Non linear term turbulence in momentum equation - fluctuating nut
    Eigen::MatrixXd ct(1, 1);
    // Tensor term PPE equation
    Eigen::MatrixXd gg(1, 1);
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    // Turbulence term temperature - average part
    Eigen::MatrixXd qt_averaged(1, 1);
    // Turbulence term temperature - fluctuating part
    Eigen::MatrixXd qt(1, 1);
    // Turbulence term - PPE
    Eigen::MatrixXd turbPPE(1, 1);
    // Momentum Term
    Eigen::VectorXd M11 = rom.pCommonMatrices->BTotal * a * rom.nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = rom.pCommonMatrices->K * b;
    // Mass Term
    Eigen::VectorXd M5 = rom.pCommonMatrices->M * a_dot;
    // Pressure Term - Laplacian
    Eigen::VectorXd M3 = rom.pPPEMatrices->D * b;
    // Buoyancy Term - Momentum equation
    Eigen::VectorXd M10 = rom.pCommonMatrices->H * c;
    // Buoyancy term - PPE Equation
    Eigen::VectorXd M12 = rom.pPPEMatrices->HP * c;
    // BC Term - PPE Equation
    Eigen::VectorXd M7 = rom.pPPEMatrices->nuBC * a * rom.nu;
    // Time dep BC term - PPE Equation
    Eigen::VectorXd M13 = rom.pPPEMatrices->timedepBC * a_dot;
    // diffusive term temperature
    Eigen::VectorXd M6 = rom.pCommonMatrices->Y * c * (rom.nu / rom.Pr);
    // Mass Term Temperature
    Eigen::VectorXd M8 = rom.pCommonMatrices->W * c_dot;
    // Penalty term velocity
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(rom.Nphi_u, rom.N_BC);
    Eigen::MatrixXd penaltyT = Eigen::MatrixXd::Zero(rom.Nphi_t, rom.N_BC_t);

    if (rom.problem->bcMethod == "penalty")
    {
        for (int l = 0; l < rom.N_BC; l++)
        {
            penaltyU.col(l) = rom.tauU(l, 0) * (rom.boundaryConditions.currentVelocityBC(l) * rom.pPenaltyMatrices->bcVelVec[l] - rom.pPenaltyMatrices->bcVelMat[l] * a);
        }
        for (int l = 0; l < rom.N_BC_t; l++)
        {
            penaltyT.col(l) = rom.tauT(l, 0) * (rom.boundaryConditions.currentTemperatureBC(l) * rom.pPenaltyMatrices->bcTempVec.col(l) - Eigen::SliceFromTensor(rom.pPenaltyMatrices->bcTempMat, 0, l) * c);
        }
    }
    for (int i = 0; i < rom.Ntest_u ; i++)
    {
        cc = a.transpose() * Eigen::SliceFromTensor(rom.pCommonMatrices->C, 0, i) * a;
        ct = rom.nut_fluct.transpose() * Eigen::SliceFromTensor(rom.pCommonMatrices->CTotal, 0, i) * a;
        caveraged = rom.nut_param.transpose() * Eigen::SliceFromTensor(rom.pCommonMatrices->CTotalAve, 0, i) * a;
        residual(i) = -M5(i) + M11(i) - cc(0, 0) - M10(i) - M2(i) + ct(0, 0) + caveraged(0, 0);

        if (rom.problem->bcMethod == "penalty")
        {
            for (int l = 0; l < rom.N_BC; l++)
            {
                residual(i) += penaltyU(i, l);
            }
        }
    }
    for (int j = 0; j < rom.Ntest_prgh; j++)
    {
        int k = j + rom.Nphi_u;
        gg = a.transpose() * Eigen::SliceFromTensor(rom.pPPEMatrices->G, 0, j) * a;
        turbPPE = rom.nut_fluct.transpose() * Eigen::SliceFromTensor(rom.pPPEMatrices->CTotalPPEFluct, 0, j) * a + rom.nut_param.transpose() * Eigen::SliceFromTensor(rom.pPPEMatrices->CTotalPPEAve, 0, j) * a;
        residual(k) = M3(j, 0) + gg(0, 0) + M12(j, 0) - M7(j, 0) - turbPPE(0, 0);

        if (rom.problem->timedepbcMethod == "yes")
        {
            residual(k) += M13(j, 0);
        }
    }
    for (int j = 0; j < rom.Ntest_t; j++)
    {
        int k = j + rom.Nphi_u + rom.Nphi_prgh;
        qq = a.transpose() * Eigen::SliceFromTensor(rom.pCommonMatrices->Q, 0, j) * c;
        qt = rom.nut_fluct.transpose() * Eigen::SliceFromTensor(rom.pCommonMatrices->YTurb, 0, j) * c;
        qt_averaged = rom.nut_param.transpose() * Eigen::SliceFromTensor(rom.pCommonMatrices->AveYTurb, 0, j) * c;
        residual(k) = -M8(j) + M6(j) - qq(0, 0) + qt(0, 0) / rom.Pr_t + qt_averaged(0, 0) / rom.Pr_t;
        if (rom.problem->bcMethod == "penalty")
        {
            for (int l = 0; l < rom.N_BC_t; l++)
            {
                residual(k) += penaltyT(j, l);
            }
        }
    }

    if (rom.problem->bcMethod == "Gunzburger")
    {
        Eigen::MatrixXd gunzburgerBCProduct = a.transpose() * rom.problem->GunzburgerBCMatrixVelocity;
        Eigen::MatrixXd gunzburgerBCProductTemp = c.transpose() * rom.problem->GunzburgerBCMatrixTemperature;
        for (int j = 0; j < rom.N_BC; j++)
        {
            int idx = j + rom.Ntest_u;
            residual(idx) = gunzburgerBCProduct(0, j) - rom.boundaryConditions.currentVelocityBC(j);
        }
        for (int j = 0; j < rom.N_BC_t; j++)
        {
            int k = j + rom.Nphi_u + rom.Nphi_prgh + rom.Ntest_t;
            residual(k) = gunzburgerBCProductTemp(0, j) - rom.boundaryConditions.currentTemperatureBC(j);
        }
    }

    else if (rom.problem->bcMethod == "lift")
    {
        for (int j = 0; j < rom.N_BC; j++)
        {
            residual(j) = state(j) - rom.boundaryConditions.currentVelocityBC(j);
        }

        for (int j = 0; j < rom.N_BC_t; j++)
        {
            int k = j + rom.Nphi_u + rom.Nphi_prgh;
            residual(k) = state(k) - rom.boundaryConditions.currentTemperatureBC(j);
        }
    }
}

// ************************************************************************ //
