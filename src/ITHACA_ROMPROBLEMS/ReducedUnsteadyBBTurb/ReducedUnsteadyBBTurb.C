/*---------------------------------------------------------------------------*\
v/*---------------------------------------------------------------------------*\
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
/// Source file of the ReducedUnsteadyBBTurb class


#include "ReducedUnsteadyBBTurb.H"

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
ReducedUnsteadyBBTurb::ReducedUnsteadyBBTurb()
{
}

ReducedUnsteadyBBTurb::ReducedUnsteadyBBTurb(UnsteadyBBTurb& FOMproblem):
    problem(&FOMproblem)
{
    N_BC = problem->inletIndex.rows();
    N_BC_t = problem->inletIndexT.rows();
    Nphi_u = problem->B_matrix.rows();
    Nphi_prgh = problem->K_matrix.cols();
    Nphi_t = problem->Y_matrix.rows();
    Nphi_nut = problem->C_total_tensor.dimension(1);
    dimA = problem->dimA;
    newton_object_sup = newton_unsteadyBBTurb_sup(Nphi_u + Nphi_prgh + Nphi_t,
        Nphi_u + Nphi_prgh + Nphi_t, FOMproblem);
    newton_object_VMB = newton_unsteadyBBTurb_VMB(Nphi_u + Nphi_t, Nphi_u + Nphi_t, FOMproblem);
    newton_object_PPE = newton_unsteadyBBTurb_PPE(Nphi_u + Nphi_prgh + Nphi_t,
        Nphi_u + Nphi_prgh + Nphi_t, FOMproblem);

    // Create locally the velocity modes, with lifting and supremizer. Now it's commented
    // but maybe this could be useful so that the FOM problem can be deleted to free memory.

    // for (int k = 0; k < problem->liftfield.size(); k++)
    // {
    //     LUmodes.append((problem->liftfield[k]).clone());
    // }

    // for (int k = 0; k < problem->NUmodes; k++)
    // {
    //     LUmodes.append((problem->Umodes[k]).clone());
    // }

    // for (int k = 0; k < problem->NSUPmodes; k++)
    // {
    //     LUmodes.append((problem->supmodes[k]).clone());
    // }

    // // Create locally the Prgh modes - Probably unnecessary
    // for (int k = 0; k < problem->NPrghmodes; k++)
    // {
    //     Prghmodes.append((problem->P_rghmodes[k]).clone());
    // }

    // // Create locally the temperature modes including BC with liftfield
    // for (int k = 0; k < problem->liftfieldT.size(); k++)
    // {
    //     LTmodes.append((problem->liftfieldT[k]).clone());
    // }

    // for (int k = 0; k < problem->NTmodes; k++)
    // {
    //     LTmodes.append((problem->Tmodes[k]).clone());
    // }
}

// * * * * * * * * * * * * * * * Newton Object Setters * * * * * * * * * * * * * * * * //

template <typename NewtonType>
void ReducedUnsteadyBBTurb::setNewtonObjectBC(NewtonType& newtonObject, BoundaryConditions& boundaryConditions)
{
    newtonObject.BC.head(N_BC) = boundaryConditions.currentVelocityBC.head(N_BC);
    newtonObject.BC_t.head(N_BC_t) = boundaryConditions.currentTemperatureBC.head(N_BC_t);
}

template <typename NewtonType>
void ReducedUnsteadyBBTurb::configureNewtonObject(NewtonType& newtonObject, const Eigen::VectorXd& currentY, BoundaryConditions& boundaryConditions)
{
    newtonObject.y_old = currentY;
    newtonObject.dt = dt;
    newtonObject.Pr = Pr;
    newtonObject.Pr_t = Pr_t;
    newtonObject.BC_t.resize(N_BC_t);
    newtonObject.BC.resize(N_BC);
    newtonObject.tauU = tauU;
    newtonObject.tauT = tauT;
    newtonObject.nu = nu;
    // newtonObject.nu_fluct = nut0;
    // newtonObject.nu_param = nut_param_0;

    setNewtonObjectBC<NewtonType>(newtonObject, boundaryConditions);
}


// * * * * * * * * * * * * * * * Operators supremizer  * * * * * * * * * * * * * //
// Operator to evaluate the residual for the supremizer approach
int newton_unsteadyBBTurb_sup::operator()(const Eigen::VectorXd& x,
    Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_prgh);
    Eigen::VectorXd c_dot(Nphi_t);
    Eigen::VectorXd c_tmp(Nphi_t);
    a_tmp = x.head(Nphi_u);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    b_tmp = x.segment(Nphi_u, Nphi_prgh);
    c_tmp = x.tail(Nphi_t);
    c_dot = (x.tail(Nphi_t) - y_old.tail(Nphi_t)) / dt;
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Turbulence term
    Eigen::MatrixXd ct(1, 1);
    // Averaged turbulence term;
    Eigen::MatrixXd caveraged(1, 1);
    // Total Diffusive Term - In the laminar version it's just B_matrix.
    Eigen::VectorXd M11 = problem->B_total_matrix * a_tmp * nu;
    // Mass Term Velocity
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Continuity
    Eigen::VectorXd M3 = problem->P_matrix * a_tmp;
    // Buoyancy Term
    Eigen::VectorXd M10 = problem->H_matrix * c_tmp;
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    // // Convective term temperature (turbulence)
    // Eigen::MatrixXd qt(1, 1);
    // // Convective term temperature averaged (turbulence)
    // Eigen::MatrixXd qt_averaged(1, 1);
    // diffusive term temperature
    Eigen::VectorXd M6 = problem->Y_matrix * c_tmp * (nu / Pr);
    // Mass Term Temperature
    Eigen::VectorXd M8 = problem->W_matrix * c_dot;
    // Penalty term velocity
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    Eigen::MatrixXd penaltyT = Eigen::MatrixXd::Zero(Nphi_t, N_BC_t);
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = tauU(l, 0) * (BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * a_tmp);
        }
        for (int l = 0; l < N_BC_t; l++)
        {
            penaltyT.col(l) = tauT(l, 0) * (BC_t(l) * problem->bcTempVec[l] - problem->bcTempMat[l] * c_tmp);
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * a_tmp;
          // ct = nu_fluct.transpose() * Eigen::SliceFromTensor(problem->C_total_tensor, 0, i) * a_tmp;
          // caveraged = nu_param.transpose() * Eigen::SliceFromTensor(problem->C_total_ave_tensor, 0, i) * a_tmp;
        fvec(i) = -M5(i) + M11(i) - M2(i) - cc(0, 0) - M10(i); // + ct(0, 0) + caveraged(0, 0);  //DEBUG: to disable turbulence

        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC; l++)
            {
                fvec(i) += penaltyU(i, l);
            }
        }
    }

    for (int j = 0; j < Nphi_prgh; j++)
    {
        int k = j + Nphi_u;
        fvec(k) = M3(j);
    }

    for (int j = 0; j < Nphi_t; j++)
    {
        int k = j + Nphi_u + Nphi_prgh;
        qq = a_tmp.transpose() * Eigen::SliceFromTensor(problem->Q_tensor, 0, j) * c_tmp;
        // qt = nu_fluct.transpose() * Eigen::SliceFromTensor(problem->YT_tensor, 0, j) * c_tmp;
        // qt_averaged = nu_param.transpose() * Eigen::SliceFromTensor(problem->YT_ave_tensor, 0, j) * c_tmp;
        fvec(k) = -M8(j) + M6(j) - qq(0, 0);// + qt(0, 0) / Pr_t + qt_averaged(0, 0) / Pr_t; //DEBUG: to disable turbulence

        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC_t; l++)
            {
                fvec(k) += penaltyT(j, l);
            }
        }
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);
        }

        for (int j = 0; j < N_BC_t; j++)
        {
            int k = j + Nphi_u + Nphi_prgh;
            fvec(k) = x(k) - BC_t(j);
        }
    }
    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyBBTurb_sup::df(const Eigen::VectorXd& x,
    Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyBBTurb_sup> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators VMB * * * * * * * * * * * * * * //
int newton_unsteadyBBTurb_VMB::operator()(const Eigen::VectorXd& x,
    Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd c_dot(Nphi_t);
    Eigen::VectorXd c_tmp(Nphi_t);
    a_tmp = x.head(Nphi_u);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    c_tmp = x.tail(Nphi_t);
    c_dot = (x.tail(Nphi_t) - y_old.tail(Nphi_t)) / dt;
    // Convective term
    Eigen::MatrixXd cc(1, 1);
    // Turbulence term
    Eigen::MatrixXd ct(1, 1);
    // Averaged turbulence term;
    Eigen::MatrixXd caveraged(1, 1);
    // Total Diffusive Term - In the laminar version it's just B_matrix.
    Eigen::VectorXd M11 = problem->B_total_matrix * a_tmp * nu;
    // Mass Term Velocity
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * a_tmp;
    // Buoyancy Term
    Eigen::VectorXd M10 = problem->H_matrix * c_tmp;
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    // Convective term temperature (turbulence)
    Eigen::MatrixXd qt(1, 1);
    // Convective term temperature averaged (turbulence)
    Eigen::MatrixXd qt_averaged(1, 1);
    // diffusive term temperature
    Eigen::VectorXd M6 = problem->Y_matrix * c_tmp * (nu / Pr);
    // Mass Term Temperature
    Eigen::VectorXd M8 = problem->W_matrix * c_dot;
    // Penalty term velocity
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    Eigen::MatrixXd penaltyT = Eigen::MatrixXd::Zero(Nphi_t, N_BC_t);
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = tauU(l, 0) * (BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * a_tmp);
        }
        for (int l = 0; l < N_BC_t; l++)
        {
            penaltyT.col(l) = tauT(l, 0) * (BC_t(l) * problem->bcTempVec[l] - problem->bcTempMat[l] * c_tmp);
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * a_tmp;
        ct = nu_fluct.transpose() * Eigen::SliceFromTensor(problem->C_total_tensor, 0, i) * a_tmp;
        caveraged = nu_param.transpose() * Eigen::SliceFromTensor(problem->C_total_ave_tensor, 0, i) * a_tmp;
        fvec(i) = -M5(i) + M11(i) - M2(i) - cc(0, 0) - M10(i) + ct(0, 0) + caveraged(0, 0);

        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC; l++)
            {
                fvec(i) += penaltyU(i, l);
            }
        }
    }

    for (int j = 0; j < Nphi_t; j++)
    {
        int k = j + Nphi_u;
        qq = a_tmp.transpose() * Eigen::SliceFromTensor(problem->Q_tensor, 0, j) * c_tmp;
        qt = nu_fluct.transpose() * Eigen::SliceFromTensor(problem->YT_tensor, 0, j) * c_tmp;
        qt_averaged = nu_param.transpose() * Eigen::SliceFromTensor(problem->YT_ave_tensor, 0, j) * c_tmp;
        fvec(k) = -M8(j) + M6(j) - qq(0, 0) + qt(0, 0) / Pr_t + qt_averaged(0, 0) / Pr_t;

        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC_t; l++)
            {
                fvec(k) += penaltyT(j, l);
            }
        }
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);
        }

        for (int j = 0; j < N_BC_t; j++)
        {
            int k = j + Nphi_u;
            fvec(k) = x(k) - BC_t(j);
        }
    }
    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyBBTurb_VMB::df(const Eigen::VectorXd& x,
    Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyBBTurb_VMB> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Operators PPE * * * * * * * * * * * * * * * //
int newton_unsteadyBBTurb_PPE::operator()(const Eigen::VectorXd& x,
    Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd a_dot(Nphi_u);
    Eigen::VectorXd a_tmp(Nphi_u);
    Eigen::VectorXd b_tmp(Nphi_prgh);
    Eigen::VectorXd c_dot(Nphi_t);
    Eigen::VectorXd c_tmp(Nphi_t);
    a_tmp = x.head(Nphi_u);
    b_tmp = x.segment(Nphi_u, Nphi_prgh);
    a_dot = (x.head(Nphi_u) - y_old.head(Nphi_u)) / dt;
    c_tmp = x.tail(Nphi_t);
    c_dot = (x.tail(Nphi_t) - y_old.tail(Nphi_t)) / dt;
    // Convective terms
    Eigen::MatrixXd cc(1, 1);
    Eigen::MatrixXd gg(1, 1);
    Eigen::MatrixXd bb(1, 1);
    // Convective term temperature
    Eigen::MatrixXd qq(1, 1);
    // Momentum Term
    Eigen::VectorXd M11 = problem->B_total_matrix * a_tmp * nu;
    // Gradient of pressure
    Eigen::VectorXd M2 = problem->K_matrix * b_tmp;
    // Mass Term
    Eigen::VectorXd M5 = problem->M_matrix * a_dot;
    // Pressure Term - Laplacian
    Eigen::VectorXd M3 = problem->D_matrix * b_tmp;
    // BC PPE
    Eigen::VectorXd M6 = problem->BC1_matrix * a_tmp * nu;
    // BC PPE
    // Buoyancy Term
    Eigen::VectorXd M10 = problem->H_matrix * c_tmp;
    Eigen::VectorXd M12 = problem->HP_matrix * c_tmp;
    Eigen::VectorXd M7 = problem->BC3_matrix * a_tmp * nu;
    // diffusive term temperature
    Eigen::VectorXd M9 = problem->Y_matrix * c_tmp * (nu / Pr);
    // Mass Term Temperature
    Eigen::VectorXd M8 = problem->W_matrix * c_dot;
    // Penalty term velocity
    Eigen::MatrixXd penaltyU = Eigen::MatrixXd::Zero(Nphi_u, N_BC);
    Eigen::MatrixXd penaltyT = Eigen::MatrixXd::Zero(Nphi_t, N_BC_t);
    if (problem->bcMethod == "penalty")
    {
        for (int l = 0; l < N_BC; l++)
        {
            penaltyU.col(l) = tauU(l, 0) * (BC(l) * problem->bcVelVec[l] - problem->bcVelMat[l] * a_tmp);
        }
        for (int l = 0; l < N_BC_t; l++)
        {
            penaltyT.col(l) = tauT(l, 0) * (BC_t(l) * problem->bcTempVec[l] - problem->bcTempMat[l] * c_tmp);
        }
    }

    for (int i = 0; i < Nphi_u; i++)
    {
        cc = a_tmp.transpose() * Eigen::SliceFromTensor(problem->C_tensor, 0, i) * a_tmp;
        fvec(i) = -M5(i) + M11(i) - cc(0, 0) - M10(i) - M2(i);

        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC; l++)
            {
                fvec(i) += penaltyU(i, l);
            }
        }
    }
    for (int j = 0; j < Nphi_prgh; j++)
    {
        int k = j + Nphi_u;
        gg = a_tmp.transpose() * Eigen::SliceFromTensor(problem->G_tensor, 0, j) * a_tmp;
        // bb = a_tmp.transpose() * Eigen::SliceFromTensor(problem->BC2_tensor, 0, j) * a_tmp;
        fvec(k) = M3(j, 0) + gg(0, 0) + M12(j, 0) - M7(j, 0);
    }
    for (int j = 0; j < Nphi_t; j++)
    {
        int k = j + Nphi_u + Nphi_prgh;
        qq = a_tmp.transpose() * Eigen::SliceFromTensor(problem->Q_tensor, 0, j) * c_tmp;
        // qt = nu_fluct.transpose() * Eigen::SliceFromTensor(problem->YT_tensor, 0, j) * c_tmp;
        // qt_averaged = nu_param.transpose() * Eigen::SliceFromTensor(problem->YT_ave_tensor, 0, j) * c_tmp;
        fvec(k) = -M8(j) + M6(j) - qq(0, 0);// + qt(0, 0) / Pr_t + qt_averaged(0, 0) / Pr_t; //DEBUG: to disable turbulence
        Info << "### PPE Temp Residual at mode " << j << " : " << fvec(k) << endl;
        if (problem->bcMethod == "penalty")
        {
            for (int l = 0; l < N_BC_t; l++)
            {
                fvec(k) += penaltyT(j, l);
                Info << "### PPE Temp Penalty at mode " << j << " from BC " << l << " : " << penaltyT(j, l) << endl;
            }
        }
    }

    if (problem->bcMethod == "lift")
    {
        for (int j = 0; j < N_BC; j++)
        {
            fvec(j) = x(j) - BC(j);
        }

        for (int j = 0; j < N_BC_t; j++)
        {
            int k = j + Nphi_u + Nphi_prgh;
            fvec(k) = x(k) - BC_t(j);
        }
    }
    return 0;
}

// Operator to evaluate the Jacobian for the supremizer approach
int newton_unsteadyBBTurb_PPE::df(const Eigen::VectorXd& x,
    Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<newton_unsteadyBBTurb_PPE> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// * * * * * * * * * * * * * * * Solve Functions  * * * * * * * * * * * * * //
void ReducedUnsteadyBBTurb::solveOnline_sup(const Eigen::MatrixXd& temperatureBC,
    const Eigen::MatrixXd& velocityBC, int NParaSet, int startSnap)
{
    Info << "################## Online solve N° " << NParaSet << " ##################" << endl;
    validateSettings();
    BoundaryConditions boundaryConditions(velocityBC, temperatureBC, "linear");

    y.resize(Nphi_u + Nphi_prgh + Nphi_t, 1);
    y.setZero();

    boundaryConditions.initializeReducedCoeffs(startSnap, y, problem, Nphi_u, Nphi_prgh, Nphi_t, N_BC_t);

    if (problem->bcMethod == "lift")
    {
        boundaryConditions.correctLiftingCoeffs(y, N_BC, N_BC_t, Nphi_u, Nphi_prgh);
    }

    // nut0 = ITHACAutilities::getCoeffs(problem->fluctNutfield[startSnap], problem->nutmodes);
    // nut_param_0 = interpolateIDW();

    int firstRBFInd;
    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFInd = problem->skipRBFIndex;
        Info << "Skipping lifting modes in the RBF evaluation. This means that the first RBF index is set to " << firstRBFInd << endl;
    } else
    {
        firstRBFInd = 0;
        Info << "## COMM - Not skipping lifting modes in the RBF evaluation. This means that the first RBF index is set to " << firstRBFInd << endl;
    }

    configureNewtonObject(newton_object_sup, y, boundaryConditions);

    // Set number of online solutions - Time related
    int numberOfStores = round((storeEvery) / dt); // Number of time steps between two stored solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores); // Total stored solutions, excluding initial condition
    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(Nphi_nut + 1, onlineSize + 3);
    time = tstart;
    int nextStore = 0;
    int timeStepCounter = 0;
    int storedSnapshotsCounter = 0;

    // Creates a vector to store the temporal solution, while saving the initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_prgh + Nphi_t + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    // // This part up to the while loop is just to compute the eddy viscosity field at time = startTime if time is not equal to zero
    // if ((time != 0) || (startFromZero == true))
    // {
    //     online_solution[timeStepCounter] = tmp_sol;
    //     timeStepCounter++;
    //     rbfCoeffMat(0, storedSnapshotsCounter) = time;
    //     rbfCoeffMat.block(1, storedSnapshotsCounter, Nphi_nut, 1) = nut0;
    //     storedSnapshotsCounter++;
    //     nextStore += numberOfStores;
    // }

    Eigen::HybridNonLinearSolver<newton_unsteadyBBTurb_sup> hnls(newton_object_sup);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // Eigen::VectorXd tv;
    // Eigen::VectorXd aDer;
    // tv.resize(dimA);

    while (time < finalTime - dt)
    {
        time = time + dt;

        boundaryConditions.updateTimeDependentBC(time);

        if (problem->timedepbcMethod == "yes")
        {
            setNewtonObjectBC(newton_object_sup, boundaryConditions);
        }

        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);

        // tv.setZero();
        // aDer.setZero();
        // aDer = (y.head(Nphi_u) - newton_object_sup.y_old.head(Nphi_u)) / dt;
        // tv << y.segment(firstRBFInd, dimA / 2), aDer.segment(firstRBFInd, dimA / 2);

        /// Now normalize each tv entry with the eigen::vectorXd problem->meanA and problem->stdA
        // for (int i = 0; i < dimA; i++)
        // {
        //     tv(i) = (tv(i) - problem->meanA(i)) / problem->stdA(i);
        // }

        // for (int j = 0; j < Nphi_nut; j++)
        // {
        //     newton_object_sup.nu_fluct(j) = problem->rbfSplines[j]->eval(tv); // * problem->stdG(j) + problem->meanG(j);
        // }

        if (problem->bcMethod == "lift")
        {
            boundaryConditions.correctLiftingCoeffs(y, N_BC, N_BC_t, Nphi_u, Nphi_prgh);
        }

        newton_object_sup.operator()(y, res);
        newton_object_sup.y_old = y;

        Info << "Time = " << time << endl;
        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl
                      << std::endl;
        } else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl
                      << std::endl;
        }

        count_online_solve += 1;
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
            // rbfCoeffMat(0, storedSnapshotsCounter) = time;
            // rbfCoeffMat.block(1, storedSnapshotsCounter, Nphi_nut, 1) = newton_object_sup.nu_fluct;
            nextStore += numberOfStores;
            storedSnapshotsCounter++;
        }
        timeStepCounter++;
    }
    if (Pstream::master())
    {
        ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
            "./ITHACAoutput/red_coeff_" + name(NParaSet) + "/");
    }
    Info << "Online solve finished, total time steps solved: " << timeStepCounter << endl;
    count_online_solve += 1;
}

void ReducedUnsteadyBBTurb::solveOnline_PPE(const Eigen::MatrixXd& temperatureBC,
    const Eigen::MatrixXd& velocityBC, int NParaSet, int startSnap)
{
    Info << "################## Online solve N° " << NParaSet << " ##################" << endl;
    validateSettings();
    BoundaryConditions boundaryConditions(velocityBC, temperatureBC, "linear");

    y.resize(Nphi_u + Nphi_prgh + Nphi_t, 1);
    y.setZero();

    boundaryConditions.initializeReducedCoeffs(startSnap, y, problem, Nphi_u, Nphi_prgh, Nphi_t, N_BC_t);
    if (problem->bcMethod == "lift")
    {
        boundaryConditions.correctLiftingCoeffs(y, N_BC, N_BC_t, Nphi_u, Nphi_prgh);
    }
    // nut0 = ITHACAutilities::getCoeffs(problem->fluctNutfield[startSnap], problem->nutmodes);
    // nut_param_0 = interpolateIDW();

    int firstRBFInd;
    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFInd = problem->skipRBFIndex;
        Info << "Skipping lifting modes in the RBF evaluation. This means that the first RBF index is set to " << firstRBFInd << endl;
    } else
    {
        firstRBFInd = 0;
        Info << "## COMM - Not skipping lifting modes in the RBF evaluation. This means that the first RBF index is set to " << firstRBFInd << endl;
    }

    configureNewtonObject(newton_object_PPE, y, boundaryConditions);

    // Set number of online solutions - Time related
    int numberOfStores = round((storeEvery) / dt); // Number of time steps between two stored solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores); // Total stored solutions, excluding initial condition
    online_solution.resize(onlineSize);
    // rbfCoeffMat.resize(Nphi_nut + 1, onlineSize + 3);
    time = tstart;
    int nextStore = 0;
    int timeStepCounter = 0;
    int storedSnapshotsCounter = 0;

    // Creates a vector to store the temporal solution, while saving the initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_prgh + Nphi_t + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    Eigen::HybridNonLinearSolver<newton_unsteadyBBTurb_PPE> hnls(newton_object_PPE);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    // Eigen::VectorXd tv;
    // Eigen::VectorXd aDer;
    // tv.resize(dimA);
    while (time < finalTime - dt)
    {
        time = time + dt;
        boundaryConditions.updateTimeDependentBC(time);

        if (problem->timedepbcMethod == "yes")
        {
            setNewtonObjectBC(newton_object_PPE, boundaryConditions);
        }

        Eigen::VectorXd res(y);
        res.setZero();
        Info << "### DEBUG: Solving the linear system ###" << endl;
        hnls.solve(y);

        // tv.setZero();
        // aDer.setZero();
        // aDer = (y.head(Nphi_u) - newton_object_sup.y_old.head(Nphi_u)) / dt;
        // tv << y.segment(firstRBFInd, dimA / 2), aDer.segment(firstRBFInd, dimA / 2);

        /// Now normalize each tv entry with the eigen::vectorXd problem->meanA and problem->stdA
        // for (int i = 0; i < dimA; i++)
        // {
        //     tv(i) = (tv(i) - problem->meanA(i)) / problem->stdA(i);
        // }

        // for (int j = 0; j < Nphi_nut; j++)
        // {
        //     newton_object_sup.nu_fluct(j) = problem->rbfSplines[j]->eval(tv); // * problem->stdG(j) + problem->meanG(j);
        // }

        if (problem->bcMethod == "lift")
        {
            boundaryConditions.correctLiftingCoeffs(y, N_BC, N_BC_t, Nphi_u, Nphi_prgh);
        }

        newton_object_PPE.operator()(y, res);
        newton_object_PPE.y_old = y;

        Info << "Time = " << time << endl;
        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl
                      << std::endl;
        } else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl
                      << std::endl;
        }

        count_online_solve += 1;
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
            // rbfCoeffMat(0, storedSnapshotsCounter) = time;
            // rbfCoeffMat.block(1, storedSnapshotsCounter, Nphi_nut, 1) = newton_object_sup.nu_fluct;
            nextStore += numberOfStores;
            storedSnapshotsCounter++;
        }
        timeStepCounter++;
    }
    if (Pstream::master())
    {
        ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
            "./ITHACAoutput/red_coeff_" + name(NParaSet) + "/");
    }
    Info << "Online solve finished, total time steps solved: " << timeStepCounter << endl;
    count_online_solve += 1;
}

void ReducedUnsteadyBBTurb::solveOnline_VMB(const Eigen::MatrixXd& temperatureBC,
    const Eigen::MatrixXd& velocityBC, int NParaSet, int startSnap)
{
    pressureMethod = "VMB"; // Set the flag to use VMB method in the post-processing
    validateSettings();
    BoundaryConditions boundaryConditions(velocityBC, temperatureBC, "linear");

    Info << "################## Online solve N° " << NParaSet << " ##################" << endl;

    y.resize(Nphi_u + Nphi_t, 1); // Solution vector
    y.setZero();

    volScalarField T_IC("T_IC", problem->Tfield[startSnap]);
    for (int j = 0; j < T_IC.boundaryField().size(); j++)
    {
        for (int i = 0; i < N_BC_t; i++)
        {
            if (j == problem->inletIndexT(i, 0))
            {
                T_IC.boundaryFieldRef()[problem->inletIndexT(i, 0)][j] = boundaryConditions.currentTemperatureBC(i);
            } else
            {
            }
        }
    }
    y.head(Nphi_u) = ITHACAutilities::getCoeffs(problem->Ufield[startSnap], problem->L_U_SUPmodes);
    y.tail(Nphi_t) = ITHACAutilities::getCoeffs(T_IC, problem->L_Tmodes);


    // Change initial condition for the lifting function
    if (problem->bcMethod == "lift")
    {
        for (int i = 0; i < N_BC; i++)
        {
            y(i) = boundaryConditions.currentVelocityBC(i);
        }
        for (int i = 0; i < N_BC_t; i++)
        {
            int k = i + Nphi_u;
            y(k) = boundaryConditions.currentTemperatureBC(i);
        }
    }

    nut0 = ITHACAutilities::getCoeffs(problem->fluctNutfield[startSnap], problem->nutmodes);
    nut_param_0 = interpolateIDW();

    int firstRBFInd;
    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFInd = problem->skipRBFIndex;
        Info << "Skipping lifting modes in the RBF evaluation. This means that the first RBF index is set to " << firstRBFInd << endl;
    } else
    {
        firstRBFInd = 0;
        Info << "## COMM - Not skipping lifting modes in the RBF evaluation. This means that the first RBF index is set to " << firstRBFInd << endl;
    }

    configureNewtonObject(newton_object_VMB, y, boundaryConditions);

    // Set number of online solutions - Time related
    int numberOfStores = round((storeEvery) / dt); // Number of time steps between two stored solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores); // Total stored solutions, excluding initial condition
    online_solution.resize(onlineSize);
    rbfCoeffMat.resize(Nphi_nut + 1, onlineSize + 3);
    time = tstart;
    int nextStore = 0;
    int timeStepCounter = 0;
    int storedSnapshotsCounter = 0;

    // Creates a vector to store the temporal solution, while saving the initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_t + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    // This part up to the while loop is just to compute the eddy viscosity field at time = startTime if time is not equal to zero
    if ((time != 0) || (startFromZero == true))
    {
        online_solution[timeStepCounter] = tmp_sol;
        timeStepCounter++;
        rbfCoeffMat(0, storedSnapshotsCounter) = time;
        rbfCoeffMat.block(1, storedSnapshotsCounter, Nphi_nut, 1) = nut0;
        storedSnapshotsCounter++;
        nextStore += numberOfStores;
    }

    Eigen::HybridNonLinearSolver<newton_unsteadyBBTurb_VMB> hnls(newton_object_VMB);
    // Set output colors for fancy output
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);

    Eigen::VectorXd tv;
    Eigen::VectorXd aDer;
    tv.resize(dimA);
    while (time < finalTime)
    {
        time = time + dt;

        boundaryConditions.updateTimeDependentBC(time);

        if (problem->timedepbcMethod == "yes")
        {
            setNewtonObjectBC(newton_object_VMB, boundaryConditions);
        }
        Eigen::VectorXd res(y);
        res.setZero();
        hnls.solve(y);

        tv.setZero();
        aDer.setZero();
        aDer = (y.head(Nphi_u) - newton_object_VMB.y_old.head(Nphi_u)) / dt;
        tv << y.segment(firstRBFInd, dimA / 2), aDer.segment(firstRBFInd, dimA / 2);

        /// Now normalize each tv entry with the eigen::vectorXd problem->meanA and problem->stdA
        // for (int i = 0; i < dimA; i++)
        // {
        //     tv(i) = (tv(i) - problem->meanA(i)) / problem->stdA(i);
        // }

        for (int j = 0; j < Nphi_nut; j++)
        {
            newton_object_VMB.nu_fluct(j) = problem->rbfSplines[j]->eval(tv); // * problem->stdG(j) + problem->meanG(j);
        }

        if (problem->bcMethod == "lift")
        {
            boundaryConditions.correctLiftingCoeffs(y, N_BC, N_BC_t, Nphi_u, 0);
        }

        newton_object_VMB.operator()(y, res);
        newton_object_VMB.y_old = y;
        Info << "######### Online solve N° " << count_online_solve << " ##########" << endl;
        Info << "Time = " << time << endl;
        if (res.norm() < 1e-5)
        {
            std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl
                      << std::endl;
        } else
        {
            std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl
                      << std::endl;
        }

        count_online_solve += 1;
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
            rbfCoeffMat.block(1, storedSnapshotsCounter, Nphi_nut, 1) = newton_object_VMB.nu_fluct;
            nextStore += numberOfStores;
            storedSnapshotsCounter++;
        }
        timeStepCounter++;
    }
    if (Pstream::master())
    {
        ITHACAstream::exportMatrix(online_solution, "red_coeff", "python",
            "./ITHACAoutput/red_coeff_" + name(NParaSet) + "/");
    }
    Info << "Online solve finished, total time steps solved: " << timeStepCounter << endl;
    count_online_solve += 1;
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
    // List<Eigen::MatrixXd> CoeffNut;
    CoeffU.resize(0);
    CoeffPrgh.resize(0);
    CoeffT.resize(0);
    // CoeffNut.resize(0);

    for (int i = 0; i < online_solution.size(); i++)
    {
        if (timeStepCounter == nextWrite)
        {
            Eigen::MatrixXd currentUCoeff;
            Eigen::MatrixXd currentPrghCoeff;
            Eigen::MatrixXd currentTCoeff;
            // Eigen::MatrixXd currentNutCoeff;

            currentUCoeff = online_solution[i].block(1, 0, Nphi_u, 1);
            if (pressureMethod != "VMB")
            {
                currentPrghCoeff = online_solution[i].block(Nphi_u + 1, 0, Nphi_prgh, 1);
            } else
            {
                currentPrghCoeff = online_solution[i].block(1, 0, Nphi_u, 1);
            }
            currentTCoeff = online_solution[i].bottomRows(Nphi_t);
            // currentNutCoeff = rbfCoeffMat.block(1, i, Nphi_nut, 1);

            CoeffPrgh.append(currentPrghCoeff);
            CoeffU.append(currentUCoeff);
            CoeffT.append(currentTCoeff);
            // CoeffNut.append(currentNutCoeff);

            nextWrite += exportEveryIndex;
        }
        timeStepCounter++;
    }
    volVectorField uRec("uRec", problem->L_U_SUPmodes[0]);
    volScalarField TRec("TRec", problem->L_Tmodes[0]);
    volScalarField prghRec("prghRec", problem->P_rghmodes[0]);
    // volScalarField nutFluctRec("nutFluctRec", problem->nutmodes[0]);

    // TODO: Implement a reconstruction of the field p_rgh using the reconstruc method instead of the messy one below

    uRecFields = problem->L_U_SUPmodes.reconstruct(uRec, CoeffU, "uRec");
    TRecFields = problem->L_Tmodes.reconstruct(TRec, CoeffT, "TRec");
    // nutFluctRecFields = problem->nutmodes.reconstruct(nutFluctRec, CoeffNut, "nutFluctRec");

    // Reconstruct the averaged eddy viscosity field as a linear combination of the nut_param_0 coefficients and the PtrList<volScalarField> avgNutfield;

    // volScalarField nutAvg(
    //     IOobject(
    //         "nutAvgRec",
    //         problem->nutmodes[0].time().timeName(),
    //         problem->nutmodes[0].mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE),
    //     problem->nutmodes[0].mesh(),
    //     dimensionedScalar("zero", problem->nutmodes[0].dimensions(), 0.0));

    // for (int k = 0; k < nut_param_0.size(); k++)
    // {
    //     nutAvg += nut_param_0(k) * problem->avgNutfield[k];
    // }

    // volScalarField nutRecField(
    //     IOobject(
    //         "nutRec",
    //         problem->nutmodes[0].time().timeName(),
    //         problem->nutmodes[0].mesh(),
    //         IOobject::NO_READ,
    //         IOobject::NO_WRITE),
    //     problem->nutmodes[0].mesh(),
    //     dimensionedScalar("zero", problem->nutmodes[0].dimensions(), 0.0));

    // PtrList<volScalarField> nutRecFields;

    // // Now sum the reconstructed fluctuating and averaged eddy viscosity fields
    // forAll(nutFluctRecFields, i)
    // {
    //     nutRecFields.append(nutFluctRecFields[i] + nutAvg);
    // }

    if (exportFields)
    {
        ITHACAstream::exportFields(uRecFields, folder, "uRec");
        ITHACAstream::exportFields(TRecFields, folder, "TRec");
        ITHACAstream::exportFields(nutFluctRecFields, folder, "nutFluctRec");
        // ITHACAstream::exportFields(nutRecFields, folder, "nutRec");
    }
    // TODO: Implement correct BC handling for shifted pressure reconstruction (fixedFluxPressure in theory requires a gradient evaluation I guess)
    // prghRecFields = problem->P_rghmodes.reconstruct(prghRec, CoeffPrgh, "prghRec");
    // if (exportFields)
    // {
    //     ITHACAstream::exportFields(prghRecFields, folder, "prghRec");
    // }
    for (int i = 0; i < CoeffPrgh.size(); i++)
    {
        volScalarField prghRec("prghRec", problem->P_rghmodes[0]);
        prghRec *= scalar(0);
        for (int j = 0; j < Nphi_prgh; j++)
        {
            prghRec += problem->P_rghmodes[j] * CoeffPrgh[i](j, 0);
        }
        ITHACAstream::exportSolution(prghRec, name(i + 1), folder);
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

// * * * * * * * * *  Inverse Distance Weighting Functions  * * * * * * * * //
void ReducedUnsteadyBBTurb::setTimeSettings(const Eigen::MatrixXd& timeMatrix, label solutionCounter)
{
    tstart = timeMatrix(solutionCounter, 0);
    finalTime = timeMatrix(solutionCounter, 1);
    dt = timeMatrix(solutionCounter, 2);
    storeEvery = timeMatrix(solutionCounter, 3);
    exportEvery = timeMatrix(solutionCounter, 4);
}

// * * * * * * * * *  Penalty Factor Estimation Functions  * * * * * * * * //
void ReducedUnsteadyBBTurb::estimatePenaltyFactorSupremizer(const Eigen::MatrixXd& velocityBC,
    const Eigen::MatrixXd& temperatureBC, int NtimeStepPenalty,
    int maxIter, double tolerancePenaltyU, double tolerancePenaltyT, int startSnap)
{
    M_Assert(tolerancePenaltyU > 0, "The tolerance for the penalty factor estimation must be positive.");
    M_Assert(tolerancePenaltyT > 0, "The tolerance for the penalty factor estimation must be positive.");
    M_Assert(tauU.size() == N_BC, "The size of the tauU vector must be equal to the number of velocity BCs.");
    M_Assert(tauT.size() == N_BC_t, "The size of the tauT vector must be equal to the number of temperature BCs.");

    Eigen::VectorXd currentTimeErrorsU = Eigen::VectorXd::Constant(N_BC, 1e6);
    Eigen::VectorXd currentTimeErrorsT = Eigen::VectorXd::Constant(N_BC_t, 1e6);
    Eigen::VectorXd currentTimeU = Eigen::VectorXd::Zero(N_BC);
    Eigen::VectorXd currentTimeT = Eigen::VectorXd::Zero(N_BC_t);
    Eigen::VectorXd toleranceVectorU = Eigen::VectorXd::Constant(N_BC, tolerancePenaltyU);
    Eigen::VectorXd toleranceVectorT = Eigen::VectorXd::Constant(N_BC_t, tolerancePenaltyT);

    BoundaryConditions boundaryConditions(velocityBC, temperatureBC, "linear");
    y.resize(Nphi_u + Nphi_prgh + Nphi_t, 1); // Solution vector
    y.setZero();

    boundaryConditions.initializeReducedCoeffs(startSnap, y, problem, Nphi_u, Nphi_prgh, Nphi_t, N_BC_t);
    // nut0 = ITHACAutilities::getCoeffs(problem->fluctNutfield[startSnap], problem->nutmodes);
    // nut_param_0 = interpolateIDW();

    if (problem->bcMethod == "lift")
    {
        boundaryConditions.correctLiftingCoeffs(y, N_BC, N_BC_t, Nphi_u, Nphi_prgh);
    }

    int firstRBFInd;
    if (skipLift == true && problem->bcMethod == "lift")
    {
        firstRBFInd = problem->skipRBFIndex;
    } else
    {
        firstRBFInd = 0;
    }

    configureNewtonObject(newton_object_sup, y, boundaryConditions);
    
    int numberOfStores = round((storeEvery) / dt); // Number of time steps between two stored solutions
    int Ntsteps = static_cast<int>((finalTime - tstart) / dt);
    int onlineSize = static_cast<int>(Ntsteps / numberOfStores); // Total stored solutions, excluding initial condition
    rbfCoeffMat.resize(Nphi_nut + 1, onlineSize + 3);
    time = tstart;

    // Creates a vector to store the temporal solution, while saving the initial condition as first solution
    Eigen::MatrixXd tmp_sol(Nphi_u + Nphi_prgh + Nphi_t + 1, 1);
    tmp_sol(0) = time;
    tmp_sol.col(0).tail(y.rows()) = y;

    Eigen::HybridNonLinearSolver<newton_unsteadyBBTurb_sup> hnls(newton_object_sup);
    Color::Modifier red(Color::FG_RED);
    Color::Modifier green(Color::FG_GREEN);
    Color::Modifier def(Color::FG_DEFAULT);
    volVectorField uRecPen("uRecPen", problem->L_U_SUPmodes[0]);
    volScalarField TRecPen("TRecPen", problem->L_Tmodes[0]);

    // Eigen::VectorXd tv;
    // Eigen::VectorXd aDer;
    // tv.resize(dimA);

    int currentTimeStep = 0;
    int currentIteration = 0;
    Info << "Penalty factor estimation - Current iteration: " << currentIteration << endl;
    while (currentTimeStep < NtimeStepPenalty)
    {
        currentTimeErrorsU.setConstant(1e6);
        currentTimeErrorsT.setConstant(1e6);
        currentIteration = 0;
        currentTimeStep += 1;
        
        bool convergedU = false;
        bool convergedT = false; // These flags are used so that, when one of the two variables has converged, its penalty factor is not updated anymore

        time = time + dt;
        boundaryConditions.updateTimeDependentBC(time);
        if (problem->timedepbcMethod == "yes")
        {
            setNewtonObjectBC(newton_object_sup, boundaryConditions);
        }

        while ((currentTimeErrorsU.maxCoeff() > tolerancePenaltyU || currentTimeErrorsT.maxCoeff() > tolerancePenaltyT) && currentIteration < maxIter)
        {
            currentIteration++;
            Eigen::VectorXd res(y);
            res.setZero();
            hnls.solve(y);

            // tv.setZero();
            // aDer.setZero();
            // aDer = (y.head(Nphi_u) - newton_object_sup.y_old.head(Nphi_u)) / dt;
            // tv << y.segment(firstRBFInd, dimA / 2), aDer.segment(firstRBFInd, dimA / 2);

            // for (int j = 0; j < Nphi_nut; j++)
            // {
            //     newton_object_sup.nu_fluct(j) = problem->rbfSplines[j]->eval(tv); // * problem->stdG(j) + problem->meanG(j);
            // }

            newton_object_sup.operator()(y, res);

            Info << "Time = " << time << endl;
            if (res.norm() < 1e-5)
            {
                std::cout << green << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl
                          << std::endl;
            } else
            {
                std::cout << red << "|F(x)| = " << res.norm() << " - Minimun reached in " << hnls.iter << " iterations " << def << std::endl
                          << std::endl;
            }

            if (convergedU == false)
            {
                uRecPen = problem->L_U_SUPmodes.reconstruct(uRecPen, y.head(Nphi_u), "uRecPen");
                for (int k = 0; k < N_BC; k++)
                {
                    currentTimeU(k) = uRecPen.boundaryFieldRef()[problem->inletIndex(k, 0)][0].component(problem->inletIndex(k, 1));
                }
                currentTimeErrorsU = (currentTimeU - boundaryConditions.currentVelocityBC).cwiseAbs();
                for (int k = 0; k < N_BC; k++)
                {
                    tauU(k, 0) = tauU(k, 0) * (currentTimeErrorsU(k) / tolerancePenaltyU);
                }
                newton_object_sup.tauU = tauU;
                if (currentTimeErrorsU.maxCoeff() <= tolerancePenaltyU)
                {
                    convergedU = true;
                }
            }

            if (convergedT == false)
            {
                TRecPen = problem->L_Tmodes.reconstruct(TRecPen, y.tail(Nphi_t), "TRecPen");
                for (int k = 0; k < N_BC_t; k++)
                {
                    currentTimeT(k) = TRecPen.boundaryFieldRef()[problem->inletIndexT(k, 0)][0];
                }
                currentTimeErrorsT = (currentTimeT - boundaryConditions.currentTemperatureBC).cwiseAbs();
                for (int k = 0; k < N_BC_t; k++)
                {
                    tauT(k, 0) = tauT(k, 0) * (currentTimeErrorsT(k) / tolerancePenaltyT);
                }
                newton_object_sup.tauT = tauT;
                if (currentTimeErrorsT.maxCoeff() <= tolerancePenaltyT)
                {
                    convergedT = true;
                }
            }
        }
        newton_object_sup.y_old = y;

    }
    Info << "Penalty factor estimation finished for timeStep " << currentTimeStep << " after " << currentIteration << " iterations." << endl;
    Info << "Final tauU values: " << tauU << endl;
    Info << "Final tauT values: " << tauT << endl;
    Info << "Final velocity BC errors: " << currentTimeErrorsU << endl;
    Info << "Final temperature BC errors: " << currentTimeErrorsT << endl;
}


// * * * * * * * * *  Validation and setup helpers  * * * * * * * * //
void ReducedUnsteadyBBTurb::validateSettings()
{
    Info << "Starting the online solve."
         << "\nThe following time settings are used:"
         << "\nExport fields every: " << exportEvery
         << "\nStore coefficients every: " << storeEvery
         << "\ndt used: " << dt << endl;
    Info << "The following number of modes were requested:"
         << "\nVelocity modes (Nphi_u): " << Nphi_u
         << "\nPressure modes (Nphi_prgh): " << Nphi_prgh
         << "\nTemperature modes (Nphi_t): " << Nphi_t << endl;
    Info << "The interpolation dimension is: " << dimA << endl;
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

// ************************************************************************ //
