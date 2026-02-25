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
Class
    ReducedUnsteadyBBTurb
Description
    A reduced problem for the unsteady, turbulent Boussinesq equations
SourceFiles
    ReducedUnsteadyBBTurb.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedUnsteadyBBTurb class.
/// \dir
/// Directory containing the header and source files for the ReducedUnsteadyBBTurb class

#ifndef ReducedUnsteadyBBTurbTurb_H
#define ReducedUnsteadyBBTurbTurb_H

#include "fvCFD.H"
#include "IOmanip.H"
#include "ReducedUnsteadyNS.H"
#include "UnsteadyBBTurb.H"
#include "BoundaryConditions.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <memory>

/// Newton object for the resolution of the reduced problem using a supremizer approach
struct NewtonUnsteadyBBTurbSup : public newton_argument<double> {
public:
    NewtonUnsteadyBBTurbSup() { }

    NewtonUnsteadyBBTurbSup(int Nx, int Ny,
        UnsteadyBBTurb& problem): newton_argument<double>(Nx, Ny),
                                  problem(&problem),
                                  Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
                                  Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
                                  N_BC_t(problem.inletIndexT.rows()),
                                  N_BC(problem.inletIndex.rows()),
                                  Nphi_prgh(problem.NPrghmodes),
                                  Nphi_nut(problem.NNutModes),
                                  nu_fluct(problem.NNutModes),
                                  nu_param(problem.mu.rows()),
                                  pCommonMatrices(problem.pCommonMatrices),
                                  pSupremizerMatrices(problem.pSupremizerMatrices),
                                  pPPEMatrices(problem.pPPEMatrices),
                                  pPenaltyMatrices(problem.pPenaltyMatrices)
    { }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;
    UnsteadyBBTurb* problem;
    std::shared_ptr<UnsteadyBBTurb::CommonMatrices> pCommonMatrices;
    std::shared_ptr<UnsteadyBBTurb::SupremizerMatrices> pSupremizerMatrices;
    std::shared_ptr<UnsteadyBBTurb::PPEMatrices> pPPEMatrices;
    std::shared_ptr<UnsteadyBBTurb::PenaltyMatrices> pPenaltyMatrices;
    int Nphi_u;
    int Nphi_t;
    int N_BC_t;
    int N_BC;
    int Nphi_prgh;
    int Nphi_nut;
    Eigen::VectorXd nu_fluct;
    Eigen::VectorXd nu_param;
    scalar nu;
    scalar Pr;
    scalar Pr_t;
    scalar dt;

    Eigen::VectorXd y_old;
    Eigen::VectorXd BC_t;
    Eigen::VectorXd BC;
    Eigen::MatrixXd tauT;
    Eigen::MatrixXd tauU;
    std::vector<SPLINTER::RBFSpline*> SPLINES;
};

struct NewtonUnsteadyBBTurbSupMomentum : public newton_argument<double> {
public:
    NewtonUnsteadyBBTurbSupMomentum() { }

    NewtonUnsteadyBBTurbSupMomentum(int Nx, int Ny,
        UnsteadyBBTurb& problem): newton_argument<double>(Nx, Ny),
                                  problem(&problem),
                                  Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
                                  Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
                                  N_BC_t(problem.inletIndexT.rows()),
                                  N_BC(problem.inletIndex.rows()),
                                  Nphi_prgh(problem.NPrghmodes),
                                  Nphi_nut(problem.NNutModes),
                                  nu_fluct(problem.NNutModes),
                                  nu_param(problem.mu.rows()),
                                  pCommonMatrices(problem.pCommonMatrices),
                                  pSupremizerMatrices(problem.pSupremizerMatrices),
                                  pPPEMatrices(problem.pPPEMatrices),
                                  pPenaltyMatrices(problem.pPenaltyMatrices)
    { }
    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;
    UnsteadyBBTurb* problem;
    std::shared_ptr<UnsteadyBBTurb::CommonMatrices> pCommonMatrices;
    std::shared_ptr<UnsteadyBBTurb::SupremizerMatrices> pSupremizerMatrices;
    std::shared_ptr<UnsteadyBBTurb::PPEMatrices> pPPEMatrices;
    std::shared_ptr<UnsteadyBBTurb::PenaltyMatrices> pPenaltyMatrices;
    int Nphi_u;
    int Nphi_t;
    int N_BC_t;
    int N_BC;
    int Nphi_prgh;
    int Nphi_nut;
    Eigen::VectorXd nu_fluct;
    Eigen::VectorXd nu_param;
    scalar nu;
    scalar Pr;
    scalar Pr_t;
    scalar dt;

    Eigen::VectorXd y_old;
    Eigen::VectorXd BC_t;
    Eigen::VectorXd BC;
    Eigen::MatrixXd tauT;
    Eigen::MatrixXd tauU;
    std::vector<SPLINTER::RBFSpline*> SPLINES;
};

struct newtonUnsteadyBBTurbSUPTemperature : public newton_argument<double> {
public:
    newtonUnsteadyBBTurbSUPTemperature() { }

    newtonUnsteadyBBTurbSUPTemperature(int Nx, int Ny,
        UnsteadyBBTurb& problem): newton_argument<double>(Nx, Ny),
                                  problem(&problem),
                                  Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
                                  Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
                                  N_BC_t(problem.inletIndexT.rows()),
                                  N_BC(problem.inletIndex.rows()),
                                  Nphi_prgh(problem.NPrghmodes),
                                  Nphi_nut(problem.NNutModes),
                                  nu_fluct(problem.NNutModes),
                                  nu_param(problem.mu.rows()),
                                  pCommonMatrices(problem.pCommonMatrices),
                                  pSupremizerMatrices(problem.pSupremizerMatrices),
                                  pPPEMatrices(problem.pPPEMatrices),
                                  pPenaltyMatrices(problem.pPenaltyMatrices)
    { }
    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;
    UnsteadyBBTurb* problem;
    std::shared_ptr<UnsteadyBBTurb::CommonMatrices> pCommonMatrices;
    std::shared_ptr<UnsteadyBBTurb::SupremizerMatrices> pSupremizerMatrices;
    std::shared_ptr<UnsteadyBBTurb::PPEMatrices> pPPEMatrices;
    std::shared_ptr<UnsteadyBBTurb::PenaltyMatrices> pPenaltyMatrices;
    int Nphi_u;
    int Nphi_t;
    int N_BC_t;
    int N_BC;
    int Nphi_prgh;
    int Nphi_nut;
    Eigen::VectorXd nu_fluct;
    Eigen::VectorXd nu_param;
    scalar nu;
    scalar Pr;
    scalar Pr_t;
    scalar dt;

    Eigen::VectorXd y_old;
    Eigen::VectorXd BC_t;
    Eigen::VectorXd BC;
    Eigen::MatrixXd tauT;
    Eigen::MatrixXd tauU;
    std::vector<SPLINTER::RBFSpline*> SPLINES;
};

struct NewtonUnsteadyBBTurbVMB : public newton_argument<double> {
public:
    NewtonUnsteadyBBTurbVMB() { }

    NewtonUnsteadyBBTurbVMB(int Nx, int Ny,
        UnsteadyBBTurb& problem): newton_argument<double>(Nx, Ny),
                                  problem(&problem),
                                  Nphi_u(problem.NUmodes + problem.liftfield.size()),
                                  Nphi_prgh(problem.NPrghmodes),
                                  Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
                                  N_BC_t(problem.inletIndexT.rows()),
                                  N_BC(problem.inletIndex.rows()),
                                  Nphi_nut(problem.NNutModes),
                                  nu_fluct(problem.NNutModes),
                                  nu_param(problem.mu.rows()),
                                  pCommonMatrices(problem.pCommonMatrices),
                                  pSupremizerMatrices(problem.pSupremizerMatrices),
                                  pPPEMatrices(problem.pPPEMatrices),
                                  pPenaltyMatrices(problem.pPenaltyMatrices)
    { }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;
    UnsteadyBBTurb* problem;
    std::shared_ptr<UnsteadyBBTurb::CommonMatrices> pCommonMatrices;
    std::shared_ptr<UnsteadyBBTurb::SupremizerMatrices> pSupremizerMatrices;
    std::shared_ptr<UnsteadyBBTurb::PPEMatrices> pPPEMatrices;
    std::shared_ptr<UnsteadyBBTurb::PenaltyMatrices> pPenaltyMatrices;
    int Nphi_u;
    int Nphi_t;
    int Nphi_prgh;
    int N_BC_t;
    int N_BC;
    int Nphi_nut;
    Eigen::VectorXd nu_fluct;
    Eigen::VectorXd nu_param;
    scalar nu;
    scalar Pr;
    scalar Pr_t;
    scalar dt;

    Eigen::VectorXd y_old;
    Eigen::VectorXd BC_t;
    Eigen::VectorXd BC;
    Eigen::MatrixXd tauT;
    Eigen::MatrixXd tauU;
    std::vector<SPLINTER::RBFSpline*> SPLINES;
};

struct NewtonUnsteadyBBTurbPPE : public newton_argument<double> {
public:
    NewtonUnsteadyBBTurbPPE() { }

    NewtonUnsteadyBBTurbPPE(int Nx, int Ny,
        UnsteadyBBTurb& problem): newton_argument<double>(Nx, Ny),
                                  problem(&problem),
                                  Nphi_u(problem.NUmodes + problem.liftfield.size()),
                                  Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
                                  N_BC_t(problem.inletIndexT.rows()),
                                  N_BC(problem.inletIndex.rows()),
                                  Nphi_prgh(problem.NPrghmodes),
                                  Nphi_nut(problem.NNutModes),
                                  nu_fluct(problem.NNutModes),
                                  nu_param(problem.mu.rows()),
                                  pCommonMatrices(problem.pCommonMatrices),
                                  pSupremizerMatrices(problem.pSupremizerMatrices),
                                  pPPEMatrices(problem.pPPEMatrices),
                                  pPenaltyMatrices(problem.pPenaltyMatrices)
    { }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;

    UnsteadyBBTurb* problem;
    std::shared_ptr<UnsteadyBBTurb::CommonMatrices> pCommonMatrices;
    std::shared_ptr<UnsteadyBBTurb::SupremizerMatrices> pSupremizerMatrices;
    std::shared_ptr<UnsteadyBBTurb::PPEMatrices> pPPEMatrices;
    std::shared_ptr<UnsteadyBBTurb::PenaltyMatrices> pPenaltyMatrices;
    int Nphi_u;
    int Nphi_t;
    int Nphi_prgh;
    int Nphi_nut;
    int N_BC_t;
    int N_BC;
    Eigen::VectorXd nu_fluct;
    Eigen::VectorXd nu_param;
    scalar nu;
    scalar Pr;
    scalar Pr_t;
    scalar dt;

    Eigen::VectorXd y_old;
    Eigen::VectorXd BC_t;
    Eigen::VectorXd BC;
    Eigen::MatrixXd tauT;
    Eigen::MatrixXd tauU;
    std::vector<SPLINTER::RBFSpline*> SPLINES;
};

/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where a reduced problem for the <b> unsteady, turbulent, buoyant Navier-stokes </b> problem is implemented.
/** In this class are implemented the methods to solve a reduced version of the unsteady, turbulent, buoyant Navier-stokes equations. This class must be used together with the FOM class */
class ReducedUnsteadyBBTurb : public reducedUnsteadyNS
{
private:
    /// Word to indicate the stabilization method. It gets read from the ITHACAdict, and can be either "supremizer", "PPE" or "VMB"
    word method;

    /// Boolean flag to indicate wheter the penalty factor will be optimized using the Star method
    bool optimizePenaltyFactor = false;

    /// Number of timesteps over which the penalty method is optimized (if any)
    int nTimestepsPenaltyAdapt = 0;

    /// Maximum number of inner iteration for the adaptive penalty method (if any)
    int penaltyMaxIter = 0;

    /// Tolerance for the adaptive penalty method (if any)
    float penaltyTolU = 1e-3;

    /// Tolerance for the adaptive penalty method (if any)
    float penaltyTolT = 1e-2;

    // Index of the first RBF mode in the reduced solution vector, used for the online interpolation of the eddy viscosity field. Used to ignore the lifting modes in the interpolation
    int firstRBFIndex = 0;

    // Index of the current online solution being solved
    int onlineSolutionIndex = 0;

    /// @brief   Carries out initial checks and assertions
    ///
    void validateSettings();

    /// @brief    Inverse Distance Weighting function to interpolate the averaged eddy viscosity coefficients
    ///           IDW in this case is a misnomer, as the implementation is specific for a vector of 0s except for the i-th component
    /// @param[in]  mu_now   The vector of the online parameter
    ///
    Eigen::VectorXd interpolateIDW();

    /// @brief      Carries out the time loop for the online solve. Templated to support different Newton object types.
    ///
    /// @param[in,out]  newtonObject      The Newton object to be used in the time loop
    /// @param[in,out]  boundaryConditions  The BoundaryConditions object containing the BCs
    /// @param[in]      firstRBFInd      The first RBF index to use
    /// @param[in]      numberOfStores   The number of stores to be used in the time loop
    template <typename NewtonType>
    void onlineTimeLoop(NewtonType& newtonObject, BoundaryConditions& boundaryConditions, int firstRBFInd, int numberOfStores);

    /// @brief      Carries out the time loop for the online solve, but uses a newton object for momentum and one for temperature. Templated to support different Newton object types.
    ///
    /// @param[in,out]  newtonObject      The Newton object to be used in the time loop
    /// @param[in,out]  boundaryConditions  The BoundaryConditions object containing the BCs
    /// @param[in]      firstRBFInd      The first RBF index to use
    /// @param[in]      numberOfStores   The number of stores to be used in the time loop
    template <typename NewtonTypeMomentum, typename NewtonTypeTemperature>
    void onlineTimeLoopSegregated(NewtonTypeMomentum& newtonObjectMomentum, NewtonTypeTemperature& newtonObjectTemperature, BoundaryConditions& boundaryConditions, int firstRBFInd, int numberOfStores);

    /// @brief      Set the boundary conditions in the (templated) Newton object
    ///
    /// @param[in,out]  newtonObject  The Newton object where to set the BCs
    /// @param[in]      boundaryConditions  The BoundaryConditions object containing the BCs
    template <typename NewtonType>
    void setNewtonObjectBC(NewtonType& newtonObject, BoundaryConditions& boundaryConditions);

    /// @brief      Method to configure a generic Newton object with the current settings. What it does is assign physical properties,
    ///             time step, boundary conditions, stabilization coefficients, etc. Templated to support different Newton object types.
    ///
    /// @param[in,out]  newtonObject  The Newton object to be configured
    /// @param[in]      currentY     The current reduced solution vector
    template <typename NewtonType>
    void configureNewtonObject(NewtonType& newtonObject, const Eigen::VectorXd& currentY, BoundaryConditions& boundaryConditions);

    /// @brief Method to perform an online solve using a supremizer stabilisation method
    /// @param[in]  boundaryConditions The boundary conditions object to be used in the online solve
    /// @param[in]  startSnap Index of the snapshot used to compute the initial coefficients
    ///
    void solveOnline_sup(BoundaryConditions& boundaryConditions, int startSnap = 0);

    /// @brief Method to perform an online solve using a PPE stabilisation method
    /// @param[in]  boundaryConditions The boundary conditions object to be used in the online solve
    /// @param[in]  startSnap Index of the snapshot used to compute the initial coefficients
    ///
    void solveOnline_PPE(BoundaryConditions& boundaryConditions, int startSnap = 0);

    /// @brief Method to perform an online solve using a VMB  method
    ///
    /// @param[in]  boundaryConditions The boundary conditions object to be used in the online solve
    /// @param[in]  startSnap Index of the snapshot used to compute the initial coefficients
    ///
    void solveOnline_VMB(BoundaryConditions& boundaryConditions, int startSnap = 0);

    /// @brief   Estimate the penalty factor using the iterative approach developed by S. Star.
    ///
    /// @param[in] boundaryConditions  BoundaryConditions object containing the BCs to be used in the online solve
    /// @param[in] startSnap        Snapshot to use as initial condition
    ///
    void estimatePenaltyFactorSupremizer(BoundaryConditions& boundaryConditions, int startSnap);

    /// @brief   Estimate the penalty factor using the iterative approach developed by S. Star.
    ///
    /// @param[in] boundaryConditions  BoundaryConditions object containing the BCs to be used in the online solve
    /// @param[in] startSnap        Snapshot to use as initial condition
    ///
    void estimatePenaltyFactorPPE(BoundaryConditions& boundaryConditions, int startSnap);

public:
    // Constructors
    /// Construct Null
    ReducedUnsteadyBBTurb();
    /// Construct Null
    ///
    /// @param[in]  problem  a full order UnsteadyBBTurb problem
    ///
    explicit ReducedUnsteadyBBTurb(UnsteadyBBTurb& problem);
    ~ReducedUnsteadyBBTurb() { };

    /// List of pointers to store the modes for temperature including lift field modes
    PtrList<volScalarModes> LTmodes;

    /// List of pointers to store the modes for velocity including lift field modes
    PtrList<volVectorModes> LUmodes;

    /// List of pointers to store the eddy viscosity modes
    PtrList<volScalarModes> nutmodes;

    /// List of pointers to store the shifted pressure (prgh) modes
    PtrList<volScalarModes> Prghmodes;

    /// Function object to call the non linear solver sup approach
    NewtonUnsteadyBBTurbSup newton_object_sup;

    /// Newton object to call the non linear solver for momentum when solving momentum and temperature equations separately
    NewtonUnsteadyBBTurbSupMomentum newton_object_sup_momentum;

    /// Newton object to call the non linear solver for temperature when solving momentum and temperature equations separately
    newtonUnsteadyBBTurbSUPTemperature newton_object_sup_temperature;

    /// Function object to call the non linear solver PPE approach
    NewtonUnsteadyBBTurbPPE newton_object_PPE;

    /// Function object to call the non linear solver VMB approach
    NewtonUnsteadyBBTurbVMB newton_object_VMB;

    /// Shared pointer to common matrices object for the reduced problem
    std::shared_ptr<UnsteadyBBTurb::CommonMatrices> pCommonMatrices;

    /// Shared pointer to supremizer matrices object for the reduced problem
    std::shared_ptr<UnsteadyBBTurb::SupremizerMatrices> pSupremizerMatrices;

    /// Shared pointer to PPE matrices object for the reduced problem
    std::shared_ptr<UnsteadyBBTurb::PPEMatrices> pPPEMatrices;

    /// Shared pointer to penalty matrices object for the reduced problem
    std::shared_ptr<UnsteadyBBTurb::PenaltyMatrices> pPenaltyMatrices;

    /// Reconstructed temperature field
    PtrList<volScalarField> TRecFields;

    /// Reconstructed pressure field
    PtrList<volScalarField> prghRecFields;

    /// Reconstructed turbulent viscosity field
    PtrList<volScalarField> nutFluctRecFields;

    /// Reconstructed total turbulent viscosity field
    PtrList<volScalarField> nutRecFields;

    /// Initial eddy viscosity reduced coefficients - fluctuating part
    Eigen::VectorXd nut0;

    /// Initial eddy viscosity reduced coefficients - parametrized part
    Eigen::VectorXd nut_param_0;

    // Parameters for the current run:
    Eigen::VectorXd mu_now;

    /// Number of temperature modes
    int Nphi_t;

    /// Number of shifted pressure modes
    int Nphi_prgh;

    /// Number of turbulent viscosity mdodes
    int Nphi_nut;

    /// Number of parametrized boundary conditions related to temperature field
    int N_BC_t;

    /// Dimensions of the interpolation independent variable (for example the (A, A_dot) pair)
    int dimA;

    /// DimensionedScalar Pr;
    scalar Pr;

    /// DimensionedScalar Prt;
    scalar Pr_t;

    /// Matrix of eddy viscosity RBF coefficients, has the first row as time values
    Eigen::MatrixXd rbfCoeffMat;

    /// Penalty term for the temperature equation
    Eigen::MatrixXd penaltyT;

    /// Penalty coefficient tauT for the temperature equation
    Eigen::MatrixXd tauT;

    /// Boolean to choose whether to skip the lifting modes in the reconstruction
    bool skipLift = true;

    /// Boolean flag to chose whether to solve the momentum and temperature equations separately using two different Newton solvers
    bool separateMomentumTemperature = false;

    /// Pointer to the FOM problem
    UnsteadyBBTurb* problem;

    // Methods

    ///@brief       Computes the reduced version of the inf_sup constant
    ///
    void inf_sup_constant();

    /// @brief    Method to set the time setting for the online solution
    ///
    /// @param[in]  timeMatrix  The matrix containing the time settings, with the following columns:
    ///                         Column 1: initial time
    ///                         Column 2: final time
    ///                         Column 3: time step
    ///                         Column 4: exportEvery
    ///                         Column 5: storeEvery
    ///
    void setTimeSettings(const Eigen::MatrixXd& timeMatrix);

    /// @brief    Method to reconstruct a solution from an online solve
    ///
    /// @param[in]  exportFields  Boolean to choose whether to export or not the reconstructed fields
    /// @param[in]  folder        The folder where you want to store the results
    ///
    void reconstructSolution(bool exportFields = false, fileName folder = "./online_rec");

    /// @brief    High level method to perform an online solve. It calls the specific method for the chosen stabilization method
    ///
    /// @param[in]  temp_now_BC   The scalar of online temperature. It is defined in
    /// with an Eigen::MatrixXd. For the time-dependent BCs, it has the first column as time,
    /// then each nth column corresponds to the value of the
    /// n-th BC at that time.
    /// @param[in] vel_now_BC   The scalar of online velocity. It is defined in
    /// with an Eigen::MatrixXd. For the time-dependent BCs, it has the first column as time,
    /// then each nth column corresponds to the value of the n-th BC at that time.
    /// @param[in]  startSnap The first snapshot taken from the offline snapshots and used to get the reduced initial condition default = 0.
    ///
    void solveOnline(const Eigen::MatrixXd& vel_now_BC, const Eigen::MatrixXd& temp_now_BC, int startSnap = 0);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
