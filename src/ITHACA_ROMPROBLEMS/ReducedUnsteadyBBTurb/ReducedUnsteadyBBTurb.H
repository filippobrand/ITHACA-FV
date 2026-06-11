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
#include "ReducedODESystem.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <memory>

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

    /// Word to select the ODE solver
    word odeSolverName;

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

    // Residual vector
    Eigen::VectorXd residual;

    /// @brief   Carries out initial checks and assertions
    ///
    void validateSettings();
    
    ///@brief       Computes the reduced version of the inf_sup constant
    ///
    void inf_sup_constant();

    /// @brief      Carries out the time loop for the online solve. Templated to support different Newton object types.
    ///
    /// @param[in]      numberOfStores   The number of stores to be used in the time loop
    template <typename ODESolverType>
    void onlineTimeLoop(ODESolverType& ode_solver, int numberOfStores);

    /// @brief    Inverse Distance Weighting function to interpolate the averaged eddy viscosity coefficients
    ///           IDW in this case is a misnomer, as the implementation is specific for a vector of 0s except for the i-th component
    /// @param[in]  mu_now   The vector of the online parameter
    ///
    Eigen::VectorXd interpolateIDW();

    /// @brief Method to perform an online solve using a PPE stabilisation method
    ///
    void solveOnline_PPE();

    /// @brief Method to carry out the RBF interpolation of nut
    ///
    template <typename ODESolverType>
    void interpolateEddyViscosity(ODESolverType& ode_solver);

public:
    // Constructors
    /// Construct Null
    ReducedUnsteadyBBTurb();
    /// Construct Null
    ///
    /// @param[in]  problem  a full order UnsteadyBBTurb problem
    ///
    explicit ReducedUnsteadyBBTurb(UnsteadyBBTurb& problem);

    /// Destructor
    ~ReducedUnsteadyBBTurb() { };

    /// Boundary Conditions object
    BoundaryConditions boundaryConditions;

    /// List of pointers to store the modes for temperature including lift field modes
    PtrList<volScalarModes> LTmodes;

    /// List of pointers to store the modes for velocity including lift field modes
    PtrList<volVectorModes> LUmodes;

    /// List of pointers to store the eddy viscosity modes
    PtrList<volScalarModes> nutmodes;

    /// List of pointers to store the shifted pressure (prgh) modes
    PtrList<volScalarModes> Prghmodes;

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

    /// Vector of eddy viscosity coefficients for the fluctuating part
    Eigen::VectorXd nut_fluct;

    /// Vector of eddy viscosity coefficients for the parametrized part
    Eigen::VectorXd nut_param;

    // Parameters for the current run:
    Eigen::VectorXd mu_now;

    /// Number of temperature modes
    int Nphi_t;

    /// Number of shifted pressure modes
    int Nphi_prgh;

    /// Number of turbulent viscosity mdodes
    int Nphi_nut;

    // Velocity test function space size
    int Ntest_u;

    // Pressure test function space size
    int Ntest_prgh;

    // Temperature test function space size
    int Ntest_t;

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

      /// Bool to use lifting
    bool liftBCs = false;
    
    /// Bool to use gunzburger
    bool gunzburgerBCs = false;

    /// Bool to use the penalty method
    bool penaltyBCs = false;

    /// Boolean flag to indicate wheter the penalty factor will be optimized using the Star method
    bool optimizePenaltyFactor = false;

    /// Boolean to indicate whther time dep BCs are used
    bool timeDepBCs = false;

    /// Pointer to the FOM problem
    UnsteadyBBTurb* problem;

    // Methods
    /// @brief    Method to set the time setting for the online solution
    ///
    /// @param[in]  timeMatrix  The matrix containing the time settings, with the following columns:
    ///                         Column 1: initial time
    ///                         Column 2: final time
    ///                         Column 3: time step
    ///                         Column 4: exportEvery
    ///                         Column 5: storeEvery
    ///
    void setTimeSettings(const Eigen::MatrixXd& timeMatrix, int index);

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

    /// @brief   Saves the reduced coefficients and the RBF ones to a specified folder
    ///
    /// @param[in]  folder  The folder where to save the results
    void saveCoefficients(word folder = "./reducedCoeffs");

    /// @brief   Resets the online solution and the reconstructed fields, to be ready for a new online solve
    void reset();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*---------------------------------------------------------------------------*\
                        Class ODESystemAdapter Declaration
\*---------------------------------------------------------------------------*/
/// Adapter class that inherits from ReducedODESystem and defines a evaluateResidual method to interface with ODE solvers

class ODEStructurePPETurb : public ReducedODESystem
{
  private:
    const ReducedUnsteadyBBTurb& rom;

    bool use_penalty_bc = rom.penaltyBCs;
    bool use_time_dep_bc = rom.timeDepBCs;
    bool use_gunz_bc = rom.gunzburgerBCs;
    bool use_lift_bc = rom.liftBCs;
  
  public:
    explicit ODEStructurePPETurb(const ReducedUnsteadyBBTurb& rom) : rom(rom) { }

    void evaluateRHS(const Eigen::VectorXd& state, Eigen::VectorXd& rhs, double t) const override
    {
        // This method is not used in the current implementation, as the residual is evaluated directly in the evaluateResidual method. However, it must be defined to inherit from ReducedODESystem.
    }

    void evaluateResidual(const Eigen::VectorXd& state, const Eigen::VectorXd& state_dot, Eigen::VectorXd& residual, double t) const override;
};

#endif
