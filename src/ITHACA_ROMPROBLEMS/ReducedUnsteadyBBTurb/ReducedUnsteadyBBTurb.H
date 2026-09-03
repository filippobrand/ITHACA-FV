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
#include "TimeManager.H"
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <memory>

/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

struct PenaltySettings
{
    bool optimizePenaltyFactor = false;
    int nTimestepsPenaltyAdapt = 0;
    int penaltyMaxIter = 0;
    float penaltyTolU = 1e-3;
    float penaltyTolT = 1e-2;
};

struct ROMSettings
{
    word method = "PPE";
    word bcMethod = "Gunzburger";
    bool hasMonitors = false;
    bool timeDependentBC = false;
    bool normalizeWithEigenvalues = false;
    int N_BC = 0;
    int N_BC_t = 0;
    word solverODE = "BDF2";
};

struct InterpolationSettings
{
    int firstRBFIndex = 0;
    bool skipLift = true;
    bool derivative = false;
    int dimInputRBF = 0;
    int avgTermOfflineSamples = 0;
};

struct NumberOfModes
{
    int velocity = 0;
    int pressure = 0;
    int temperature = 0;
    int nut = 0;
};

struct FluidProperties
{
    scalar nu = 0;
    scalar Pr = 0;
    scalar Pr_t = 0;
};

// Forward declaration of the ReducedUnsteadyBBTurbSystem class, which will be defined later
class ReducedUnsteadyBBTurbSystem;

/// Class where a reduced problem for the <b> unsteady, turbulent, buoyant Navier-stokes </b> problem is implemented.
/** In this class are implemented the methods to solve a reduced version of the unsteady, turbulent, buoyant Navier-stokes equations. This class must be used together with the FOM class */
class ReducedUnsteadyBBTurb : public reducedUnsteadyNS
{
        friend class ReducedUnsteadyBBTurbSystem;
        friend class SystemPPEGunzburger;
        friend class SystemSupremizerGunzburger;
    private:
        // General settings
        ROMSettings romSettings_;

        // Settings about the penalty method
        PenaltySettings penaltySettings_;

        // Settings about the interpolation method
        InterpolationSettings interpolationSettings_;

        // Number of trial modes (Might be different from the POD modes)
        NumberOfModes testModes_;

        // Number of modes used to express the solution. (Is the POD modes)
        NumberOfModes expressionModes_;

        // Boundary conditions object, used to handle all things BC, including lift and time-dep BC
        BoundaryConditions boundaryConditions_;

        /// Initial eddy viscosity reduced coefficients - fluctuating part
        Eigen::VectorXd nutCurrentCoeffs_;

        /// Initial eddy viscosity reduced coefficients - parametrized part
        Eigen::VectorXd nutAvgCurrentCoeffs_;

        /// Scaling vector to handle normalization
        Eigen::VectorXd scalingVector_;

        /// Velocity POD eigenvalues
        Eigen::VectorXd uEigenvalues_;

        /// Pressure POD eigenvalues
        Eigen::VectorXd pEigenvalues_;

        /// Temperature POD eigenvalues
        Eigen::VectorXd tEigenvalues_;

        /// Eddy viscosity POD eigenvalues
        Eigen::VectorXd nutEigenvalues_;

        /// Solution of the online ODE system
        Eigen::MatrixXd onlineSolution;

        /// Storing nut coefficients
        Eigen::MatrixXd nutCoeffMat;

        /// Pointer to the ODE system to solve
        std::unique_ptr<ReducedUnsteadyBBTurbSystem> odeSystem_;

        /// Pointer to the ODE solver to be used
        std::unique_ptr<ODESolver> odeSolver_;

        /// @brief    Initializes the settings attributes
        ///
        void initializeSettings();
        
        /// @brief    Initializes the ODE systems
        ///
        void initializeODESystem();

        /// @brief    Initialize the ODE solver
        ///
        void initializeODESolver();

        /// @brief    Initializes the normalization with eigenvalues
        ///
        void initializeEigenvalues();
        
        /// @brief    Collect relevant matrices
        void collectMatrices();

        /// @brief    Reads the matrix dimensions
        void initializeDimensions();

        /// @brief    Does the nut interpolation
        void interpolateNutCoeffs();
        
        /// @brief    Reads the market vectors from disk and fills the eigenvalues vectors
        void readEigenvalues();

        /// @brief    Inverse Distance Weighting function to interpolate the averaged eddy viscosity coefficients
        ///           IDW in this case is a misnomer, as the implementation is specific for a vector of 0s except for the i-th component
        /// @param[in]  input_parameters   The vector of the online parameter
        ///
        Eigen::VectorXd interpolateIDW(const Eigen::VectorXd& input_parameters);

    public:
        // Constructors
        /// Construct Null
        ReducedUnsteadyBBTurb();

        /// @param[in]  problem  a full order UnsteadyBBTurb problem
        ///
        explicit ReducedUnsteadyBBTurb(UnsteadyBBTurb& problem);

        // Default destructor
        virtual ~ReducedUnsteadyBBTurb();

        // Fluid Properties
        FluidProperties fluidProperties;

        /// Shared pointer to common matrices object for the reduced problem
        std::shared_ptr<UnsteadyBBTurb::CommonMatrices> pCommonMatrices;

        /// Shared pointer to supremizer matrices object for the reduced problem
        std::shared_ptr<UnsteadyBBTurb::SupremizerMatrices> pSupremizerMatrices;

        /// Shared pointer to PPE matrices object for the reduced problem
        std::shared_ptr<UnsteadyBBTurb::PPEMatrices> pPPEMatrices;

        /// Shared pointer to penalty matrices object for the reduced problem
        std::shared_ptr<UnsteadyBBTurb::PenaltyMatrices> pPenaltyMatrices;

        /// Shared pointer to the Gunzburger matrices object for the reduced problem
        std::shared_ptr<UnsteadyBBTurb::GunzburgerMatrices> pGunzburgerMatrices;

        /// Reconstructed temperature field
        PtrList<volScalarField> TRecFields;

        /// Reconstructed pressure field
        PtrList<volScalarField> prghRecFields;

        /// Reconstructed turbulent viscosity field
        PtrList<volScalarField> nutFluctRecFields;

        /// Reconstructed total turbulent viscosity field
        PtrList<volScalarField> nutRecFields;

        /// Matrix of eddy viscosity RBF coefficients, has the first row as time values
        Eigen::MatrixXd rbfCoeffMat;

        /// Penalty related matrices
        Eigen::MatrixXd tauU, tauT;

        /// Pointer to the FOM problem
        UnsteadyBBTurb* problem;

        // Methods
        ///@brief       Computes the reduced version of the inf_sup constant
        ///
        void inf_sup_constant();

        /// @brief    Sets the fluid properties for the reduced problem
        ///
        /// @param[in]  nu  Kinematic viscosity
        /// @param[in]  Pr  Prandtl number
        /// @param[in]  Pr_t  Turbulent Prandtl number
        ///
        void setFluidProperties(scalar nu, scalar Pr, scalar Pr_t);

        /// @brief    Method that, given a boundary condition object, solves the online problem
        ///
        void solveOnline(const Eigen::MatrixXd& vel_now_BC,
                         const Eigen::MatrixXd& temp_now_BC, int startSnap, TimeManager& timeManager);

        /// @brief    Method to reconstruct a solution from an online solve
        ///
        /// @param[in]  exportFields  Boolean to choose whether to export or not the reconstructed fields
        /// @param[in]  folder        The folder where you want to store the results
        ///
        void reconstructSolution(TimeManager& time_manager, bool exportFields = false,
                                 fileName folder = "./online_rec");

        /// @brief   Saves the reduced coefficients and the RBF ones to a specified folder
        ///
        /// @param[in]  folder  The folder where to save the results
        /// @param[in]  modeIdentifier  The identifier for the mode numbers. Used when comparing different runs
        void saveCoefficients(word folder = "./reducedCoeffs",
                              word modeIdentifier = "default");

        /// @brief   Const method that returns a reference to the boundary conditions object, to be used in the ODE
        const BoundaryConditions& getBoundaryConditions() const
        {
            return boundaryConditions_;
        }

        /// @brief Reset
        ///
        void reset();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
