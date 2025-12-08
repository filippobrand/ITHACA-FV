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
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Newton object for the resolution of the reduced problem using a supremizer approach
struct newton_unsteadyBBTurb_sup : public newton_argument<double> {
public:
    newton_unsteadyBBTurb_sup() { }

    newton_unsteadyBBTurb_sup(int Nx, int Ny,
        UnsteadyBBTurb& problem): newton_argument<double>(Nx, Ny),
                                  problem(&problem),
                                  Nphi_u(problem.NUmodes + problem.liftfield.size() + problem.NSUPmodes),
                                  Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
                                  N_BC_t(problem.inletIndexT.rows()),
                                  N_BC(problem.inletIndex.rows()),
                                  Nphi_prgh(problem.NPrghmodes),
                                  Nphi_nut(problem.Nnutmodes),
                                  nu_fluct(problem.Nnutmodes),
                                  nu_param(problem.mu.rows())
    { }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;
    UnsteadyBBTurb* problem;
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
    std::vector<SPLINTER::RBFSpline*> SPLINES;
};

struct newton_unsteadyBBTurb_PPE : public newton_argument<double> {
public:
    newton_unsteadyBBTurb_PPE() { }

    newton_unsteadyBBTurb_PPE(int Nx, int Ny,
        UnsteadyBBTurb& problem): newton_argument<double>(Nx, Ny),
                                  problem(&problem),
                                  Nphi_u(problem.NUmodes + problem.liftfield.size()),
                                  Nphi_t(problem.NTmodes + problem.liftfieldT.size()),
                                  N_BC_t(problem.inletIndexT.rows()),
                                  N_BC(problem.inletIndex.rows()),
                                  Nphi_p(problem.NPmodes),
                                  Nphi_prgh(problem.NPrghmodes),
                                  Nphi_nut(problem.Nnutmodes)
    { }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;

    UnsteadyBBTurb* problem;
    int Nphi_u;
    int Nphi_p;
    int Nphi_t;
    int Nphi_prgh;
    int Nphi_nut;
    int N_BC_t;
    int N_BC;
    scalar nu;
    scalar Pr;
    scalar Pr_t;
    scalar dt;
    int writeEvery;
    word BCmethod;

    Eigen::VectorXd y_old;
    Eigen::VectorXd BC_t;
    Eigen::VectorXd BC;
    std::vector<SPLINTER::RBFSpline*> SPLINES;
};

/*---------------------------------------------------------------------------*\
                        Class reducedProblem Declaration
\*---------------------------------------------------------------------------*/

/// Class where a reduced problem for the <b> unsteady, turbulent, buoyant Navier-stokes </b> problem is implemented.
/** In this class are implemented the methods to solve a reduced version of the stationary Navier-stokes
equations. This class must be used together with the steadyNS class */
class ReducedUnsteadyBBTurb : public reducedUnsteadyNS
{
private:

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
    PtrList<volScalarField> LTmodes;

    /// List of pointers to store the modes for velocity including lift field modes
    PtrList<volVectorField> LUmodes;

    /// List of pointers to store the eddy viscosity modes
    PtrList<volScalarField> nutmodes;

    /// List of pointers to store the shifted pressure (prgh) modes 
    PtrList<volScalarField> Prghmodes;

    /// Function object to call the non linear solver sup approach
    newton_unsteadyBBTurb_sup newton_object_sup;

    /// Function object to call the non linear solver PPE approach
    newton_unsteadyBBTurb_PPE newton_object_PPE;

    /// Reconstructed temperature field
    PtrList<volScalarField> TRecFields;

    /// Reconstructed pressure field
    PtrList<volScalarField> prghRecFields;

    /// Reconstructed turbulent viscosity field
    PtrList<volScalarField> nutFluctRecFields;

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

    /// DimensionedScalar Pr;
    scalar Pr;

    /// DimensionedScalar Prt;
    scalar Pr_t;

    /// Dimensions of the interpolation independent variable (for example the (A, A_dot) pair)
    int dimA;

    /// Matrix of eddy viscosity RBF coefficients, has the first row as time values
    Eigen::MatrixXd rbfCoeffMat;

    /// Boolean to choose whether to skip the lifting modes in the reconstruction
    bool skipLift = false;

    /// Pointer to the FOM problem
    UnsteadyBBTurb* problem;

    // Functions

    /// Method to perform an online solve using a supremizer stabilisation method
    ///
    /// @param[in]  temp_now_BC   The scalar of online temperature. It is defined in
    /// with an Eigen::MatrixXd and must have one col and as many rows as the number
    /// of parametrized boundary conditions.
    /// @param[in] vel_now_BC   The scalar of online velocity. It is defined in
    /// with an Eigen::MatrixXd and must have one col and as many rows as the number
    /// of parametrized boundary conditions.
    /// @param[in]  NParaSet The parameter sets number for which the online solution is solved.
    /// @param[in]  startSnap The first snapshot taken from the offline snapshots
    /// and used to get the reduced initial condition default = 0.
    ///
    void solveOnline_sup(Eigen::MatrixXd& temp_now_BC,
        Eigen::MatrixXd& vel_now_BC, int NParaSet = 0, int startSnap = 0);

    /// @brief    Method to reconstruct a solution from an online solve
    ///
    /// @param[in]  exportFields  Boolean to choose whether to export or not the reconstructed fields
    /// @param[in]  folder        The folder where you want to store the results
    ///
    void reconstructSolution(bool exportFields = false, fileName folder = "./online_rec");

    /// @brief    Inverse Distance Weighting function to interpolate the averaged eddy viscosity coefficients
    ///           IDW in this case is a misnomer, as the implementation is specific for a vector of 0s except for the i-th component
    /// @param[in]  mu_now   The vector of the online parameter
    ///
    Eigen::VectorXd interpolateIDW();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
