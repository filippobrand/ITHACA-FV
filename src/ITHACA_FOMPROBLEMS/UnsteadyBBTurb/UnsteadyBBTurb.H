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
    UnsteadyBBTurb
Description
    Reduction class for a Stationary Navier-Stokes problem.
SourceFiles
    UnsteadyBBTurb.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the UnsteadyBBTurb class.
/// \dir
/// Directory containing the header and source files for the UnsteadyBBTurb class.

#ifndef UnsteadyBBTurb_H
#define UnsteadyBBTurb_H
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "UnsteadyBB.H"
#include <iostream>


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                        Class SteadyBB Declaration
\*---------------------------------------------------------------------------*/

/// Implementation of a parametrized full order <b> unsteady Boussinesq problem </b> and preparation of the the reduced matrices for the online solve.
/** In this class are implemented the methods for the offline solve of a unsteady Boussinesq problem (two way coupling between momentum and energy equation)
and the for the generation of the reduced matrices for subsequent online solve, this class is a son
of the reduction problem class */
class UnsteadyBBTurb : public UnsteadyBB
{
public:
    // Constructors
    /// Null constructor
    UnsteadyBBTurb();

    /// Construct with argc and argv
    UnsteadyBBTurb(int argc, char* argv[]);
    ~UnsteadyBBTurb() { };

    /// List of snapshots for the solution for eddy viscosity
    PtrList<volScalarField> Nutfield;

    /// List of snapshots for the solution for eddy viscosity
    PtrList<volScalarField> liftfieldnut;

    /// Number of viscosity modes used for the projection
    label Nnutmodes;

    /// List of snapshots to form the homogeneous eddy viscosity fields
    PtrList<volScalarField> nutomFields;

    /// List of POD modes for eddy viscosity
    PtrList<volScalarField> nuTmodes;

    /// Eddy diffusivity field
    PtrList<volScalarField> alphatfield;

    /// Scalar to store the turbulent Prandtl number
    scalar Prt;

    /// Scalar to store the Prandtl number
    scalar Pr;

    /// Number of nut field modes
    label Nphi_nut;

    /// Create a Rbf splines for interpolation
    std::vector<SPLINTER::DataTable*> SAMPLES;

    /// Create a SAMPLES for interpolation
    std::vector<SPLINTER::RBFSpline*> rbfsplines;

    /** @name Reduced Matrices
     *
     */
    ///@{
    //
    /// Turbulent viscosity term
    Eigen::MatrixXd BT_matrix;

    /// Turbulent viscosity term
    List<Eigen::MatrixXd> CT2_matrix;

    /// Turbulent viscosity term
    List<Eigen::MatrixXd> CT1_matrix;

    /// Turbulent diffusivity term
    List<Eigen::MatrixXd> S_matrix;

    /// Total C Matrix
    List<Eigen::MatrixXd> C_total_matrix;

    /// Total B Matrix
    Eigen::MatrixXd B_total_matrix;
    ///@}

    /// Eddy viscosity field
    autoPtr<volScalarField> _nut;

    /// Eddy viscosity for the Spalart-Allmaras turbulence model
    autoPtr<volScalarField> _nuTilda;

    /// Turbulent Prandtl number
    autoPtr<dimensionedScalar> _Prt;

    /// Prandtl number
    autoPtr<dimensionedScalar> _Pr;

    /// Turbulent thermal diffusivity
    autoPtr<volScalarField> _alphat;

    // Functions

    /// @brief   Solves the CFD problem for parameters mu_now. Differently
    /// from UnsteadyBB, here the eddy viscosity field is also computed and saved.
    /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
    void truthSolve(List<scalar> mu_now);

    /// @brief   Solves the CFD problem for parameters mu_now. Differently
    /// from UnsteadyBB, here the eddy viscosity field is also computed and saved.
    void truthSolve(fileName folder = "./ITHACAOutput/Offline");

    //--------------------------------------------------------------------------
    /// Project using a supremizer approach
    ///
    /// @param[in]  folder     The folder used to save the reduced matrices.
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NPmodes    The number of pressure modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    ///
    void projectSUP(fileName folder, label NUmodes, label NPrghmodes, label NTmodes,
        label NSUPmodes, label Nnutmodes);

    /// @brief      BT added matrix for the turbulence treatement
    ///
    /// @param[in]  NUmodes    The N of velocity modes.s
    /// @param[in]  NSUPmodes  The N of supremizer modes.
    ///
    /// @return     BT matrix for turbulence treatment
    ///
    Eigen::MatrixXd BTturbulence(label NUmodes, label NSUPmodes);

    ///
    /// @brief      S added matrix for the thermal turbulence treatement
    ///
    /// @param[in]  NTmodes    The number of temperature modes.
    /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
    ///
    /// @return     S matrix for turbulence treatment
    ///
    List<Eigen::MatrixXd> temperatureTurbulenceTerm(label NTmodes, label Nnutmodes);
  
    /// @brief      CT1 added matrix for the turbulence treatement
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
    ///
    /// @return     CT1 matrix for turbulence treatment
    ///
    List < Eigen::MatrixXd > turbulenceTerm1(label NUmodes, label NSUPmodes,
            label Nnutmodes);

    ///
    /// @brief      CT2 added matrix for the turbulence treatement
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    /// @param[in]  Nnutmodes  The number of eddy viscosity modes.
    ///
    /// @return     CT2 matrix for turbulence treatment
    ///
    List<Eigen::MatrixXd> turbulenceTerm2(label NUmodes, label NSUPmodes, label Nnutmodes);
};

#endif
