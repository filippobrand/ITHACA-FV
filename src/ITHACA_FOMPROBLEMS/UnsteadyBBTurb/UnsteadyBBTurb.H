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
#include "Modes.H"
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

    // Lifted BC matrix. [NsampleRun x NBoundaryPatches]. Holds the values of the lifted BCs for each sample and each boundary patch
    Eigen::MatrixXd liftBCMatrix;

    // Shifted pressure modes
    volScalarModes P_rghmodes;

    // Temperature modes
    volScalarModes L_Tmodes;

    // Eddy viscosity modes
    volScalarModes nutmodes;

    /// List of snapshots for the solution for eddy viscosity
    PtrList<volScalarField> Nutfield;

    /// List of fluctuating eddy viscosity fields. One for each parameter sample
    PtrList<volScalarField> fluctNutfield;

    /// List of averaged eddy viscosity fields. One for each parameter sample
    PtrList<volScalarField> avgNutfield;

    /// Number of viscosity modes used for the projection
    label Nnutmodes;

    // Dimension of the interpolation matrix A
    label dimA;

    /// sample object for RBF interpolation
    std::vector<SPLINTER::DataTable*> samples;

    /// RBF spline objects for interpolation
    std::vector<SPLINTER::RBFSpline*> rbfSplines;

    /// RBF shape parameters vector
    Eigen::VectorXd radii;

    /// RBF Function radius
    double radius = 0.1;

    /// Vector with column means for normalization during RBF for matrix A
    Eigen::VectorXd meanA;

    /// Vector with column std for normalization during RBF for matrix A
    Eigen::VectorXd stdA;

    /// Vector with column means for normalization during RBF for matrix G
    Eigen::VectorXd meanG;

    /// Vector with column std for normalization during RBF for matrix G
    Eigen::VectorXd stdG;

    /// Create a Eigen vector with the time instant for each snapshot
    List<Eigen::VectorXd> timeSnapshots;

    /** @name Reduced Matrices
     *
     */
    ///@{
    //
    /// Turbulent viscosity term
    Eigen::MatrixXd BT_matrix;

    /// Turbulent viscosity tensor CT1 for the momentum equation
    Eigen::Tensor<double, 3> CT1_tensor;

    /// Turbulent viscosity tensor CT2 for the momentum equation
    Eigen::Tensor<double, 3> CT2_tensor;

    /// Turbulent viscosity tensor CT1_ave for the momentum equation
    Eigen::Tensor<double, 3> CT1_ave_tensor;

    /// Turbulent viscosity tensor CT2_ave for the momentum equation
    Eigen::Tensor<double, 3> CT2_ave_tensor;

    /// Turbulent diffusivity tensor for the temperature equation YT
    Eigen::Tensor<double, 3> YT_tensor;

    /// Turbulent diffusivity tensor for the temperature equation YT average
    Eigen::Tensor<double, 3> YT_ave_tensor;

    /// Total CT tensor
    Eigen::Tensor<double, 3> C_total_tensor;

    /// Total CT average tensor
    Eigen::Tensor<double, 3> C_total_ave_tensor;

    /// Convective tensor Q
    Eigen::Tensor<double, 3> Q_tensor;

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

    /// bool for debug mode
    bool DEBUG_MODE = false;

    /// @brief   Solves the CFD problem for parameters mu_now. Differently
    /// from UnsteadyBB, here the eddy viscosity field is also computed and saved.
    /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
    /// @param[in]  nSample  The sample number for which the truth solve is performed.
    void truthSolve(List<scalar> mu_now, label nSample);

    /// @brief   Solves the CFD problem for parameters mu_now. Differently
    /// from UnsteadyBB, here the eddy viscosity field is also computed and saved.
    void truthSolve(fileName folder = "./ITHACAOutput/Offline");

    //--------------------------------------------------------------------------
    /// @brief      Method to solve the supremizer either using the snapshot method or the pressure modes
    ///
    /// @param[in]  type       The type of supremizer approach ('snapshot' or 'modes')
    void solvesupremizer(word type = "snapshot");

    //--------------------------------------------------------------------------
    /// Project using a supremizer approach
    ///
    /// @param[in]  folder     The folder used to save the reduced matrices.
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.
    /// @param[in]  NSUP  The number of supremizer modes.
    ///
    void projectSUP(fileName folder, label NU, label NPrgh, label NT, label NSUP, label Nnut);

    //--------------------------------------------------------------------------
    /// Gradient of pressure
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NPmodes    The number of pressure modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    ///
    /// @return     reduced matrix in Eigen::MatrixXd format for the Gradient of pressure term.
    ///
    Eigen::MatrixXd pressure_gradient_term(label NUmodes, label NPrghmodes, label NSUPmodes);

    //--------------------------------------------------------------------------
    /// Diffusive Term Energy Equation
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NTmodes    The number of temperature modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    ///
    /// @return     reduced matrix in Eigen::MatrixXd format for the diffusion term.
    ///
    Eigen::MatrixXd diffusive_term_temperature(label NUmodes, label NTmodes, label NSUPmodes);
    //--------------------------------------------------------------------------

    /// Divergence Term (supremizer approach)
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NPmodes    The number of pressure modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    ///
    /// @return     reduced matrix in Eigen::MatrixXd format for the divergence term.
    ///
    Eigen::MatrixXd divergence_term(label NUmodes, label NPrghmodes, label NSUPmodes);

    //--------------------------------------------------------------------------
    /// Mass Term Energy Equation
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NTmodes    The number of temperature modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    ///
    /// @return     reduced matrix in Eigen::MatrixXd format for the mass matrix.
    ///
    Eigen::MatrixXd mass_term_temperature(label NUmodes, label NTmodes, label NSUPmodes);

    //--------------------------------------------------------------------------
    /// Buoyant Term Momentum Equation
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NTmodes    The number of temperature modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    ///
    /// @return     reduced third order tensor in List <Eigen::MatrixXd> format for the convective term.
    ///
    Eigen::MatrixXd buoyant_term(label NUmodes, label NTmodes, label NSUPmodes);

    //--------------------------------------------------------------------------
    /// @brief      BT added matrix for the turbulence treatement
    ///
    /// @param[in]  NU    The N of velocity modes.s
    /// @param[in]  NSUP  The N of supremizer modes.
    ///
    /// @return     BT matrix for turbulence treatment
    ///
    Eigen::MatrixXd BTturbulence(label NU, label NSUP);

    //--------------------------------------------------------------------------
    /// @brief      S added tensor for the thermal turbulence treatement
    ///
    /// @param[in]  NT    The number of temperature modes.
    /// @param[in]  Nnut  The number of eddy viscosity modes.
    ///
    /// @return     S matrix for turbulence treatment
    ///
    Eigen::Tensor<double, 3> temperatureTurbulenceTensor(label NT, label Nnut);

    //--------------------------------------------------------------------------
    /// @brief      ct1 tensor for the turbulence treatement
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NSUP  The number of supremizer modes.
    /// @param[in]  Nnut  The number of eddy viscosity modes.
    ///
    /// @return     ct1 tensor for the turbulence treatement
    ///
    Eigen::Tensor<double, 3> turbulenceTensor1(label NU, label NSUP, label Nnut);

    //--------------------------------------------------------------------------
    /// @brief      Method to compute one of the turbulence eddy viscosity tensors
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NSUP  The number of supremizer modes.
    /// @param[in]  Nnut  The number of eddy viscosity modes.
    ///
    /// @return     ct2 tensor for the turbulence treatement
    ///
    Eigen::Tensor<double, 3> turbulenceTensor2(label NU, label NSUP, label Nnut);

    //--------------------------------------------------------------------------
    /// @brief      Method to compute the convective tensor for the temperature equation
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NT    The number of temperature modes.
    /// @param[in]  NSUP  The number of supremizer modes.
    ///
    /// @return     ct2 tensor for the turbulence treatement
    ///
    Eigen::Tensor<double, 3> convective_tensor_temperature(label NU, label NT, label NSUP);

    //--------------------------------------------------------------------------
    /// @brief      ct1Ave added tensor for approximation of the averaged part of the eddy viscosity
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NSUP  The number of supremizer modes.
    ///
    /// @return     ct1Ave tensor corresponding to the approximation of the averaged part
    ///  of the eddy viscosity for the turbulence treatement
    ///
    Eigen::Tensor<double, 3> turbulenceAveTensor1(label NU, label NSUP);

    //--------------------------------------------------------------------------
    /// @brief      ct2Ave added tensor for approximation of the averaged part of the eddy viscosity
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NSUP  The number of supremizer modes.
    ///
    /// @return     ct2Ave tensor corresponding to the approximation of the averaged part
    ///  of the eddy viscosity for the turbulence treatement
    ///
    Eigen::Tensor<double, 3> turbulenceAveTensor2(label NU, label NSUP);

    //--------------------------------------------------------------------------
    /// @brief      Y_{T2} added tensor for approximation of the averaged part of the eddy viscosity in the temperature equation
    ///
    /// @param[in]  NT    The number of temperature modes.
    ///
    /// @return     Y_{T2} tensor corresponding to the approximation of the averaged part
    ///  of the eddy viscosity for the turbulence treatement in the energy equation
    ///
    Eigen::Tensor<double, 3> turbulenceTemperatureAveTensor(label NT);

    //--------------------------------------------------------------------------
    /// @brief      A method to compute the two matrices needed for the RBF interpolation by combining
    /// the velocity L2 projection coefficients and their time derivatives [a, a_dot]. The time derivatives
    /// are computed using a finite difference scheme according to the list of time instants provided snapshotTimes.
    ///
    /// @param[in]  A            The velocity L2 projection coefficients
    /// @param[in]  G            The eddy viscoisty L2 projection coefficients
    /// @param[in]  snapshotTimes A list of Eigen vectors containing the time instants for each snapshot
    ///
    /// @return     The matrix of the combined velocity L2 projection coefficients starting from the second
    /// snapshots for each parameter sample together with their time derivatives and also the matrix of the
    /// eddy viscosity L2 projection coefficients associated with the first matrix
    ///
    List<Eigen::MatrixXd> velDerivativeCoeff(const Eigen::MatrixXd& A, const Eigen::MatrixXd& G, const List<Eigen::VectorXd>& snapshotTimes);

    //--------------------------------------------------------------------------
    /// @brief     Splits the eddy viscosity field snapshots in average and fluctuating part
    ///
    void splitEddyViscositySnapshots();

    //--------------------------------------------------------------------------
    /// @brief     Carries out the RBF interpolation
    ///
    void offlineRBFInterpolation();

    //--------------------------------------------------------------------------
    /// @brief     Debugging tool to print the dimension of a given Eigen::MatrixXd
    ///
    /// @param[in]  matrix    The matrix to print the size of
    /// @param[in]  nameMatrix  The name of the matrix to print
    ///
    void printMatrixSize(const Eigen::MatrixXd& matrix, const std::string& nameMatrix)
    {
        Info << "--- Matrix " << nameMatrix << " size: " << matrix.rows() << " x " << matrix.cols() << endl;
    }
};

#endif
