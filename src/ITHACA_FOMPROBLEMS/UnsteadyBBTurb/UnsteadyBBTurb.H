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
#include "ithacaInterpolator.H"
#include "Modes.H"
#include <memory>
#include <iostream>
#include <functional>

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
private:
    /// bool for debug mode
    bool DEBUG_MODE = false;

    /// Dictionary to hold details about computing the nut
    dictionary viscDict;

    /// Supremizer approach (either "snapshot" or "modes")
    word supremizerApproach;

    //--------------------------------------------------------------------------
    /// @brief     Method to correctly prepare the modes (supremizer and lifting if necessary), be it PPE, Supremizer or VMB
    ///
    void prepareModes();

    /// @brief     Helper to load from disk or compute a matrix
    ///
    /// @param matrix
    /// @param folder
    /// @param nameBase
    /// @param suffix
    /// @param generator
    void loadOrCompute(Eigen::MatrixXd& matrix, const word& folder, const word& nameBase, const word& suffix, std::function<Eigen::MatrixXd()> generator);

    /// @brief     Helper to load from disk or compute a tensor
    ///
    /// @param tensor
    /// @param folder
    /// @param nameBase
    /// @param suffix
    /// @param generator
    void loadOrCompute(Eigen::Tensor<double, 3>& tensor, const word& folder, const word& nameBase, const word& suffix, std::function<Eigen::Tensor<double, 3>()> generator);

    void assembleCommonMatrices();

    void assembleSupremizerMatrices();

    void assemblePPEMatrices();

    void assemblePenaltyMatrices();

public:
    // Constructors
    /// Null constructor
    UnsteadyBBTurb();

    /// Construct with argc and argv
    UnsteadyBBTurb(int argc, char* argv[]);
    ~UnsteadyBBTurb() { };

    /// Wether the velocity coefficients derivative are used in the RBF interpolaltion for nut
    bool derivativeInRBF = true;

    // Lifted BC matrix. [NsampleRun x NBoundaryPatches]. Holds the values of the lifted BCs for each sample and each boundary patch
    Eigen::MatrixXd liftBCMatrix;

    // Shifted pressure modes - weird naming problems occur if we use Prghmodes, clashing with BBUnsteady, which has a PtrList<volScalarField>. Ideally
    // I would change to correct data type in UnsteadyBB (change to volScalarModes), but this requires more extensive changes.
    volScalarModes P_rghmodes;

    // Temperature modes
    volScalarModes L_Tmodes;

    // Eddy viscosity modes
    volScalarModes nutmodes;

    // Source term modes - Not implemented yet
    // volScalarModes Sourcemodes;

    // Source term snapshots - Not implemented yet
    // PtrList<volScalarField> Sourcefield;

    /// Number of modes used for the source term - Not implemented yet
    // label Nsourcemodes;

    /// List of snapshots for the solution for eddy viscosity
    PtrList<volScalarField> Nutfield;

    /// List of fluctuating eddy viscosity fields. One for each parameter sample
    PtrList<volScalarField> fluctNutfield;

    /// List of averaged eddy viscosity fields. One for each parameter sample
    PtrList<volScalarField> avgNutfield;

    // Dimension of the interpolation matrix A
    label dimA;

    /// Interpolators for nut coefficient interpolation
    std::vector<std::shared_ptr<ithacaInterpolator>> rbfSplines;

    /// Create a Eigen vector with the time instant for each snapshot
    List<Eigen::VectorXd> timeSnapshots;

    /// Common matrices to all methods (VMB/PPE/Supremizer)
    struct CommonMatrices {
        // Momentum equation
        Eigen::MatrixXd M, B, BTurb, K, H;
        Eigen::Tensor<double, 3> C, CTurb1, CTurb2, AveCTurb1, AveCTurb2;
        // Temperature equation
        Eigen::MatrixXd W, Y;
        Eigen::Tensor<double, 3> Q, YTurb, AveYTurb;
        // Aggregates
        Eigen::MatrixXd BTotal;
        Eigen::Tensor<double, 3> CTotal, CTotalAve;
    };
    struct PenaltyMatrices {
        // Penalty method terms
        Eigen::MatrixXd bcTempVec;
        List<Eigen::MatrixXd> bcVelVec, bcVelMat;
        Eigen::Tensor<double, 3> bcTempMat;
    };
    // Supremizer matrices
    struct SupremizerMatrices {
        // Supremizer matrices
        Eigen::MatrixXd P;
    };
    // PPE matrices
    struct PPEMatrices {
        // PPE matrices
        Eigen::MatrixXd D, HP, nuBC, timedepBC;
        Eigen::Tensor<double, 3> G, PPEAveCT1, PPECT1, PPEAveCT2, PPECT2, CTotalPPEAve, CTotalPPEFluct;
    };

    /// Shared pointer to the commin matrices. A shared pointer is used so its ownership can be shared, and thus the memory heavy UnsteadyBBTurb object can be destroyed. Still same treatment needed for the fields/modes
    std::shared_ptr<CommonMatrices> pCommonMatrices;
    /// Shared pointer to the penalty matrices
    std::shared_ptr<PenaltyMatrices> pPenaltyMatrices;
    /// Shared pointer to the supremizer matrices
    std::shared_ptr<SupremizerMatrices> pSupremizerMatrices;
    /// Shared pointer to the PPE matrices
    std::shared_ptr<PPEMatrices> pPPEMatrices;

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

    /// List of pointers for the online validation of the eddy viscosity solution
    PtrList<volScalarField> nutfield_on;

    /// List of pointers for the online validation of the shifted pressure solution
    PtrList<volScalarField> Prghfield_on;

    /// Number of velocity coefficient to skip in the RBF interpolation during the offline phase
    label skipRBFIndex = 0;

    /// boolean flag for time-dependent BCs
    bool timeDependentBC = false;

    /// @brief   Solves the CFD problem for parameters mu_now. Differently
    /// from UnsteadyBB, here the eddy viscosity field is also computed and saved.
    /// @param[in]  mu_now  The actual value of the parameter for this truthSolve. Used only
    /// @param[in]  nSample  The sample number for which the truth solve is performed.
    void truthSolve(const List<scalar> mu_now, label nSample);

    /// @brief   Solves the CFD problem for parameters mu_now. Differently
    /// from UnsteadyBB, here the eddy viscosity field is also computed and saved.
    void truthSolve(fileName folder = "./ITHACAOutput/Offline");

    //--------------------------------------------------------------------------
    /// @brief      Method to solve the supremizer either using the snapshot method or the pressure modes
    ///
    /// @param[in]  type       The type of supremizer approach ('snapshot' or 'modes')
    void solvesupremizer(word type = "snapshot");

    //--------------------------------------------------------------------------
    /// Gradient of pressure
    ///
    /// @param[in]  NU       The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.
    /// @param[in]  NSUP     The number of supremizer modes.
    ///
    /// @return     reduced matrix in Eigen::MatrixXd format for the Gradient of pressure term.
    ///
    Eigen::MatrixXd pressureGradientTerm(label NU, label NPrgh, label NSUP);

    //--------------------------------------------------------------------------
    /// Diffusive Term Energy Equation
    ///
    /// @param[in]  NUmodes    The number of velocity modes.
    /// @param[in]  NTmodes    The number of temperature modes.
    /// @param[in]  NSUPmodes  The number of supremizer modes.
    ///
    /// @return     reduced matrix in Eigen::MatrixXd format for the diffusion term.
    ///
    Eigen::MatrixXd diffusiveTermTemperature(label NUmodes, label NTmodes, label NSUPmodes);
    //--------------------------------------------------------------------------

    /// Divergence Term (supremizer approach)
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.
    /// @param[in]  NSUP  The number of supremizer modes.
    ///
    /// @return     reduced matrix in Eigen::MatrixXd format for the divergence term.
    ///
    Eigen::MatrixXd divergenceTerm(label NU, label NPrgh, label NSUP);
    //--------------------------------------------------------------------------
    /// Mass Term Energy Equation
    ///
    /// @param[in]  NT    The number of temperature modes.
    ///
    /// @return     reduced matrix in Eigen::MatrixXd format for the mass matrix.
    ///
    Eigen::MatrixXd massTermTemperature(label NT);

    //--------------------------------------------------------------------------
    /// Buoyant Term Momentum Equation
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NT    The number of temperature modes.
    /// @param[in]  NSUP  The number of supremizer modes.
    ///
    /// @return     reduced third order tensor in List <Eigen::MatrixXd> format for the convective term.
    ///
    Eigen::MatrixXd buoyancyTerm(label NU, label NT, label NSUP);

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
    /// @brief      Laplacian term of the PPE approach
    ///
    /// @param[in]  NPrgh    The number of pressure modes.
    Eigen::MatrixXd laplacianPPE(label NPrgh);

    //--------------------------------------------------------------------------
    /// @brief       Buoyancy term for the PPE equation
    ///
    /// @param[in]  NPrgh    The number of pressure modes.
    /// @param[in]  NT    The number of temperature modes.
    Eigen::MatrixXd buoyancyTermPPE(label NPrgh, label NT);

    //--------------------------------------------------------------------------
    /// @brief      Pressure BC1 matrix for the PPE approach
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.
    Eigen::MatrixXd BC1PPE(label NU, label NPrgh);

    //--------------------------------------------------------------------------
    /// @brief      Pressure BC2 matrix for the PPE approach
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.
    Eigen::Tensor<double, 3> BC2PPE(label NU, label NPrgh);

    //--------------------------------------------------------------------------
    /// @brief      Pressure BC3 matrix for the PPE approach
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.
    Eigen::MatrixXd BC3PPE(label NU, label NPrgh);

    //--------------------------------------------------------------------------
    /// @brief      PPE BC4 matrix, takes into account the time-dependent BC
    ///
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.
    Eigen::MatrixXd BC4PPE(label NU, label NPrgh);

    //--------------------------------------------------------------------------
    /// @brief      Divergence of the momentum equation for the PPE approach
    Eigen::Tensor<double, 3> divMomentum(label NU, label NPrgh);

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
    Eigen::Tensor<double, 3> convectiveTensorTemperature(label NU, label NT, label NSUP);

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
    /// @brief      PPE Average nut related tensor for the split turbulence treatment
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.

    Eigen::Tensor<double, 3> turbulencePPEAveTensor1(label NU, label NPrgh);

    //--------------------------------------------------------------------------
    /// @brief      PPE Fluctuating nut related tensor for the split turbulence treatment
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh   The number of pressure modes.
    /// @param[in]  Nnut  The number of eddy viscosity modes.

    Eigen::Tensor<double, 3> turbulencePPETensor1(label NU, label NPrgh, label Nnut);

    //--------------------------------------------------------------------------
    /// @brief      PPE nut related tensor for the split turbulence treatment
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.
    /// @param[in]  Nnut  The number of eddy viscosity modes.

    Eigen::Tensor<double, 3> turbulencePPETensor2(label NU, label NPrgh, label Nnut);

    //--------------------------------------------------------------------------
    /// @brief      PPE nut related tensor for the split turbulence treatment
    /// @param[in]  NU    The number of velocity modes.
    /// @param[in]  NPrgh    The number of pressure modes.

    Eigen::Tensor<double, 3> turbulencePPEAveTensor2(label NU, label NPrgh);

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
    /// @brief      Resets the fields after a truthSolve, reading from disk the initial conditions.
    ///
    void restart();

    //--------------------------------------------------------------------------
    /// @brief      Prepares the BC matrices for the penalty method for the temperature equation
    ///
    /// @param[in]  NT    The number of temperature modes.
    Eigen::MatrixXd bcTemperatureVec(const label NT);

    //--------------------------------------------------------------------------
    /// @brief      Prepares the BC matrices for the penalty method for the temperature equation
    ///
    /// @param[in]  NT    The number of temperature modes.
    Eigen::Tensor<double, 3> bcTemperatureMat(const label NT);

    //--------------------------------------------------------------------------
    /// @brief      Resizes the modes before the Online phase, according to the number of modes chosen
    ///
    void resizeModes();

    //--------------------------------------------------------------------------
    /// @brief      General projection method call. It then calls the specific method according to the choice of stabilization method.
    ///             Probably there is a cleaner way to implement this, but for now it is implemented with if statements.
    ///
    void project();
  };

#endif
