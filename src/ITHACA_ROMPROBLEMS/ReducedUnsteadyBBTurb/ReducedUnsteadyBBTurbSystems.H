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
    ReducedUnsteadyBBTurbSystems
Description
    Contains different ODE systems for the reduced unsteady, turbulent, buoyant Navier-stokes problem
SourceFiles
    ReducedUnsteadyBBTurbSystems.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedUnsteadyBBTurbSystems class.
/// \dir
/// Directory containing the header and source files for the ReducedUnsteadyBBTurb class

#ifndef ReducedUnsteadyBBTurbSystems_H
#define ReducedUnsteadyBBTurbSystems_H

#include <Eigen/Dense>
#include "ReducedODESystem.H"
#include "ReducedUnsteadyBBTurb.H"
#include "BoundaryConditions.H"
#include "UnsteadyBBTurb.H"


class ReducedUnsteadyBBTurbSystem : public ReducedODESystem
{
    protected:
        const ReducedUnsteadyBBTurb& rom;

    public:
        explicit ReducedUnsteadyBBTurbSystem(const ReducedUnsteadyBBTurb& rom): rom(
                rom) { };
        virtual ~ReducedUnsteadyBBTurbSystem() = default;

        virtual void evaluateRHS(
            const Eigen::VectorXd& state,
            Eigen::VectorXd& rhs,
            double t) const = 0;

        virtual void evaluateResidual(
            const Eigen::VectorXd& state,
            const Eigen::VectorXd& state_dot,
            Eigen::VectorXd& residual,
            double t) const = 0;
};

class SystemPPEGunzburger : public ReducedUnsteadyBBTurbSystem
{
    private:
        const UnsteadyBBTurb::CommonMatrices& commonMat_;
        const UnsteadyBBTurb::PPEMatrices& ppeMat_;
        const UnsteadyBBTurb::GunzburgerMatrices& gunzMat_;
        const BoundaryConditions& bc_;

    public:
        explicit SystemPPEGunzburger(const ReducedUnsteadyBBTurb& rom)
            :
            ReducedUnsteadyBBTurbSystem(rom),
            commonMat_(*rom.pCommonMatrices),
            ppeMat_(*rom.pPPEMatrices),
            gunzMat_(*rom.pGunzburgerMatrices),
            bc_(rom.getBoundaryConditions()) { };

        virtual ~SystemPPEGunzburger() = default;


        /// @brief Evaluates the right-hand side of the ODE system
        /// @param state Current state vector
        /// @param rhs  Vector to store the computed right-hand side
        /// @param t Current time
        void evaluateRHS(
            const Eigen::VectorXd& state,
            Eigen::VectorXd& rhs,
            double t) const override;

        /// @brief Evaluates the residual of the ODE system for implicit time integration
        /// @param state Current state vector
        /// @param state_dot Time derivative of the state vector
        /// @param residual Vector to store the computed residual
        /// @param t Current time
        void evaluateResidual(
            const Eigen::VectorXd& state,
            const Eigen::VectorXd& state_dot,
            Eigen::VectorXd& residual,
            double t) const override;
};

#endif // ReducedUnsteadyBBTurbSystems_H
