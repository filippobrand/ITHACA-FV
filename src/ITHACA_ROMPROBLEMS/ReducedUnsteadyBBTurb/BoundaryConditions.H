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
    BoundaryConditions
Description
    A class to handle various boundary conditions for Reduced Order Models (time-dependent, temperature, ...).
    At the moment this is only implemented for linear time-dependent BCs, where both temperature and velocity are given.
    The input BCs are given as Eigen::MatrixXd, where the first column is a series of time steps, and the other columns are the BC values at those time steps.
SourceFiles
    BoundaryConditions.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the BoundaryConditions class.
/// \dir
/// Directory containing the header and source files for the BoundaryConditions class

#ifndef BoundaryConditions_H
#define BoundaryConditions_H

#include "fvCFD.H"
#include "IOmanip.H"
#include <Eigen/Dense>

/*---------------------------------------------------------------------------*\
                        Class BoundaryConditions Declaration
\*---------------------------------------------------------------------------*/

/// Class where a general treatment of boundary conditions for ROMs is implemented
/** The idea is to have the possibility to group temperature BC, velocity BC and have methods to interpolate them
 * in case of time-dependent BCs. Inside the ROM solution, an update of the BC is called via the method updateTimeDependentBC.
 *  In the future, maybe different time-dep expression can be implemented. */
class BoundaryConditions
{
private:
    /// Matrix that stores the boundary conditions for temperature
    Eigen::MatrixXd temperatureBCMatrix_;

    /// Matrix that stores the boundary conditions for velocity
    Eigen::MatrixXd velocityBCMatrix_;

    /// Vector that stores the time instances at which the temperature BCs change
    Eigen::VectorXd timestepsTempBC_;

    /// Vector that stores the time instances at which the velocity BCs change
    Eigen::VectorXd timestepsVelBC_;
    
    /// Word that stores how to treat time-dependent BCs. At the moment only "linear" is implemented
    word timeDepMethod_;

    /// Boolean to know if there are temperature BCs
    bool hasTemperatureBC_ = false;

    /// Boolean to know if there are velocity BCs
    bool hasVelocityBC_ = false;

    /// Boolean to know if temperature BCs are time-dependent
    bool temperatureTimeDep_ = false;

    /// Boolean to know if velocity BCs are time-dependent
    bool velocityTimeDep_ = false;

public:
    // Constructors
    BoundaryConditions(const Eigen::MatrixXd& velocityBC, const Eigen::MatrixXd& temperatureBC, word timeDepMethod);
    BoundaryConditions(const Eigen::VectorXd& velocityBC, const Eigen::VectorXd& temperatureBC);

    // Attributes
    /// Vector that stores the current temperature BCs
    Eigen::VectorXd currentTemperatureBC;

    /// Vector that stores the current velocity BCs
    Eigen::VectorXd currentVelocityBC;

    // Methods

    /// @brief    Method to update time-dependent boundary conditions
    ///
    /// @param[in]  currentTime  Scalar of the current time
    ///
    void updateTimeDependentBC(const scalar currentTime);

    /// @brief    Method to linearly interpolate a value between two values
    ///
    double linearInterpolate(const double t0, const double t1, const double v0, const double v1, const scalar currentTime);

  };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
