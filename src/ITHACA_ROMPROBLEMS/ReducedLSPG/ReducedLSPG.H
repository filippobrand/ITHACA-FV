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
    ReducedLSPG
Description
    Class where the base of a Least-Squares Petrov-Galerkin reduced problem for the unsteady flow is implemented.
SourceFiles
    ReducedLSPG.C
\*---------------------------------------------------------------------------*/

/// \file
/// Header file of the ReducedLSPG class.
/// \dir
/// Directory containing the header and source files for the ReducedLSPG class

#ifndef ReducedLSPG_H
#define ReducedLSPG_H

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "ITHACAcontext.H"
#include <Eigen/Dense>
#include <memory>

/*
Class implementation for a residual object, used to fill, print, manipulate the residual in the FOM
*/
class ResidualLSPG
{
   public:
     //- No default constructor
      ResidualLSPG() = delete;

      //- Constructor that takes the size of the residual vector
      explicit ResidualLSPG(const int size)
      : current_residual_(Eigen::VectorXd::Zero(size))
      {}

      virtual ~ResidualLSPG() = default;

      //- Set the residual to zero
      void setZero()
      {
          current_residual_.setZero();
      }

      //- Norm the residual
      double norm() const
      {
          return current_residual_.norm();
      }

   private:
      //- Current residual vector, stored as an Eigen vector
      Eigen::VectorXd current_residual_;
};

/*
Base class for the Reduced Least-Squares Petrov-Galerkin (LSPG) problem
Here we only implement the constructor and the destructor, which initialize quantities that are common to all (most) reduced problems.
Detailed implementation, such as field creation, is left to the derived classes.
*/
class ReducedLSPG
{
  public:
    //- Empty constructor - Cannot be used to create an object of this class
    ReducedLSPG() = delete;
    
    //- Constructor that takes, other than the FOM, argc and argv
    explicit ReducedLSPG(std::shared_ptr<ITHACAcontext> context);

    //- Default destructor
    virtual ~ReducedLSPG() = default;

  protected:
    //- Shared pointer to the ITHACAcontext object, which contains the mesh, time, fvOptions, and other necessary objects for the FOM problem.
    std::shared_ptr<ITHACAcontext> context_;

    

    argList& args() const
    {
        return context_->args();
    }
    
    fvMesh& mesh() const
    {
        return context_->mesh();
    }

    Time& runTime() const
    {
        return context_->runTime();
    }

    pimpleControl& pimple() const
    {
        return context_->pimple();
    }

    fv::options& fvOptions() const
    {
        return context_->fvOptions();
    }

    IOMRFZoneList& MRF() const
    {
        return context_->MRF();
    }

    // *** Time step management - These are called by openfoam Time object *** //
    //- Whether to adjust the time step based on the Courant number (or custom criteria?)
    bool adjustTimeStep;

    //- Maximum Courant number allowed
    scalar maxCo;

    //- Maximum time step allowed
    scalar maxDeltaT;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
