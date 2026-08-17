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

\*---------------------------------------------------------------------------*/

/// \file
/// Source file of the ReducedODESystem class, containing the implementation of the TimeManager class and the ODEStructurePPETurb adapter class.

#include "ReducedODESystem.H"
#include "colormod.H"
#include "IOstreams.H"

/// ----- Implicit Euler Functor -----

int EigenImplicitEulerFunctor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd y_dot = (x - (*y_old)) / dt;
    system.evaluateResidual(x, y_dot, fvec, t_current);
    return 0;
}

int EigenImplicitEulerFunctor::df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<EigenImplicitEulerFunctor> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

/// ----- BDF2 Functors -----

int EigenBDF2Functor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
{
    Eigen::VectorXd y_dot = (3 * x - 4 * (*y_old) + (*y_older)) / (2 * dt);
    system.evaluateResidual(x, y_dot, fvec, t_current);
    return 0;
}

int EigenBDF2Functor::df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
{
    Eigen::NumericalDiff<EigenBDF2Functor> numDiff(*this);
    numDiff.df(x, fjac);
    return 0;
}

// ----- Solvers -----

void ImplicitEulerSolver::solveStep(Eigen::VectorXd& y, double t, double dt, bool print_residual)
{
    y_old = y; // Store the old state before updating it
    functor.y_old = &y_old;
    functor.t_current = t;
    functor.dt = dt;
    Eigen::VectorXd y_guess = y;
    solver.solve(y_guess);
    y = y_guess;

    if (print_residual)
    {
        // Assuming your functor can calculate its own residual vector
        printResidual(y, t, dt);
    }
}

void ImplicitEulerSolver::printResidual(Eigen::VectorXd& y, double t, double dt)
{
    // static const Color::Modifier red(Color::FG_RED);
    // static const Color::Modifier green(Color::FG_GREEN);
    // static const Color::Modifier def(Color::FG_DEFAULT);
    Eigen::VectorXd residual = Eigen::VectorXd::Zero(y.size());

    functor.operator()(y, residual);

    if (residual.norm() > 1e-5)
    {
        Foam::Info << "Time: " << t << " - Residual norm: " << residual.norm() << " - Reached in " << solver.iter << " iterations " << "\n"
             << "\n";
    } else
    {
        Foam::Info << "Time: " << t << " - Residual norm: " << residual.norm() << " - Reached in " << solver.iter << " iterations" << "\n" << "\n";
    }
}

void BDF2Solver::solveStep(Eigen::VectorXd& y, double t, double dt, bool print_residual)
{
    if (!hasPreviousSteps)
    {
        y_old = y;

        first_step_functor.y_old = &y_old;
        first_step_functor.t_current = t;
        first_step_functor.dt = dt;

        y_current = y;
        first_step_solver.solve(y_current);

        y_older = y_old;
        y_old = y_current;
        y = y_current;

        hasPreviousSteps = true;
    } else
    {
        functor.y_old = &y_old;
        functor.y_older = &y_older;
        functor.t_current = t;
        functor.dt = dt;

        y_current = y;
        solver.solve(y_current);

        if (print_residual)
        {
            printResidual(y_current, t, dt);
        }
        y_older = y_old;
        y_old = y_current;
        y = y_current;
    }
}

void BDF2Solver::printResidual(Eigen::VectorXd& y, double t, double dt)
{
    // static const Color::Modifier red(Color::FG_RED);
    // static const Color::Modifier green(Color::FG_GREEN);
    // static const Color::Modifier def(Color::FG_DEFAULT);
    Eigen::VectorXd residual = Eigen::VectorXd::Zero(y.size());

    functor.operator()(y, residual);

    if (residual.norm() > 1e-5)
    {
        Foam::Info << "Time: " << t << " - Residual norm: " << residual.norm() << " - Reached in " << solver.iter << " iterations. Average residual: " << residual.mean() << "\n"
             << "\n";
    } else
    {
        Foam::Info << "Time: " << t << " - Residual norm: " << residual.norm() << " - Reached in " << solver.iter << " iterations. Average residual: " << residual.mean() << "\n"
             << "\n";
    }
}
