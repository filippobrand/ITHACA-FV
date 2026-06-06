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
/// Header for the general interface of ReducedODESystem class. Guarded under the namespace ITHACA_ODE.
/// This defines how any dynamical system should expose its equations to a generic ODE solver.
/// Ideally we should prepare both implicit and explicit virtual methods

#ifndef REDUCEDODESYSTEM_H
#define REDUCEDODESYSTEM_H

#include <Eigen/Dense>
#include "fvCFD.H"
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <unsupported/Eigen/AutoDiff>

class ReducedODESystem
{
public:
    // Constructor
    ReducedODESystem() = default;

    // Destructor
    virtual ~ReducedODESystem() = default;

    // Explicit virtual, evaluates the RHS
    virtual void evaluateRHS(const Eigen::VectorXd& state, Eigen::VectorXd& rhs, double t) const = 0;

    // Implicit virtual, evaluates the residual
    virtual void evaluateResidual(const Eigen::VectorXd& state, const Eigen::VectorXd& state_dot, Eigen::VectorXd& residual, double t) const = 0;
};

struct EigenImplicitEulerFunctor : public newton_argument<double> {
    const ReducedODESystem& system;
    const Eigen::VectorXd* y_old;
    double t_current;
    double dt;

    EigenImplicitEulerFunctor(const ReducedODESystem& sys, int N): newton_argument<double>(N, N), system(sys), y_old(nullptr) { }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
    {
        Eigen::VectorXd y_dot = (x - *y_old) / dt;
        system.evaluateResidual(x, y_dot, fvec, t_current);
        return 0;
    }

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
    {
        Eigen::NumericalDiff<EigenImplicitEulerFunctor> numDiff(*this);
        numDiff.df(x, fjac);
        return 0;
    }
};

struct EigenBDF2Functor : public newton_argument<double> {
    const ReducedODESystem& system;
    const Eigen::VectorXd* y_old;
    const Eigen::VectorXd* y_older;
    double t_current;
    double dt;

    EigenBDF2Functor(const ReducedODESystem& sys, int N): newton_argument<double>(N, N), system(sys), y_old(nullptr), y_older(nullptr) { }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const
    {
        Eigen::VectorXd y_dot = (3 * x - 4 * (*y_old) + (*y_older)) / (2 * dt);
        system.evaluateResidual(x, y_dot, fvec, t_current);
        return 0;
    }

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const
    {
        Eigen::NumericalDiff<EigenBDF2Functor> numDiff(*this);
        numDiff.df(x, fjac);
        return 0;
    }
};

class ImplicitEulerSolver
{
private:
    EigenImplicitEulerFunctor functor;
    Eigen::HybridNonLinearSolver<EigenImplicitEulerFunctor> solver;

public:
    ImplicitEulerSolver(const ReducedODESystem& system, int N): functor(system, N), solver(functor), y_old(Eigen::VectorXd::Zero(N)) { }

    Eigen::VectorXd y_old; // This is stored so that old state can be accessed to compute derivatives

    void solveStep(Eigen::VectorXd& y, double t, double dt, bool print_residual = false)
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

    void printResidual(Eigen::VectorXd& y, double t, double dt)
    {
        static const Color::Modifier red(Color::FG_RED);
        static const Color::Modifier green(Color::FG_GREEN);
        static const Color::Modifier def(Color::FG_DEFAULT);
        Eigen::VectorXd residual = Eigen::VectorXd::Zero(y.size());

        functor.operator()(y, residual);

        if (residual.norm() > 1e-5)
        {
            Info << red << "Time: " << t << " - Residual norm: " << residual.norm() << " - Reached in " << solver.iter << " iterations " << def << endl
                 << endl;
        } else
        {
            Info << green << "Time: " << t << " - Residual norm: " << residual.norm() << " - Reached in " << solver.iter << " iterations" <<def << endl
                 << endl;
        }
    }
};

class BDF2Solver
{
private:
    EigenBDF2Functor functor;
    Eigen::HybridNonLinearSolver<EigenBDF2Functor> solver;
    EigenImplicitEulerFunctor first_step_functor; // We will use this for the first step
    Eigen::HybridNonLinearSolver<EigenImplicitEulerFunctor> first_step_solver; // We will use this for the first step
    bool hasPreviousSteps;

public:
    BDF2Solver(const ReducedODESystem& system, int N): functor(system, N), solver(functor),
                                                       hasPreviousSteps(false), first_step_functor(system, N),
                                                       first_step_solver(first_step_functor),
                                                       y_old(Eigen::VectorXd::Zero(N)),
                                                       y_older(Eigen::VectorXd::Zero(N)) { }

    Eigen::VectorXd y_old; // This is stored so that old state can be accessed
    Eigen::VectorXd y_older; // Kept as a buffer
    Eigen::VectorXd y_current; // Kept as a buffer

    void solveStep(Eigen::VectorXd& y, double t, double dt, bool print_residual = false)
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

    void printResidual(Eigen::VectorXd& y, double t, double dt)
    {
        static const Color::Modifier red(Color::FG_RED);
        static const Color::Modifier green(Color::FG_GREEN);
        static const Color::Modifier def(Color::FG_DEFAULT);
        Eigen::VectorXd residual = Eigen::VectorXd::Zero(y.size());

        functor.operator()(y, residual);

        if (residual.norm() > 1e-5)
        {
            Info << red << "Time: " << t << " - Residual norm: " << residual.norm() << " - Reached in " << solver.iter << " iterations. Average residual: " << residual.mean() << def << "\n"
                 << "\n";
        } else
        {
            Info << green << "Time: " << t << " - Residual norm: " << residual.norm() << " - Reached in " << solver.iter << " iterations. Average residual: " << residual.mean() << def << "\n"
                 << "\n";
        }
    }
};

class TimeManager
{
private:
    double current_time;
    double dt;
    double initial_time;
    double final_time;
    double dt_save_coefficients;
    double dt_export_fields;
    int steps_between_saves;
    int steps_between_exports;
    int number_of_steps_to_solve;
    int number_of_saved_steps;
    int number_of_exported_steps;
    int number_of_steps_taken;
    bool adaptive_time_stepping;

public:
    // Constructors

    TimeManager() = default;

    TimeManager(double initial, double final, double timestep, double dt_save, double dt_export): current_time(initial), dt(timestep), initial_time(initial), final_time(final), dt_save_coefficients(dt_save), dt_export_fields(dt_export), adaptive_time_stepping(false)
    {
        steps_between_saves = static_cast<int>(dt_save_coefficients / dt);
        steps_between_exports = static_cast<int>(dt_export_fields / dt);
        number_of_steps_to_solve = static_cast<int>((final_time - initial_time) / dt);
        number_of_saved_steps = static_cast<int>((final_time - initial_time) / dt_save_coefficients);
        number_of_exported_steps = static_cast<int>((final_time - initial_time) / dt_export_fields);
        number_of_steps_taken = 0;
    }
    // Destructor
    ~TimeManager() = default;

    // Getters and utility functions
    double getCurrentTime() const { return current_time; }
    double getFinalTime() const { return final_time; }
    double getTimeStep() const { return dt; }
    bool isFinished() const { return current_time >= final_time - dt; }
    void advanceTime()
    {
        current_time += dt;
        number_of_steps_taken++;
    }
    int getNumberOfStepsTaken() const { return number_of_steps_taken; }
    bool shouldSaveCoefficients() const
    {
        return (number_of_steps_taken % steps_between_saves == 0);
    };
    // bool checkExportFields(double time) const;
};

#endif
