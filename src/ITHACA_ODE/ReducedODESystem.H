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
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

/// Minimal Eigen-compatible nonlinear functor base.
/// Keeping this local avoids dependencies on internal ITHACA/OpenFOAM headers.
template<typename ScalarT, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
class EigenNonlinearFunctor
{
    public:
        typedef ScalarT Scalar;
        enum
        {
            InputsAtCompileTime = NX,
            ValuesAtCompileTime = NY
        };
        typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
        typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
        typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime>
        JacobianType;

        EigenNonlinearFunctor() : m_inputs(InputsAtCompileTime),
            m_values(ValuesAtCompileTime) { }

        EigenNonlinearFunctor(int inputs, int values) : m_inputs(inputs),
            m_values(values) { }

        int inputs() const
        {
            return m_inputs;
        }

        int values() const
        {
            return m_values;
        }

    private:
        int m_inputs;
        int m_values;
};

class ODESolver
{
    public:
        virtual ~ODESolver() = default;

        virtual void solveStep(Eigen::VectorXd& y, double t, double dt,
                               bool print_residual = false) = 0;

        virtual void printResidual(Eigen::VectorXd& y, double t, double dt) = 0;

        virtual Eigen::VectorXd getLatestDerivative(const Eigen::VectorXd& x) const = 0;

        virtual void reset() = 0;
};

class ReducedODESystem
{
    public:
        // Constructor
        ReducedODESystem() = default;

        // Destructor
        virtual ~ReducedODESystem() = default;

        // Explicit virtual, evaluates the RHS
        virtual void evaluateRHS(const Eigen::VectorXd& state, Eigen::VectorXd& rhs,
                                 double t) const = 0;

        // Implicit virtual, evaluates the residual
        virtual void evaluateResidual(const Eigen::VectorXd& state,
                                      const Eigen::VectorXd& state_dot, Eigen::VectorXd& residual,
                                      double t) const = 0;
};

struct EigenImplicitEulerFunctor : public EigenNonlinearFunctor<double>
{
    const ReducedODESystem& system;
    const Eigen::VectorXd* y_old;
    double t_current;
    double dt;

    EigenImplicitEulerFunctor(const ReducedODESystem& sys,
                              int N): EigenNonlinearFunctor<double>(N, N), system(sys), y_old(nullptr) { }

    Eigen::VectorXd computeTimeDerivative(const Eigen::VectorXd& x) const;

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;
};

struct EigenBDF2Functor : public EigenNonlinearFunctor<double>
{
    const ReducedODESystem& system;
    const Eigen::VectorXd* y_old;
    const Eigen::VectorXd* y_older;
    double t_current;
    double dt;

    EigenBDF2Functor(const ReducedODESystem& sys,
                     int N): EigenNonlinearFunctor<double>(N, N), system(sys), y_old(nullptr),
        y_older(nullptr) { }

    Eigen::VectorXd computeTimeDerivative(const Eigen::VectorXd& x) const;

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const;

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const;
};

class ImplicitEulerSolver : public ODESolver
{
    private:
        EigenImplicitEulerFunctor functor;
        Eigen::HybridNonLinearSolver<EigenImplicitEulerFunctor> solver;

    public:
        ImplicitEulerSolver(const ReducedODESystem& system, int N): functor(system, N),
            solver(functor), y_old(Eigen::VectorXd::Zero(N)) { }

        Eigen::VectorXd
        y_old; // This is stored so that old state can be accessed to compute derivatives

        void solveStep(Eigen::VectorXd& y, double t, double dt,
                       bool print_residual = false) override;

        void printResidual(Eigen::VectorXd& y, double t, double dt) override;

        Eigen::VectorXd getLatestDerivative(const Eigen::VectorXd& x) const;

        void reset() override
        {
            y_old.setZero();
        }
};

class BDF2Solver : public ODESolver
{
    private:
        EigenBDF2Functor functor;
        Eigen::HybridNonLinearSolver<EigenBDF2Functor> solver;
        EigenImplicitEulerFunctor
        first_step_functor; // We will use this for the first step
        Eigen::HybridNonLinearSolver<EigenImplicitEulerFunctor>
        first_step_solver; // We will use this for the first step
        bool hasPreviousSteps;

    public:
        BDF2Solver(const ReducedODESystem& system, int N): functor(system, N),
            solver(functor),
            hasPreviousSteps(false), first_step_functor(system, N),
            first_step_solver(first_step_functor),
            y_old(Eigen::VectorXd::Zero(N)),
            y_older(Eigen::VectorXd::Zero(N)) { }

        Eigen::VectorXd y_old; // This is stored so that old state can be accessed
        Eigen::VectorXd y_older; // Kept as a buffer
        Eigen::VectorXd y_current; // Kept as a buffer

        void solveStep(Eigen::VectorXd& y, double t, double dt,
                       bool print_residual = false) override;

        void printResidual(Eigen::VectorXd& y, double t, double dt) override;

        Eigen::VectorXd getLatestDerivative(const Eigen::VectorXd& x) const;

        void reset() override
        {
            hasPreviousSteps = false;
            y_old.setZero();
            y_older.setZero();
            y_current.setZero();
        }
};

#endif
