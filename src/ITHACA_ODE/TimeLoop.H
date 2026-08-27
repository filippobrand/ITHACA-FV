// This header and its implementation file define the classes and functions to manage the ODE Time loop.

#ifndef TIMELOOP_H
#define TIMELOOP_H

#include <iostream>
#include <Eigen/Dense>
#include "ReducedODESystem.H"
#include "TimeManager.H"

template <typename ROMType>
class TimeLoop
{
        // This class takes in a TimeManager instance and a Reactor instance.
        // The TimeManager contains details about time steps, initial and final times and so on
        // The Reactor instance contains the physics of the model, and the methods such as evaluateResidual().
    private:
        ODESolver& solver_;
        ROMType& rom_;
        TimeManager time_manager_;

        Eigen::MatrixXd state_history_;
        Eigen::MatrixXd monitor_history_;
        Eigen::VectorXd current_state_;
        Eigen::VectorXd initial_state_;

        int num_of_variables_;

    public:
        TimeLoop(ODESolver &solver, ROMType &rom, TimeManager time_manager)
            : solver_(solver), rom_(rom), time_manager_(time_manager)
        {
            num_of_variables_ = rom_.getStateSize();
            state_history_.resize(time_manager_.getNumberOfStepsToSave(),
                                  1 + num_of_variables_);
            state_history_.setZero();

            if (rom_.hasMonitors())
            {
                monitor_history_.resize(time_manager_.getNumberOfStepsToExport(),
                                        1 + rom_.getMonitorSize());
                monitor_history_.setZero();
            }
            else
            {
                monitor_history_.resize(0, 0); // No monitors, so empty matrix
            }

            std::cout << "TimeLoop initialized with " << num_of_variables_ << " variables, "
                      << time_manager_.getNumberOfStepsToSave() << " steps to save, and "
                      << time_manager_.getNumberOfStepsToExport() << " steps to export." << std::endl;
        }

        Eigen::MatrixXd getStateHistory() const
        {
            return state_history_;
        }

        Eigen::MatrixXd getMonitorHistory() const
        {
            return monitor_history_;
        }

        void assignInitialState(const Eigen::VectorXd &initial_state)
        {
            if (initial_state.size() != num_of_variables_)
            {
                throw std::invalid_argument("Initial state size does not match the number of variables in the rom.");
            }

            initial_state_ = initial_state;
            current_state_ = initial_state;
        }

        void run();

        void updateHistory(int& saved_steps);

        void printLatestState() const
        {
            std::cout << "Latest state at time " << time_manager_.getCurrentTime() << ": "
                      << current_state_.transpose() << std::endl;
        }
};

#endif // TIMELOOP_H
