#include "TimeLoop.H"
#include "ReducedUnsteadyBBTurb.H"
#include <iostream>

template <typename ROMType>
void TimeLoop<ROMType>::run() // Runs an entire ODE time loop,
{
    int saved_steps = 0;

    while (!time_manager_.isFinished())
    {
        time_manager_.advanceTime();
        solver_.solveStep(current_state_, time_manager_.getCurrentTime(),
                          time_manager_.getActualTimeStep(), false);

        if (time_manager_.shouldSaveCoefficients())
        {
            updateHistory(saved_steps);
            saved_steps++;
            time_manager_.updateAfterSave();
        }
    }
}

//-Update the history
template <typename ROMType>
void TimeLoop<ROMType>::updateHistory(int& saved_steps)
{
    state_history_(saved_steps, 0) = time_manager_.getCurrentTime();
    state_history_.block(saved_steps, 1, 1,
                         num_of_variables_) = current_state_.transpose();

    if (rom_.hasMonitors())
    {
        monitor_history_(saved_steps, 0) = time_manager_.getCurrentTime();
        monitor_history_.block(saved_steps, 1, 1,
                               rom_.getMonitorSize()) = rom_.computeMonitorQuantities(current_state_,
                                   time_manager_.getCurrentTime()).transpose();
    }

    solver_.printResidual(current_state_, time_manager_.getCurrentTime(),
                          time_manager_.getActualTimeStep());
}


template class
TimeLoop<ReducedUnsteadyBBTurb>; // Explicit instantiation for Reactor class
