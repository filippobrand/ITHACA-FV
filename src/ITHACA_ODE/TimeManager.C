#include <iostream>
#include <fstream>
#include "TimeManager.H"
// Here I try minimal dependence on ITHACA/OpenFOAM, as I steal this code
// elsewhere

const double TOLERANCE = 1e-10; // Tolerance for floating-point comparisons

TimeManager::TimeManager(double initial, double final, double timestep,
                         double dt_save, double dt_export): initial_time(initial), final_time(final),
    dt(timestep), dt_save_coefficients(dt_save), dt_export_fields(dt_export),
    current_time(initial), adaptive_time_stepping(false)
{
    number_of_steps_taken = 0;
    next_time_to_save = initial_time + dt_save_coefficients;
    actual_dt = dt;
    steps_between_saves = static_cast<int>(dt_save_coefficients / dt);
    steps_between_exports = static_cast<int>(dt_export_fields / dt);
    steps_saved = 0;
    total_steps_to_solve = static_cast<int>((final_time - initial_time) / dt);
    total_steps_to_save = static_cast<int>((final_time - initial_time) /
                                           dt_save_coefficients);
    total_steps_to_export = static_cast<int>((final_time - initial_time) /
        dt_export_fields);
    export_every_saved = static_cast<int>(total_steps_to_save /
                                          total_steps_to_export);

    // Ensure that every value makes sense and is positive
    if (dt <= 0 || dt_save_coefficients <= 0 || dt_export_fields <= 0
            || final_time <= initial_time)
    {
        throw std::invalid_argument("TimeManager: All time parameters must be positive and final_time must be greater than initial_time.");
    }

    if (steps_between_saves <= 0 || steps_between_exports <= 0)
    {
        throw std::invalid_argument("TimeManager: dt_save_coefficients and dt_export_fields must be greater than dt.");
    }

    if (total_steps_to_solve <= 0 || total_steps_to_save <= 0
            || total_steps_to_export <= 0)
    {
        throw std::invalid_argument("TimeManager: The total number of steps to solve, save, and export must be positive.");
    }

    if (export_every_saved <= 0)
    {
        std::cout << "The number of export every saved step is: " << export_every_saved
                  << std::endl;
        throw std::invalid_argument("TimeManager: The number of exports per saved step must be positive.");
    }

    std::cout << "###TIME - TimeManager initialized with the following parameters:" << std::endl;
    std::cout << "---Initial time: " << initial_time << std::endl;
    std::cout << "---Final time: " << final_time << std::endl;
    std::cout << "---Time step: " << dt << std::endl;
    std::cout << "---Time step to save coefficients: " << dt_save_coefficients << std::endl;
    std::cout << "---Time step to export fields: " << dt_export_fields << std::endl;
    std::cout << "---Total steps to solve: " << total_steps_to_solve << std::endl;
    std::cout << "---Total steps to save: " << total_steps_to_save << std::endl;
    std::cout << "---Total steps to export: " << total_steps_to_export << std::endl;
    std::cout << "---Steps between saves: " << steps_between_saves << std::endl;
    std::cout << "---Steps between exports: " << steps_between_exports << std::endl;
    std::cout << "---Exports every saved step: " << export_every_saved << std::endl;
}

TimeManager::TimeManager(const std::string& filename)
{
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        throw std::runtime_error("Could not open time settings file");
    }

    throw std::runtime_error("TimeManager constructor from file is not implemented yet. One day...");
}

bool TimeManager::isFinished()
{
    // Must solve even last step, when current_time == final_time - dt, so we check for >=
    return current_time >= final_time - TOLERANCE;
}

void TimeManager::advanceTime()
{
    if (current_time + dt > final_time - TOLERANCE)
    {
      actual_dt = final_time - current_time;
    }
    if (current_time + dt >= next_time_to_save - TOLERANCE)
    {
        save_coefficients = true;
        actual_dt = next_time_to_save - current_time;
    }
    else
    {
        save_coefficients = false;
        actual_dt = dt;
    }

    current_time += actual_dt;
    number_of_steps_taken++;
}

void TimeManager::advanceTime(double new_dt)
{
    dt = new_dt;
    advanceTime();
}

void TimeManager::updateAfterSave()
{
    steps_saved++;
    next_time_to_save = initial_time + (steps_saved + 1) * dt_save_coefficients;
}
