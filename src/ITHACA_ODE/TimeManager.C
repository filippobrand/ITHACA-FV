#include <iostream>
#include <fstream>
#include "TimeManager.H"
// Here I try minimal dependence on ITHACA/OpenFOAM, as I steal this code
// elsewhere

const double TOLERANCE = 1e-10; // Tolerance for floating-point comparisons

TimeManager::TimeManager(double initial, double final, double timestep, double dt_save, double dt_export): initial_time(initial), final_time(final), dt(timestep), dt_save_coefficients(dt_save), dt_export_fields(dt_export), current_time(initial), adaptive_time_stepping(false)
{
    number_of_steps_taken = 0;
    next_time_to_save = initial_time + dt_save_coefficients;
    actual_dt = dt;
    steps_between_saves = static_cast<int>(dt_save_coefficients / dt);
    steps_between_exports = static_cast<int>(dt_export_fields / dt);
    steps_saved = 0;
    total_steps_to_solve = static_cast<int>((final_time - initial_time) / dt);
    total_steps_to_save = static_cast<int>((final_time - initial_time) / dt_save_coefficients);
    total_steps_to_export = static_cast<int>((final_time - initial_time) / dt_export_fields);
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
{ // Must solve even last step, when current_time == final_time - dt, so we check for >=
    return current_time >= final_time + TOLERANCE;
}

void TimeManager::advanceTime()
{
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