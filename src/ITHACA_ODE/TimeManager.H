#ifndef TIME_MANAGER_H
#define TIME_MANAGER_H

#include <string>

class TimeManager
{
private:
    double initial_time;
    double final_time;
    double dt;
    double actual_dt;

    double current_time;

    double dt_save_coefficients;
    double dt_export_fields;
    double next_time_to_save;

    int steps_between_saves;
    int steps_between_exports;

    int number_of_steps_taken;
    int steps_saved;

    int total_steps_to_solve;
    int total_steps_to_save;
    int total_steps_to_export;
    bool adaptive_time_stepping;
    bool save_coefficients = false;

public:
    // Constructors
    TimeManager() = default;

    // Constructor with specific time parameters
    TimeManager(double initial, double final, double timestep, double dt_save, double dt_export);

    // Constructor that reads required entries from a file
    TimeManager(const std::string& filename);

    // Destructor
    ~TimeManager() = default;

    // Allocation related methods - So you can call these to resize matrices before running the time loop
    int getNumberOfStepsToSolve() const { return total_steps_to_solve; }
    int getNumberOfStepsToSave() const { return total_steps_to_save; }
    int getNumberOfStepsToExport() const { return total_steps_to_export; }

    // Time management
    int getNumberOfStepsTaken() const { return number_of_steps_taken; }

    double getCurrentTime() const { return current_time; }
    double getTimeStep() const { return dt; }
    double getActualTimeStep() const { return actual_dt; }
    
    bool isFinished();
    bool shouldSaveCoefficients() const { return save_coefficients; }
    void advanceTime();
    void advanceTime(double new_dt);
    void updateAfterSave();
    
};

#endif // TIME_MANAGER_H
