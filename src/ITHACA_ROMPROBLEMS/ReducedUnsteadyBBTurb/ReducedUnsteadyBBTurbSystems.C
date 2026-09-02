#include "ReducedUnsteadyBBTurbSystems.H"

void SystemPPEGunzburger::evaluateRHS(
    const Eigen::VectorXd& state,
    Eigen::VectorXd& rhs, double t) const
{
    Info << "RHS not implemented for SystemPPEGunzburger" << endl;
}

void SystemPPEGunzburger::evaluateResidual(
    const Eigen::VectorXd& state,
    const Eigen::VectorXd& state_dot,
    Eigen::VectorXd& residual,
    double t) const
{
    residual.setZero(state.size());
    // Cache dimensions for clarity and safety
    const int nVel = rom.expressionModes_.velocity;
    const int nPress = rom.expressionModes_.pressure;
    const int nTemp = rom.expressionModes_.temperature;
    // Test modes
    const int tVel = rom.testModes_.velocity;
    const int tPress = rom.testModes_.pressure;
    const int tTemp = rom.testModes_.temperature;

    // Extract and scale state vectors avoiding unnecessary copies where possible
    Eigen::VectorXd a_tmp = state.head(nVel).cwiseQuotient(rom.scalingVector_.head(nVel));
    Eigen::VectorXd b_tmp = state.segment(nVel, nPress).cwiseQuotient(rom.scalingVector_.segment(nVel, nPress));
    Eigen::VectorXd c_tmp = state.tail(nTemp).cwiseQuotient(rom.scalingVector_.tail(nTemp));

    Eigen::VectorXd a_dot = state_dot.head(nVel).cwiseQuotient(rom.scalingVector_.head(nVel));
    Eigen::VectorXd c_dot = state_dot.tail(nTemp).cwiseQuotient(rom.scalingVector_.tail(nTemp));

    
    // Momentum equation terms
    const Eigen::VectorXd M11 = commonMat_.BTotal * a_tmp *
                                rom.fluidProperties.nu; // Viscous term
    const Eigen::VectorXd M2 = commonMat_.K * b_tmp; // Pressure gradient
    const Eigen::VectorXd M5 = commonMat_.M * a_dot; // Mass term
    const Eigen::VectorXd M10 = commonMat_.H * c_tmp; // Buoyancy term
    // Temperature equation terms
    const Eigen::VectorXd M6 = commonMat_.Y * c_tmp * (rom.fluidProperties.nu /
        rom.fluidProperties.Pr); // Diffusive term
    const Eigen::VectorXd M8 = commonMat_.W * c_dot; // Mass term
    // PPE equation terms
    const Eigen::VectorXd M3 = ppeMat_.D * b_tmp; // Laplacian term
    const Eigen::VectorXd M12 = ppeMat_.HP * c_tmp; // Buoyancy term
    const Eigen::VectorXd M7 = ppeMat_.nuBC * a_tmp *
                               rom.fluidProperties.nu; // BC term
    const Eigen::VectorXd M13 = ppeMat_.timedepBC * a_dot; // Time-dependent BC
    // Gunzburger matrices 
    const Eigen::MatrixXd gunzburgerBCProduct = a_tmp.transpose() *
        gunzMat_.bcVelMat;
    const Eigen::MatrixXd gunzburgerBCProductTemp = c_tmp.transpose() *
        gunzMat_.bcTempMat;

    // Assemble the residual vector - Momentum equation
    for (int i = 0; i < tVel; i++)
    {
        scalar cc = a_tmp.transpose() * Eigen::SliceFromTensor(commonMat_.C, 0,
            i) * a_tmp;
        scalar ct = rom.nutCurrentCoeffs_.transpose() * Eigen::SliceFromTensor(
                        commonMat_.CTotal, 0, i) * a_tmp;
        scalar caveraged = rom.nutAvgCurrentCoeffs_.transpose() *
                           Eigen::SliceFromTensor(commonMat_.CTotalAve, 0, i) * a_tmp;
        residual(i) = -M5(i) + M11(i) - cc - M10(i) - M2(i) + ct + caveraged;
    }

    // Fill the BC condition equation for the gunzburger method
    for (int i = 0; i < rom.romSettings_.N_BC; i++)
    {
        residual(tVel + i) = gunzburgerBCProduct(0, i) - bc_.currentVelocityBC(i);
    }

    // Fill the Pressure poisson equation
    for (int j = 0; j < tPress; j++)
    {
        int idx = nVel + j;
        scalar gg = a_tmp.transpose() * Eigen::SliceFromTensor(ppeMat_.G, 0, j) * a_tmp;
        scalar turbPPE = rom.nutCurrentCoeffs_.transpose() * Eigen::SliceFromTensor(
                             ppeMat_.CTotalPPEFluct, 0, j) * a_tmp;
        scalar turbPPE_ave = rom.nutAvgCurrentCoeffs_.transpose() *
                             Eigen::SliceFromTensor(ppeMat_.CTotalPPEAve, 0, j) * a_tmp;
        residual(idx) = M3(j) + gg + M12(j) - M7(j) - turbPPE - turbPPE_ave;
    }

    // Fill the Temperature equation
    const int nTempOffset = nVel + nPress;
    for (int j = 0; j < tTemp; j++)
    {
        const int idx = nTempOffset + j;
        scalar qq = a_tmp.transpose() * Eigen::SliceFromTensor(commonMat_.Q, 0,
            j) * c_tmp;
        scalar qt = rom.nutCurrentCoeffs_.transpose() * Eigen::SliceFromTensor(
                        commonMat_.YTurb, 0, j) * c_tmp;
        scalar qaveraged = rom.nutAvgCurrentCoeffs_.transpose() *
                           Eigen::SliceFromTensor(commonMat_.AveYTurb, 0, j) * c_tmp;
        residual(idx) = -M8(j) + M6(j) - qq + qt / rom.fluidProperties.Pr_t +
                        qaveraged / rom.fluidProperties.Pr_t;
    }

    const int tempBCOffset = nTempOffset + tTemp;
    for (int j = 0; j < rom.romSettings_.N_BC_t; j++)
    {
        residual(tempBCOffset + j) = gunzburgerBCProductTemp(0, j) - bc_.currentTemperatureBC(j);
    }
}


void SystemSupremizerGunzburger::evaluateRHS(
    const Eigen::VectorXd& state,
    Eigen::VectorXd& rhs, double t) const
{
    Info << "RHS not implemented for ReducedUnsteadyBBTurbSystem" << endl;
}


void SystemSupremizerGunzburger::evaluateResidual(
    const Eigen::VectorXd& state,
    const Eigen::VectorXd& state_dot,
    Eigen::VectorXd& residual,
    double t) const
{
    residual.setZero(state.size());
    Eigen::VectorXd a_tmp = state.head(rom.expressionModes_.velocity);
    Eigen::VectorXd b_tmp = state.segment(rom.expressionModes_.velocity, rom.expressionModes_.pressure);
    Eigen::VectorXd c_tmp = state.tail(rom.expressionModes_.temperature);
    Eigen::VectorXd a_dot = state_dot.head(rom.expressionModes_.velocity);
    Eigen::VectorXd c_dot = state_dot.tail(rom.expressionModes_.temperature);

    a_tmp = a_tmp.cwiseQuotient(rom.scalingVector_.head(rom.expressionModes_.velocity));
    b_tmp = b_tmp.cwiseQuotient(rom.scalingVector_.segment(rom.expressionModes_.velocity, rom.expressionModes_.pressure));
    c_tmp = c_tmp.cwiseQuotient(rom.scalingVector_.tail(rom.expressionModes_.temperature));
    a_dot = a_dot.cwiseQuotient(rom.scalingVector_.head(rom.expressionModes_.velocity));
    c_dot = c_dot.cwiseQuotient(rom.scalingVector_.tail(rom.expressionModes_.temperature));
    
    // Momentum equation terms
    const Eigen::VectorXd M11 = commonMat_.BTotal * a_tmp *
                                rom.fluidProperties.nu; // Viscous term
    const Eigen::VectorXd M2 = commonMat_.K * b_tmp; // Pressure gradient
    const Eigen::VectorXd M5 = commonMat_.M * a_dot; // Mass term
    const Eigen::VectorXd M10 = commonMat_.H * c_tmp; // Buoyancy term
    // Temperature equation terms
    const Eigen::VectorXd M6 = commonMat_.Y * c_tmp * (rom.fluidProperties.nu /
        rom.fluidProperties.Pr); // Diffusive term
    const Eigen::VectorXd M8 = commonMat_.W * c_dot; // Mass term
    // Supremizer term. It's just one
    const Eigen::VectorXd M14 = supMat_.P * a_tmp; // Supremizer term
      
    // Gunzburger matrices 
    const Eigen::MatrixXd gunzburgerBCProduct = a_tmp.transpose() *
        gunzMat_.bcVelMat;
    const Eigen::MatrixXd gunzburgerBCProductTemp = c_tmp.transpose() *
        gunzMat_.bcTempMat;

    // Assemble the residual vector - Momentum equation
    for (int i = 0; i < rom.testModes_.velocity; i++)
    {
        scalar cc = a_tmp.transpose() * Eigen::SliceFromTensor(commonMat_.C, 0,
            i) * a_tmp;
        scalar ct = rom.nutCurrentCoeffs_.transpose() * Eigen::SliceFromTensor(
                        commonMat_.CTotal, 0, i) * a_tmp;
        scalar caveraged = rom.nutAvgCurrentCoeffs_.transpose() *
                           Eigen::SliceFromTensor(commonMat_.CTotalAve, 0, i) * a_tmp;
        residual(i) = -M5(i) + M11(i) - cc - M10(i) - M2(i) + ct + caveraged;
    }

    // Fill the BC condition equation for the gunzburger method
    for (int i = 0; i < rom.romSettings_.N_BC; i++)
    {
        int idx = rom.testModes_.velocity + i;
        residual(idx) = gunzburgerBCProduct(0, i) - bc_.currentVelocityBC(i);
    }

    // Fill the Pressure poisson equation
    for (int j = 0; j < rom.testModes_.pressure; j++)
    {
        int idx = rom.expressionModes_.velocity + j;
        residual(idx) = M14(j);
    }

    // Fill the Temperature equation
    for (int j = 0; j < rom.testModes_.temperature; j++)
    {
        int idx = j + rom.expressionModes_.velocity + rom.expressionModes_.pressure;
        scalar qq = a_tmp.transpose() * Eigen::SliceFromTensor(commonMat_.Q, 0,
            j) * c_tmp;
        scalar qt = rom.nutCurrentCoeffs_.transpose() * Eigen::SliceFromTensor(
                        commonMat_.YTurb, 0, j) * c_tmp;
        scalar qaveraged = rom.nutAvgCurrentCoeffs_.transpose() *
                           Eigen::SliceFromTensor(commonMat_.AveYTurb, 0, j) * c_tmp;
        residual(idx) = -M8(j) + M6(j) - qq + qt / rom.fluidProperties.Pr_t +
                        qaveraged / rom.fluidProperties.Pr_t;
    }

    for (int j = 0; j < rom.romSettings_.N_BC_t; j++)
    {
        int idx = rom.expressionModes_.velocity 
                + rom.expressionModes_.pressure 
                + rom.testModes_.temperature 
                + j;
        residual(idx) = gunzburgerBCProductTemp(0, j) - bc_.currentTemperatureBC(j);
    }
}