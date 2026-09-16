#include "ReducedLSPGUnsteadyBBTurb.H"
#include "Foam2Eigen.H" // Needed to convert OpenFOAM fields to Eigen vectors and matrices

ReducedLSPGUnsteadyBBTurb::ReducedLSPGUnsteadyBBTurb(UnsteadyBBTurb& problem, int argc, char* argv[])
: ReducedLSPG(argc, argv)
{
  #include "createFields.H"
  romSettings_ = ROMSettings
  {
      problem.ITHACAdict->lookupOrDefault<word>("method", "PPE"),
      problem.ITHACAdict->lookupOrDefault<word>("bcMethod", "Gunzburger"),
      problem.ITHACAdict->lookupOrDefault<bool>("hasMonitors", false),
      problem.ITHACAdict->lookupOrDefault<bool>("timeDependentBC", false),
      problem.ITHACAdict->lookupOrDefault<bool>("normalizeWithEigenvalues", false),
      (int)problem.inletIndex.rows(),
      (int)problem.inletIndexT.rows(),
      problem.ITHACAdict->lookupOrDefault<word>("solverODE", "BDF2")
  };
  interpolationSettings_ = InterpolationSettings
  {
      problem.ITHACAdict->lookupOrDefault<int>("firstRBFIndex", 0),
      problem.dimA,
      problem.ITHACAdict->lookupOrDefault<bool>("derivativeInRBF", false),
      problem.ITHACAdict->lookupOrDefault<label>("dimInputRBF", 0),
      problem.mu.cols()
  };
  // Create the number of modes object
  numberOfModes_ = NumberOfModes
  {
    problem.ITHACAdict->lookupOrDefault<label>("NUmodes", 0),
    problem.ITHACAdict->lookupOrDefault<label>("NPmodes", 0),
    problem.ITHACAdict->lookupOrDefault<label>("NTmodes", 0),
    problem.ITHACAdict->lookupOrDefault<label>("NNutModes", 0)
  };
  currentState_ = Eigen::VectorXd::Zero(numberOfModes_.velocity + numberOfModes_.pressure + numberOfModes_.temperature);
  problem_ = &problem;
}

void ReducedLSPGUnsteadyBBTurb::setTime(const scalar startTime, const scalar endTime, const scalar deltaT)
{
  _runTime->setTime(startTime, 0);
  _runTime->setEndTime(endTime);
  _runTime->setDeltaT(deltaT);
}

void ReducedLSPGUnsteadyBBTurb::solveOnline(const Eigen::MatrixXd& vel_now_BC,
                                             const Eigen::MatrixXd& temp_now_BC, 
                                             int startSnap)
{
  boundaryConditions_ = BoundaryConditions
  {
    vel_now_BC,
    temp_now_BC,
    "linear"
  };
  boundaryConditions_.initializeReducedCoeffs(startSnap, currentState_, problem_,
    numberOfModes_.velocity, numberOfModes_.pressure, numberOfModes_.temperature
  );

  /* Here we create some references to OpenFOAM objects.
  This is for convenience and compatibility with OpenFOAM code.*/
  Time& runTime = _runTime();
  fvMesh& mesh = _mesh();
  #include "initContinuityErrs.H"
  fv::options& fvOptions = _fvOptions();
  singlePhaseTransportModel& laminarTransport = _laminarTransport();
  volScalarField& p = _p();
  volVectorField& U = _U();
  volScalarField& p_rgh = _p_rgh();
  volScalarField& T = _T();
  volScalarField& nut = _nut();
  volScalarField& alphat = _alphat();
  volScalarField& rhok = _rhok();
  volScalarField& gh = _gh();
  surfaceScalarField& ghf = _ghf();
  surfaceScalarField& phi = _phi();
  pimpleControl& pimple = _pimple();
  IOMRFZoneList& MRF = _MRF();
  dimensionedScalar& beta = _beta();
  dimensionedScalar& TRef = _TRef();
  dimensionedScalar& Pr = _Pr();
  dimensionedScalar& Prt = _Prt();
  
  #include "readTimeControls.H"
  Eigen::VectorXd residual = Eigen::VectorXd::Zero(U.internalField().size() + T.internalField().size() + p_rgh.internalField().size());
  
  Info << "Starting the LSPG time loop" << endl;
  while (runTime.run())
  {
    runTime++;
    boundaryConditions_.updateTimeDependentBC(runTime.time().value());
    U = problem_->L_U_SUPmodes.reconstruct(U, currentState_.head(numberOfModes_.velocity), "U");
    phi = fvc::flux(U);
    p_rgh = problem_->P_rghmodes.reconstruct(p_rgh, currentState_.segment(numberOfModes_.velocity, numberOfModes_.pressure), "p_rgh");
    T = problem_->L_Tmodes.reconstruct(T, currentState_.tail(numberOfModes_.temperature), "T");
    // TODO: Add nut RBF interpolation here
    
    // Here we create the residual for the LSPG problem leveraging the OpenFOAM using the reduced fields
    #include "UEqn.H" // Assembles UEqn
    #include "TEqn.H"
    #include "pEqn.H"
    Field<vector> Ru = UEqn.residual();
    Field<scalar> RT = TEqn.residual();
    Field<scalar> Rp = p_rghEqn.residual();

    Eigen::VectorXd Ru_eigen = Foam2Eigen::field2Eigen(Ru);
    Eigen::VectorXd RT_eigen = Foam2Eigen::field2Eigen(RT);
    Eigen::VectorXd Rp_eigen = Foam2Eigen::field2Eigen(Rp);
  
    stackResiduals(Ru, RT, Rp, residual);
  }
}

void ReducedLSPGUnsteadyBBTurb::stackResiduals(const Field<vector>& Ru, const Field<scalar>& RT, const Field<scalar>& Rp, Eigen::VectorXd& residual)
{
  Eigen::VectorXd Ru_eigen = Foam2Eigen::field2Eigen(Ru);
  Eigen::VectorXd RT_eigen = Foam2Eigen::field2Eigen(RT);
  Eigen::VectorXd Rp_eigen = Foam2Eigen::field2Eigen(Rp);
  residual << Ru_eigen, RT_eigen, Rp_eigen;
}

void ReducedLSPGUnsteadyBBTurb::reconstructReducedFields()
{
  // Reconstruct the reduced solution fields from the currentState_ vector
  // Ur_ = problem_->L_U_SUPmodes.reconstruct(Ur_, currentState_.head(numberOfModes_.velocity), "Ur");
  // Pr_ = problem_->P_rghmodes.reconstruct(Pr_, currentState_.segment(numberOfModes_.velocity, numberOfModes_.pressure), "Pr");
  // Tr_ = problem_->L_Tmodes.reconstruct(Tr_, currentState_.segment(numberOfModes_.velocity + numberOfModes_.pressure, numberOfModes_.temperature), "Tr");
  // Nutr_ = problem_->nutmodes.reconstruct(Nutr_, currentState_.tail(numberOfModes_.nut), "Nutr");
}
