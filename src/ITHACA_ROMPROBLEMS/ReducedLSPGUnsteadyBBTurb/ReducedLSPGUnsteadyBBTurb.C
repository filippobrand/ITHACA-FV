#include "ReducedLSPGUnsteadyBBTurb.H"
#include "Foam2Eigen.H" // Needed to convert OpenFOAM fields to Eigen vectors and matrices

ReducedLSPGUnsteadyBBTurb::ReducedLSPGUnsteadyBBTurb(UnsteadyBBTurb& problem, int argc, char* argv[])
: ReducedLSPG(argc, argv)
{
  #include "createFields.H"
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
    problem.ITHACAdict->lookupOrDefault<label>("NmodesUproj", 0),
    problem.ITHACAdict->lookupOrDefault<label>("NmodesPrghproj", 0),
    problem.ITHACAdict->lookupOrDefault<label>("NmodesTproj", 0),
    problem.ITHACAdict->lookupOrDefault<label>("NmodesNutproj", 0)
  };
  currentState_ = Eigen::VectorXd::Zero(numberOfModes_.velocity + numberOfModes_.pressure + numberOfModes_.temperature);
  currentNutCoeffs_ = Eigen::VectorXd::Zero(numberOfModes_.nut);
  problem_ = &problem;
  copyNutAvgFieldsToROMMesh();
  readEigenvalues();
}

void ReducedLSPGUnsteadyBBTurb::setTime(
  const scalar startTime,
  const scalar endTime,
  const scalar deltaT)
{
  _runTime->setTime(startTime, 0);
  _runTime->setEndTime(endTime);
  _runTime->setDeltaT(deltaT);
}

void ReducedLSPGUnsteadyBBTurb::solveOnline(const Eigen::MatrixXd& vel_now_BC,
                                             const Eigen::MatrixXd& temp_now_BC, 
                                             int startSnap)
{
  // The boundaryConditions_ object does not directly interact with the ROM,
  // except for the initialization of the reduced coeff. It simply stores the BCs values as a function of time
  // and can be interrogated to get the values at a given time
  boundaryConditions_ = BoundaryConditions{vel_now_BC, temp_now_BC, "linear"};
  boundaryConditions_.initializeReducedCoeffs(startSnap, currentState_, problem_,
    numberOfModes_.velocity, numberOfModes_.pressure, numberOfModes_.temperature
  );
  currentNutCoeffs_ = ITHACAutilities::getCoeffs(
    problem_->fluctNutfield[startSnap], problem_->nutmodes
  );
  currentNutAvgCoeffs_ = interpolateIDW(
    boundaryConditions_.getCurrentBCs().head(2)
  );

  GaussNewtonSettings gaussnewton_settings = GaussNewtonSettings
  {
    problem_->ITHACAdict->lookupOrDefault<int>("gaussNewtonMaxIter", 5),
    problem_->ITHACAdict->lookupOrDefault<float>("gaussNewtonTol", 1e-4)
  }; 
  
  Time& runTime = _runTime();
  // Maybe here we need: #include "initContinuityErrs.H" --- Check later
  #include "readTimeControls.H"
  while (runTime.run())
  {
    runTime++;
    Info << "Time = " << runTime.timeName() << nl << endl;
    boundaryConditions_.updateTimeDependentBC(runTime.time().value());
    for (int gnIter = 0; gnIter < gaussnewton_settings.maxIter; gnIter++)
    {
      Eigen::VectorXd residual = assembleResidual(currentState_, runTime);
      double residualNorm = residual.norm();
      if (residualNorm < gaussnewton_settings.tol)
      {
          Info << "Gauss-Newton converged at iteration " << gnIter
              << " with residual norm " << residualNorm << endl;
          break;
      }

      Eigen::MatrixXd jacobian = assembleJacobian(currentState_, residual, runTime);
      Eigen::VectorXd dq = jacobian.colPivHouseholderQr().solve(-residual);

      // Backtracking line search
      double alpha = 1.0;
      const int maxLineSearch = 5;
      bool accepted = false;
      for (int m = 0; m < maxLineSearch; m++)
      {
          Eigen::VectorXd trialState = currentState_ + alpha * dq;
          double trialNorm = assembleResidual(trialState, runTime).norm();
          if (trialNorm < residualNorm)
          {
              currentState_ = trialState;
              accepted = true;
              break;
          }
          alpha *= 0.5;
      }
      if (!accepted)
      {
          Info << "Line search failed to reduce residual at GN iter " << gnIter
              << " — accepting smallest trial step anyway" << endl;
          currentState_ += alpha * dq;   // last (smallest) alpha tried
      }
    }
    Info << "Final residual norm: " << assembleResidual(currentState_, runTime).norm() << endl;
    Info << "Current state vector: " << currentState_.transpose() << endl;
    runTime.write();
  }
}

Eigen::VectorXd ReducedLSPGUnsteadyBBTurb::assembleResidual(const Eigen::VectorXd& state, const Time& runTime)
{
  fvMesh& mesh = _mesh();
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

  reconstructReducedFields(state, U, p_rgh, T);
  phi = fvc::flux(U);
  interpolateNutCoeffs(nut, problem_->nutmodes, state);
  #include "UEqn.H"
  #include "TEqn.H"
  #include "pEqn.H"

  Eigen::VectorXd residual = Eigen::VectorXd::Zero(
    3 * U.internalField().size() + T.internalField().size() + p_rgh.internalField().size()
  );
  stackResiduals(UEqn, TEqn, p_rghEqn, residual);
  return residual;
}

void ReducedLSPGUnsteadyBBTurb::stackResiduals(
  const fvVectorMatrix& UEqn,
  const fvScalarMatrix& TEqn,
  const fvScalarMatrix& pEqn,
  Eigen::VectorXd& residual_)
{
  Field<vector> Ru = UEqn.residual();
  Field<scalar> RT = TEqn.residual();
  Field<scalar> Rp = pEqn.residual();

  Eigen::VectorXd Ru_eigen = Foam2Eigen::field2Eigen(Ru);
  Eigen::VectorXd RT_eigen = Foam2Eigen::field2Eigen(RT);
  Eigen::VectorXd Rp_eigen = Foam2Eigen::field2Eigen(Rp);
  
  residual_ << Ru_eigen, Rp_eigen, RT_eigen;
}

void ReducedLSPGUnsteadyBBTurb::reconstructReducedFields(
  const Eigen::VectorXd& state,
  volVectorField& velocity_field, 
  volScalarField& pressure_field,
  volScalarField& temperature_field)
{
  velocity_field = problem_->L_U_SUPmodes.reconstruct(velocity_field, state.head(numberOfModes_.velocity), "U");
  pressure_field = problem_->P_rghmodes.reconstruct(pressure_field, state.segment(numberOfModes_.velocity, numberOfModes_.pressure), "p_rgh");
  temperature_field = problem_->L_Tmodes.reconstruct(temperature_field, state.tail(numberOfModes_.temperature), "T");
}

Eigen::MatrixXd ReducedLSPGUnsteadyBBTurb::assembleJacobian(
  const Eigen::VectorXd& state,
  const Eigen::VectorXd& residual,
  const Time& runTime)
{
  int n = state.size();
  Eigen::MatrixXd J = Eigen::MatrixXd::Zero(residual.size(), n);
  double eps = 1e-5;
  double absoluteFloor_ = 1e-8; // To avoid division by zero in the finite difference approximation
  
  for (int i = 0; i < n; i++)
  {
    double h = std::sqrt(eps) * std::max({std::abs(state[i]), eigenvalues_(i), absoluteFloor_});
    Eigen::VectorXd state_plus = state;
    state_plus[i] += h;
    Eigen::VectorXd residual_plus = assembleResidual(state_plus, runTime);
    J.col(i) = (residual_plus - residual) / h; // forward difference approximation - Maybe use central
  }
  return J;
}

void ReducedLSPGUnsteadyBBTurb::interpolateNutCoeffs(volScalarField& nut_field, volScalarModes& nut_modes, const Eigen::VectorXd& state)
{
  // const label inputRBFSize = interpolationSettings_.dimInputRBF;
  const label inputRBFSize =
    interpolationSettings_.dimInputRBF > 0
  ? interpolationSettings_.dimInputRBF
  : numberOfModes_.velocity;
  Eigen::VectorXd velocity_coeffs = state.head(inputRBFSize);
  if (interpolationSettings_.derivative)
  {
    FatalErrorInFunction << "Derivative interpolation for nut coefficients is not implemented yet." << exit(FatalError);
  }
  else
  {
    for (int j=0; j < numberOfModes_.nut; j++)
    {
      currentNutCoeffs_(j) = problem_->rbfSplines[j]->predict(velocity_coeffs);
    }
  }
  nut_field = nut_modes.reconstruct(nut_field, currentNutCoeffs_, "nut");
  volScalarField nutAvg(
    IOobject(
        "nutAvgRec",
        nut_field.time().timeName(),
        nut_field.mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    nut_field.mesh(),
    dimensionedScalar("zero", nut_modes[0].dimensions(), 0.0));

  for (int k = 0; k < currentNutAvgCoeffs_.size(); k++)
  {
      nutAvg += currentNutAvgCoeffs_(k) * avgNutFields_[k];
  }
  nut_field += nutAvg;
}

Eigen::VectorXd ReducedLSPGUnsteadyBBTurb::interpolateIDW(const Eigen::VectorXd&
        input_parameters)
{
    const label n_samples = interpolationSettings_.avgTermOfflineSamples;
    Eigen::VectorXd weights(n_samples);

    for (label i = 0; i < n_samples; i++)
    {
        weights(i) = 1.0 / ((input_parameters - problem_->mu.col(i)).norm() + 1e-10);
    }

    const double weightSum = weights.sum();
    if (weightSum > 1e-10)
    {
        weights /= weightSum;
    }
    else
    {
        weights.setConstant(1.0 / n_samples);
    }

    // Try with a diagonal matrix with ones
    const Eigen::MatrixXd nutAvgCoeffs = Eigen::MatrixXd::Identity(n_samples, n_samples);

    Eigen::VectorXd interpolatedNutCoeffs = nutAvgCoeffs * weights;

    for (int j = 0; j < interpolatedNutCoeffs.size(); j++)
    {
        if (std::abs(interpolatedNutCoeffs(j)) < 1e-6)
        {
            interpolatedNutCoeffs(j) = 0.0;
        }
    }
    return interpolatedNutCoeffs;
}

void ReducedLSPGUnsteadyBBTurb::copyNutAvgFieldsToROMMesh()
{
  int num_of_avgNutFields = problem_->avgNutfield.size();
  avgNutFields_.setSize(0);
  for (int i = 0; i < num_of_avgNutFields; i++)
  {
    volScalarField avgNutFieldOnROMMesh(
      IOobject(
        "avgNut" + Foam::name(i),
        _runTime().timeName(),
        _mesh(),
        IOobject::NO_READ,
        IOobject::NO_WRITE),
      _mesh(),
      dimensionedScalar("zero", problem_->avgNutfield[i].dimensions(), 0.0));
    
    avgNutFieldOnROMMesh.primitiveFieldRef() =
    problem_->avgNutfield[i].primitiveField();

    forAll(avgNutFieldOnROMMesh.boundaryFieldRef(), patchI)
    {
        scalarField& destination =
            avgNutFieldOnROMMesh.boundaryFieldRef()[patchI];

        const scalarField& source =
            problem_->avgNutfield[i].boundaryField()[patchI];

        forAll(destination, faceI)
        {
            destination[faceI] = source[faceI];
        }
    }
    avgNutFields_.append(avgNutFieldOnROMMesh.clone());
  }
}

void ReducedLSPGUnsteadyBBTurb::readEigenvalues()
{
    // We read the eigenvalues from the files in the ./ITHACAoutput/POD folder
    // These files follow this form:
    // %%MatrixMarket matrix array real general
    // 2720 1
    // 0.34428033233490540344
    // 0.12360804880946965612
    // 0.08097235508005383442
    // ...

    Eigen::VectorXd uEigenvalues_, pEigenvalues_, tEigenvalues_, nutEigenvalues_;
    uEigenvalues_ = Eigen::VectorXd::Zero(numberOfModes_.velocity);
    pEigenvalues_ = Eigen::VectorXd::Zero(numberOfModes_.pressure);
    tEigenvalues_ = Eigen::VectorXd::Zero(numberOfModes_.temperature);
    nutEigenvalues_ = Eigen::VectorXd::Zero(numberOfModes_.nut);
    if (Pstream::master())
    {
      std::ifstream uFile("ITHACAoutput/POD/Eigenvalues_U");
      M_Assert(uFile.is_open(),
              "Could not open file ITHACAoutput/POD/Eigenvalues_U. Please make sure the file exists and is readable.");
      std::string line;
      std::getline(uFile, line); // Ignore first line
      std::getline(uFile, line); // Ignore second line
      uEigenvalues_.resize(numberOfModes_.velocity);

      for (int i = 0; i < numberOfModes_.velocity; i++)
      {
          std::getline(uFile, line);
          uEigenvalues_(i) = std::stod(line);
      }

      std::ifstream pFile("ITHACAoutput/POD/Eigenvalues_p_rgh");
      M_Assert(pFile.is_open(),
              "Could not open file ITHACAoutput/POD/Eigenvalues_p_rgh. Please make sure the file exists and is readable.");
      std::getline(pFile, line); // Ignore first line
      std::getline(pFile, line); // Ignore second line
      pEigenvalues_.resize(numberOfModes_.pressure);

      for (int i = 0; i < numberOfModes_.pressure; i++)
      {
          std::getline(pFile, line);
          pEigenvalues_(i) = std::stod(line);
      }

      std::ifstream tFile("ITHACAoutput/POD/Eigenvalues_T");
      M_Assert(tFile.is_open(),
              "Could not open file ITHACAoutput/POD/Eigenvalues_T. Please make sure the file exists and is readable.");
      std::getline(tFile, line); // Ignore first line
      std::getline(tFile, line); // Ignore second line
      tEigenvalues_.resize(numberOfModes_.temperature);

      for (int i = 0; i < numberOfModes_.temperature; i++)
      {
          std::getline(tFile, line);
          tEigenvalues_(i) = std::stod(line);
      }

      std::ifstream nutFile("ITHACAoutput/POD/Eigenvalues_fluctNut");
      M_Assert(nutFile.is_open(),
              "Could not open file ITHACAoutput/POD/Eigenvalues_fluctNut. Please make sure the file exists and is readable.");
      std::getline(nutFile, line); // Ignore first line
      std::getline(nutFile, line); // Ignore second line
      nutEigenvalues_.resize(numberOfModes_.nut);

      for (int i = 0; i < numberOfModes_.nut; i++)
      {
          std::getline(nutFile, line);
          nutEigenvalues_(i) = std::stod(line);
      }
    }

    if (Pstream::parRun())
    {
        reduce(uEigenvalues_, sumOp<Eigen::VectorXd>());
        reduce(pEigenvalues_, sumOp<Eigen::VectorXd>());
        reduce(tEigenvalues_, sumOp<Eigen::VectorXd>());
        reduce(nutEigenvalues_, sumOp<Eigen::VectorXd>());
    }    

    Info << "### EIGS - Velocity eigenvalues: " << uEigenvalues_.transpose() <<
         endl;
    Info << "### EIGS - Pressure eigenvalues: " << pEigenvalues_.transpose() <<
         endl;
    Info << "### EIGS - Temperature eigenvalues: " << tEigenvalues_.transpose() <<
         endl;
    Info << "### EIGS - FluctNut eigenvalues: " << nutEigenvalues_.transpose() <<
         endl;
        
    eigenvalues_ = Eigen::VectorXd::Zero(
        numberOfModes_.velocity + numberOfModes_.pressure + numberOfModes_.temperature);
    eigenvalues_ << uEigenvalues_, pEigenvalues_, tEigenvalues_;
}
