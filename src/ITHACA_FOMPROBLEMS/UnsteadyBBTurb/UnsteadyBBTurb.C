/*---------------------------------------------------------------------------*\
     ██╗████████╗██╗  ██╗ █████╗  ██████╗ █████╗       ███████╗██╗   ██╗
     ██║╚══██╔══╝██║  ██║██╔══██╗██╔════╝██╔══██╗      ██╔════╝██║   ██║
     ██║   ██║   ███████║███████║██║     ███████║█████╗█████╗  ██║   ██║
     ██║   ██║   ██╔══██║██╔══██║██║     ██╔══██║╚════╝██╔══╝  ╚██╗ ██╔╝
     ██║   ██║   ██║  ██║██║  ██║╚██████╗██║  ██║      ██║      ╚████╔╝
     ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝      ╚═╝       ╚═══╝

 * In real Time Highly Advanced Computational Applications for Finite Volumes
 * Copyright (C) 2017 by the ITHACA-FV authors
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
/// Source file of the UnsteadyBBTurb class.

#include "UnsteadyBBTurb.H"
#include "viscosityModel.H"
#include "alphatJayatillekeWallFunctionFvPatchScalarField.H" // Used to implement BCs
#include "calculatedFvPatchField.H" // Used to implement BCs
#include <cmath>

// * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * * //

// Constructor
UnsteadyBBTurb::UnsteadyBBTurb() { }

UnsteadyBBTurb::UnsteadyBBTurb(int argc, char* argv[])
{
    _args = autoPtr<argList>(
        new argList(argc, argv));

    if (!_args->checkRootCase())
    {
        Foam::FatalError.exit();
    }

    argList& args = _args();
#include "createTime.H"
#include "createMesh.H"
    _pimple = autoPtr<pimpleControl>(
        new pimpleControl(
            mesh));
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "createFields.H"
#pragma GCC diagnostic pop
#include "createFvOptions.H"
    ITHACAdict = new IOdictionary(
        IOobject(
            "ITHACAdict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE));
    // Are these calls necessary? Some of these also present in the constructor of unsteadyNS. Check please.
    method = ITHACAdict->lookupOrDefault<word>("method", "supremizer");
    bcMethod = ITHACAdict->lookupOrDefault<word>("bcMethod", "lift");
    timedepbcMethod = ITHACAdict->lookupOrDefault<word>("timedepbcMethod", "no");

    M_Assert(!(bcMethod == "lift" && timedepbcMethod == "yes"),
        "The lifting method is not compatible with time-dependent BCs. Please choose the penalty method.");
    M_Assert(method == "supremizer" || method == "PPE",
        "A method must be set to either 'supremizer' or 'PPE' in ITHACAdict");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
        "The BC method must be set to lift or penalty in ITHACAdict");
    M_Assert(timedepbcMethod == "yes" || timedepbcMethod == "no",
        "The BC method can be set to yes or no");
    turbulence->validate();
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
}


// * * * * * * * * * * * * * * Full Order Methods * * * * * * * * * * * * * * //
#include "fvCFD.H"

// Method to performa a truthSolve
void UnsteadyBBTurb::truthSolve(List<scalar> mu_now, label nSample)
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    Info << "Created mesh and runTime references." << nl << endl;
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
    // volScalarField& Sourceterm = _Sourceterm();
    surfaceScalarField& ghf = _ghf();
    surfaceScalarField& phi = _phi();
    pimpleControl& pimple = _pimple();
    IOMRFZoneList& MRF = _MRF();
    dimensionedScalar& beta = _beta();
    dimensionedScalar& TRef = _TRef();
    dimensionedScalar& Pr = _Pr();
    dimensionedScalar& Prt = _Prt();

    // Here we handle the time settings. In other instances in ITHACA-FV,
    // runTime.setTime() uses as argument (Times[1],0). Don't know the reason, ask for clarification.
    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(startTime, 0);
    runTime.setDeltaT(timeStep);

    nextWrite = startTime + writeEvery;
    label nSavedTimesteps = (finalTime - startTime) / writeEvery;
    M_Assert(timeSnapshots.size() > nSample,
        "The timeSnapshots list does not have enough space for the current sample index.");
    timeSnapshots[nSample].resize(nSavedTimesteps);
    label stepCounter = 0;

    Info << "Starting time loop." << nl << endl;
    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime++;
        // runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;
        while (pimple.loop())
        {
#include "UEqn.H"
#include "TEqn.H"
            while (pimple.correct())
            {
#include "pEqn.H"
            }
            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        if (checkWrite(runTime))
        {
            // nut = turbulence->nut();
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p_rgh, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
            // ITHACAstream::exportSolution(nut, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" + runTime.timeName());

            timeSnapshots[nSample](stepCounter) = runTime.value();
            stepCounter++;

            Ufield.append(U.clone());
            Prghfield.append(p_rgh.clone());
            Tfield.append(T.clone());
            // Nutfield.append(nut.clone());
            nextWrite += writeEvery;
            writeMu(mu_now);
            counter++;
        }
    }
}

void UnsteadyBBTurb::truthSolve(fileName folder)
{
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

    instantList Times = runTime.times();
    runTime.setEndTime(finalTime);
    runTime.setTime(startTime, 0);
    runTime.setDeltaT(timeStep);

    nextWrite = startTime + writeEvery;

    while (runTime.run())
    {
#include "readTimeControls.H"
#include "CourantNo.H"
#include "setDeltaT.H"
        runTime++;
        // runTime.setEndTime(finalTime + timeStep);
        Info << "Time = " << runTime.timeName() << nl << endl;
        while (pimple.loop())
        {
#include "UEqn.H"
#include "TEqn.H"
            while (pimple.correct())
            {
#include "pEqn.H"
            }
            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
        if (checkWrite(runTime))
        {
            nut = turbulence->nut();
            ITHACAstream::exportSolution(U, name(counter), folder);
            ITHACAstream::exportSolution(p, name(counter), folder);
            ITHACAstream::exportSolution(p_rgh, name(counter), folder);
            ITHACAstream::exportSolution(T, name(counter), folder);
            ITHACAstream::exportSolution(nut, name(counter), folder);
            std::ofstream of(folder + "/" + name(counter) + "/" + runTime.timeName());

            Ufield.append(U.clone());
            Prghfield.append(p_rgh.clone());
            Tfield.append(T.clone());
            Nutfield.append(nut.clone());
            nextWrite += writeEvery;
            counter++;
        }
    }
}

void UnsteadyBBTurb::solvesupremizer(word type)
{
    M_Assert(type == "modes" || type == "snapshots",
        "You must specify the variable type with either snapshots or modes");
    PtrList<volScalarField> P_sup;

    if (type == "snapshots")
    {
        P_sup = Prghfield;
    } else
    {
        P_sup = P_rghmodes.toPtrList();
    }

    if (supex == 1)
    {
        volVectorField U = _U();
        volVectorField Usup(
            IOobject(
                "Usup",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero));

        if (type == "snapshots")
        {
            ITHACAstream::read_fields(supfield, Usup, "./ITHACAoutput/supfield/");
        } else
        {
            ITHACAstream::read_fields(supmodes, Usup, "./ITHACAoutput/supremizer/");
        }
    } else
    {
        volVectorField U = _U();
        volVectorField Usup(
            IOobject(
                "Usup",
                U.time().timeName(),
                U.mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            U.mesh(),
            dimensionedVector("zero", U.dimensions(), vector::zero));
        dimensionedScalar nu_fake(
            "nu_fake",
            dimensionSet(0, 2, -1, 0, 0, 0, 0),
            scalar(1));
        Vector<double> v(0, 0, 0);

        for (label i = 0; i < Usup.boundaryField().size(); i++)
        {
            if (Usup.boundaryField()[i].type() != "processor")
            {
                ITHACAutilities::changeBCtype(Usup, "fixedValue", i);
                assignBC(Usup, i, v);
                assignIF(Usup, v);
            }
        }

        if (type == "snapshots")
        {
            for (label i = 0; i < P_sup.size(); i++)
            {
                fvVectorMatrix u_sup_eqn(
                    -fvm::laplacian(nu_fake, Usup));
                solve(
                    u_sup_eqn == fvc::grad(P_sup[i]));
                supfield.append(Usup.clone());
                ITHACAstream::exportSolution(Usup, name(i + 1), "./ITHACAoutput/supfield/");
            }
            ITHACAutilities::createSymLink("./ITHACAoutput/supfield");
        } else
        {
            for (label i = 0; i < P_sup.size(); i++)
            {
                fvVectorMatrix u_sup_eqn(
                    -fvm::laplacian(nu_fake, Usup));
                solve(
                    u_sup_eqn == fvc::grad(P_sup[i]));
                supmodes.append(Usup.clone());
                ITHACAstream::exportSolution(Usup, name(i + 1), "./ITHACAoutput/supremizer/");
            }
            ITHACAutilities::createSymLink("./ITHACAoutput/supremizer");
        }
    }
}

// * * * * * * * * * * * * * * Projection Methods * * * * * * * * * * * * * * //

void UnsteadyBBTurb::projectSUP(fileName folder, label NU, label NPrgh, label NT,
    label NSUP, label Nnut, float RBFShapeParameter)
{
    NUmodes = NU;
    NTmodes = NT;
    NSUPmodes = NSUP;
    NPrghmodes = NPrgh;
    Nnutmodes = Nnut;
    L_U_SUPmodes.resize(0);
    if (liftfield.size() != 0)
    {
        for (label k = 0; k < liftfield.size(); k++)
        {
            L_U_SUPmodes.append(liftfield[k].clone());
        }
    }
    if (NU != 0)
    {
        for (label k = 0; k < NU; k++)
        {
            L_U_SUPmodes.append(Umodes[k].clone());
        }
    }
    if (NSUP != 0)
    {
        for (label k = 0; k < NSUP; k++)
        {
            L_U_SUPmodes.append(supmodes[k].clone());
        }
    }

    L_Tmodes.resize(0);
    if (liftfieldT.size() != 0)
    {
        for (label k = 0; k < liftfieldT.size(); k++)
        {
            L_Tmodes.append(liftfieldT[k].clone());
        }
    }
    if (NT != 0)
    {
        for (label k = 0; k < NT; k++)
        {
            L_Tmodes.append(Tmodes[k].clone());
        }
    }

    if (ITHACAutilities::check_folder("./ITHACAoutput/Matrices/"))
    {
        word M_str = "M_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP);
        word B_str = "B_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP);
        word BT_str = "BT_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP);
        word C_tstr = "C_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_t";
        word K_str = "K_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(NPrgh);
        word H_str = "H_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(liftfieldT.size()) + "_" + name(NT);
        word W_str = "W_" + name(liftfieldT.size()) + "_" + name(NT);
        word Y_str = "Y_" + name(liftfieldT.size()) + "_" + name(NT);
        word P_str = "P_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(NPrgh);
        word Q_str = "Q_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(liftfieldT.size()) + "_" + name(NT) + "_t";
        word CT1_str = "CT1_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(Nnut) + "_t";
        word CT2_str = "CT2_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(Nnut) + "_t";
        word CT1_ave_str = "CT1_ave_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(Nnut) + "_t";
        word CT2_ave_str = "CT2_ave_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(Nnut) + "_t";
        word YT_str = "YT_" + name(liftfield.size()) + "_" + name(NT) + "_" + name(Nnut) + "_t";
        word YT_ave_str = "YT_ave_" + name(liftfield.size()) + "_" + name(NT) + "_" + name(Nnut) + "_t";

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + M_str))
        {
            ITHACAstream::ReadDenseMatrix(M_matrix, "./ITHACAoutput/Matrices/", M_str);
        } else
        {
            M_matrix = mass_term(NU, NPrgh, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + B_str))
        {
            ITHACAstream::ReadDenseMatrix(B_matrix, "./ITHACAoutput/Matrices/", B_str);
        } else
        {
            B_matrix = diffusive_term(NU, NPrgh, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + BT_str))
        {
            ITHACAstream::ReadDenseMatrix(BT_matrix, "./ITHACAoutput/Matrices/", BT_str);
        } else
        {
            BT_matrix = BTturbulence(NU, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + C_tstr))
        {
            ITHACAstream::ReadDenseTensor(C_tensor, "./ITHACAoutput/Matrices/", C_tstr);
        } else
        {
            C_tensor = convective_term_tens(NU, NPrgh, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + K_str))
        {
            ITHACAstream::ReadDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/", K_str);
        } else
        {
            K_matrix = pressure_gradient_term(NU, NPrgh, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + H_str))
        {
            ITHACAstream::ReadDenseMatrix(H_matrix, "./ITHACAoutput/Matrices/", H_str);
        } else
        {
            H_matrix = buoyant_term(NU, NPrgh, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + W_str))
        {
            ITHACAstream::ReadDenseMatrix(W_matrix, "./ITHACAoutput/Matrices/", W_str);
        } else
        {
            W_matrix = mass_term_temperature(NU, NT, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Q_str))
        {
            ITHACAstream::ReadDenseTensor(Q_tensor, "./ITHACAoutput/Matrices/", Q_str);
        } else
        {
            Q_tensor = convective_tensor_temperature(NU, NT, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + CT1_str))
        {
            ITHACAstream::ReadDenseTensor(CT1_tensor, "./ITHACAoutput/Matrices/", CT1_str);
        } else
        {
            CT1_tensor = turbulenceTensor1(NU, NSUP, Nnut);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + CT2_str))
        {
            ITHACAstream::ReadDenseTensor(CT2_tensor, "./ITHACAoutput/Matrices/", CT2_str);
        } else
        {
            CT2_tensor = turbulenceTensor2(NU, NSUP, Nnut);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + CT1_ave_str))
        {
            ITHACAstream::ReadDenseTensor(CT1_ave_tensor, "./ITHACAoutput/Matrices/", CT1_ave_str);
        } else
        {
            CT1_ave_tensor = turbulenceAveTensor1(NU, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + CT2_ave_str))
        {
            ITHACAstream::ReadDenseTensor(CT2_ave_tensor, "./ITHACAoutput/Matrices/", CT2_ave_str);
        } else
        {
            CT2_ave_tensor = turbulenceAveTensor2(NU, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Y_str))
        {
            ITHACAstream::ReadDenseMatrix(Y_matrix, "./ITHACAoutput/Matrices/", Y_str);
        } else
        {
            Y_matrix = diffusive_term_temperature(NU, NT, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + P_str))
        {
            ITHACAstream::ReadDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/", P_str);
        } else
        {
            P_matrix = divergence_term(NU, NPrgh, NSUP);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + YT_str))
        {
            ITHACAstream::ReadDenseTensor(YT_tensor, "./ITHACAoutput/Matrices/", YT_str);
        } else
        {
            YT_tensor = temperatureTurbulenceTensor(NT, Nnut);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + YT_ave_str))
        {
            ITHACAstream::ReadDenseTensor(YT_ave_tensor, "./ITHACAoutput/Matrices/", YT_ave_str);
        } else
        {
            YT_ave_tensor = turbulenceTemperatureAveTensor(NT);
        }

    } else
    {
        M_matrix = mass_term(NU, NPrgh, NSUP);
        B_matrix = diffusive_term(NU, NPrgh, NSUP);
        K_matrix = pressure_gradient_term(NU, NPrgh, NSUP);
        H_matrix = buoyant_term(NU, NT, NSUP);
        W_matrix = mass_term_temperature(NU, NT, NSUP);
        Y_matrix = diffusive_term_temperature(NU, NT, NSUP);
        P_matrix = divergence_term(NU, NPrgh, NSUP);
        BT_matrix = BTturbulence(NU, NSUP);
        C_tensor = convective_term_tens(NU, NPrgh, NSUP);
        Q_tensor = convective_tensor_temperature(NU, NT, NSUP);
        CT1_tensor = turbulenceTensor1(NU, NSUP, Nnut);
        CT2_tensor = turbulenceTensor2(NU, NSUP, Nnut);
        CT1_ave_tensor = turbulenceAveTensor1(NU, NSUP);
        CT2_ave_tensor = turbulenceAveTensor2(NU, NSUP);
        YT_tensor = temperatureTurbulenceTensor(NT, Nnut);
        YT_ave_tensor = turbulenceTemperatureAveTensor(NT);
    }

    if (bcMethod == "penalty")
    {
        Info << "Assembling BC penalty matrices." << endl;
        bcVelVec = bcVelocityVec(NUmodes, NSUPmodes);
        bcVelMat = bcVelocityMat(NUmodes, NSUPmodes);
        bcTempVec = bcTemperatureVec(NTmodes);
        bcTempMat = bcTemperatureMat(NTmodes);
    }

    B_total_matrix = B_matrix + BT_matrix;
    label cSize = NU + NSUP + liftfield.size();
    C_total_tensor.resize(cSize, Nnut, cSize);
    C_total_tensor = CT1_tensor + CT2_tensor;
    C_total_ave_tensor.resize(cSize, avgNutfield.size(), cSize);
    C_total_ave_tensor = CT1_ave_tensor + CT2_ave_tensor;

    offlineRBFInterpolation(float(RBFShapeParameter));
}

// * * * * * * * * * * * * * * RBF Prep Methods * * * * * * * * * * * * * * //
void UnsteadyBBTurb::offlineRBFInterpolation(float RBFshapeParam) // TODO: Make this runnable in parallel
{
    samples.resize(Nnutmodes);
    rbfSplines.resize(Nnutmodes);
    Eigen::MatrixXd weights;
    Eigen::MatrixXd coeffL2nut = ITHACAutilities::getCoeffs(fluctNutfield, nutmodes, Nnutmodes);
    Eigen::MatrixXd coeffL2vel;
    coeffL2vel.resize(0, 0);
    if (bcMethod == "lift")
    {
        // TODO: Strong doubts about this. Check both Hijazi paper and PhD thesis. He/she says to use the Homogeneous field. Also says to drop the supremizer modes
        // Should the firts liftfield.size() rows be dropped from coeffL2vel later when using the RBFs? Or should they forcefully be set to the BC values?
        coeffL2vel = ITHACAutilities::getCoeffs(Uomfield, L_U_SUPmodes, NUmodes + liftfield.size()); // Returns a [modes x snapshots]
        skipRBFIndex = liftfield.size();
    } else if (bcMethod == "penalty")
    {
        coeffL2vel = ITHACAutilities::getCoeffs(Ufield, L_U_SUPmodes, NUmodes); // Returns a [modes x snapshots]
    }
    Info << "Obtained the L2 coefficients for velocity and eddy viscosity." << endl;
    Info << "Shape of the L2 velocity coeff matrix: " << coeffL2vel.rows() << " x " << coeffL2vel.cols() << endl;
    // // Substitute the BC rows with the appropriate values when using lifting approach
    // if (liftfield.size() > 0)
    // {
    //     M_Assert(coeffL2vel.rows() >= liftfield.size(), "coeffL2vel has fewer rows than liftfield.size()");
    //     // precompute total snapshots and check
    //     label totalSnapshots = 0;
    //     for (label jj = 0; jj < timeSnapshots.size(); ++jj)
    //         totalSnapshots += timeSnapshots[jj].size();
    //     M_Assert(totalSnapshots == coeffL2vel.cols(),
    //         "Total timeSnapshots columns do not match coeffL2vel.columns()");

    //     for (label i = 0; i < liftfield.size(); i++) // For each BC
    //     {
    //         label blockStart = 0;
    //         for (label j = 0; j < timeSnapshots.size(); j++)
    //         {
    //             label n = timeSnapshots[j].size();
    //             double BCValue = liftBCMatrix(j, i);
    //             // bounds check and set the correct block (columns [blockStart, blockStart+n))
    //             M_Assert(blockStart + n <= coeffL2vel.cols(),
    //                 "blockStart + n exceeds coeffL2vel columns");
    //             coeffL2vel.row(i).segment(blockStart, n).setConstant(BCValue);
    //             blockStart += n;
    //         }
    //     }
    // }

    List<Eigen::MatrixXd> velDerCoeff = velDerivativeCoeff(coeffL2vel.transpose(), coeffL2nut.transpose(), timeSnapshots);
    dimA = velDerCoeff[0].cols();
    meanA.resize(velDerCoeff[0].cols());
    stdA.resize(velDerCoeff[0].cols());
    // Normalize and scale the data for better RBF performance. We save the mean and std for later use during the online phase.
    for (label i = 0; i < velDerCoeff[0].cols(); i++)
    {
        meanA(i) = velDerCoeff[0].col(i).mean();
        stdA(i) = std::sqrt((velDerCoeff[0].col(i).array() - meanA(i)).square().sum() / (velDerCoeff[0].rows() - 1));
        if (stdA(i) == 0)
        {
            stdA(i) = 1; // To avoid division by zero in case of constant mode
        }
        // velDerCoeff[0].col(i) = (velDerCoeff[0].col(i).array() - meanA(i)) / stdA(i);
    }
    meanG.resize(velDerCoeff[1].cols());
    stdG.resize(velDerCoeff[1].cols());
    for (label i = 0; i < velDerCoeff[1].cols(); i++)
    {
        meanG(i) = velDerCoeff[1].col(i).mean();
        stdG(i) = std::sqrt((velDerCoeff[1].col(i).array() - meanG(i)).square().sum() / (velDerCoeff[1].rows() - 1));
        if (stdG(i) == 0)
        {
            stdG(i) = 1; // To avoid division by zero in case of constant mode
        }
        // velDerCoeff[1].col(i) = (velDerCoeff[1].col(i).array() - meanG(i)) / stdG(i);
    }

    if (Pstream::master())
    {
        ITHACAutilities::createSymLink("./ITHACAoutput/Debug");
        ITHACAstream::exportMatrix(velDerCoeff[0], "A_RBF", "python", "./ITHACAoutput/Debug/");
        ITHACAstream::exportMatrix(velDerCoeff[1], "G_RBF", "python", "./ITHACAoutput/Debug/");
    }

    // Guesstimation of the shapeParam of the RBFs based on the average distance between points
    if (ITHACAutilities::check_file("./ITHACAoutput/shapeParameters"))
    {
        shapeParameters = ITHACAstream::readMatrix("./ITHACAoutput/shapeParameters");
        M_Assert(shapeParameters.size() == Nnutmodes,
            "Thes size of the shape parameters vector must be equal to the number of eddy viscosity modes nNutModes");
        Info << "RBF shape parameters (shapeParameters) read from file as: " << shapeParameters << endl;
    } else
    {
        shapeParameters = Eigen::MatrixXd::Ones(Nnutmodes, 1) * RBFshapeParam;
    }

    if (estimateShapeParam)
    {
        Info << "Estimating RBF shape parameters using LOOCV..." << endl;

        // Define search range around the provided guess
        double minEps = 0.1 * RBFshapeParam;
        double maxEps = 10.0 * RBFshapeParam;
        int nSteps = 10;

        for (label i = 0; i < Nnutmodes; i++)
        {
            double bestEps = RBFshapeParam;
            double minErr = 1.0e30;

            for (int s = 0; s < nSteps; ++s)
            {
                double eps = minEps + (maxEps - minEps) * double(s) / double(nSteps - 1);
                double errSum = 0.0;
                label nSamples = velDerCoeff[0].rows();

                for (label k = 0; k < nSamples; ++k)
                {
                    SPLINTER::DataTable looSamples(true, true);
                    for (label j = 0; j < nSamples; ++j)
                    {
                        if (j == k)
                            continue;
                        looSamples.addSample(Eigen::VectorXd(velDerCoeff[0].row(j)), velDerCoeff[1](j, i));
                    }
                    SPLINTER::RBFSpline spline(looSamples, SPLINTER::RadialBasisFunctionType::GAUSSIAN, false, eps);

                    double pred = spline.eval(Eigen::VectorXd(velDerCoeff[0].row(k)));
                    double act = velDerCoeff[1](k, i);
                    errSum += std::pow(pred - act, 2);
                }

                if (errSum < minErr)
                {
                    minErr = errSum;
                    bestEps = eps;
                }
            }
            shapeParameters(i) = bestEps;
            Info << "Mode " << i << " optimized shape param: " << bestEps << " (MSE: " << minErr / velDerCoeff[0].rows() << ")" << endl;
        }

        Info << "RBF shape parameters (shapeParameters) estimated as: " << shapeParameters << endl;
        if (Pstream::master())
        {
            ITHACAstream::SaveDenseMatrix(shapeParameters, "./ITHACAoutput/", "RBFshapeParameters");
        }
    }

    // Interpolation of the eddy viscosity proceeds per-mode
    for (label i = 0; i < Nnutmodes; i++)
    {
        word weightName = "wRBF_N" + name(i + 1) + "_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes);
        if (ITHACAutilities::check_file("./ITHACAoutput/weightsRBF/" + weightName))
        {
            samples[i] = new SPLINTER::DataTable(true, true);
            for (label j = 0; j < velDerCoeff[1].rows(); j++)
            {
                samples[i]->addSample(velDerCoeff[0].row(j), velDerCoeff[1](j, i));
            }
            ITHACAstream::ReadDenseMatrix(weights, "./ITHACAoutput/weightsRBF/", weightName);
            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                SPLINTER::RadialBasisFunctionType::GAUSSIAN, weights, shapeParameters(i));
        } else
        {
            samples[i] = new SPLINTER::DataTable(true, true);
            for (label j = 0; j < velDerCoeff[1].rows(); j++)
            {
                samples[i]->addSample(velDerCoeff[0].row(j), velDerCoeff[1](j, i));
            }

            rbfSplines[i] = new SPLINTER::RBFSpline(*samples[i],
                SPLINTER::RadialBasisFunctionType::GAUSSIAN, false, shapeParameters(i));
            if (Pstream::master())
            {
                ITHACAstream::SaveDenseMatrix(rbfSplines[i]->weights,
                    "./ITHACAoutput/weightsRBF/", weightName);
                ITHACAstream::exportMatrix(rbfSplines[i]->weights,
                    "weightsRBF", "python", "./ITHACAoutput/weightsRBF/");
            }
        }
    }
}

List<Eigen::MatrixXd> UnsteadyBBTurb::velDerivativeCoeff(const Eigen::MatrixXd& A,
    const Eigen::MatrixXd& G, const List<Eigen::VectorXd>& snapshotTimes)
{
    List<Eigen::MatrixXd> newCoeffs;
    newCoeffs.setSize(2);
    const label velCoeffsNum = A.cols();
    const label snapshotsNum = A.rows();
    const label parsSamplesNum = snapshotTimes.size();

    Eigen::VectorXi timeSnapshotsPerSampleVec(parsSamplesNum);
    for (label i = 0; i < parsSamplesNum; i++)
    {
        timeSnapshotsPerSampleVec(i) = snapshotTimes[i].size();
    }
    const label newColsNum = 2 * velCoeffsNum;
    const label newRowsNum = timeSnapshotsPerSampleVec.sum() - 1 * parsSamplesNum;
    newCoeffs[0].resize(newRowsNum, newColsNum);
    newCoeffs[1].resize(newRowsNum, G.cols());
    label outOffset = 0;
    for (label j = 0; j < parsSamplesNum; j++)
    {
        const Eigen::VectorXd& timeSnap = snapshotTimes[j];
        // Create shifted blocks to compute differences. Remember that
        // A has all the snapshots stacked for all parameter samples, with the column
        // indicating the different time-varying coefficients (a1, a2, ..., an)
        label rowsPerBlock = timeSnapshotsPerSampleVec(j) - 1;
        label blockStart = timeSnapshotsPerSampleVec.head(j).sum();
        Eigen::MatrixXd b0 = A.middleRows(blockStart, rowsPerBlock);
        Eigen::MatrixXd b2 = A.middleRows(blockStart + 1, rowsPerBlock);
        Eigen::VectorXd deltaT = timeSnap.tail(rowsPerBlock) - timeSnap.head(rowsPerBlock);
        Eigen::MatrixXd derivative = (b2 - b0).array().colwise() / deltaT.array();
        Eigen::MatrixXd bNew(rowsPerBlock, newColsNum);
        bNew.leftCols(velCoeffsNum) = b2;
        bNew.rightCols(velCoeffsNum) = derivative;
        newCoeffs[0].block(outOffset, 0, rowsPerBlock, newColsNum) = bNew;
        newCoeffs[1].middleRows(outOffset, rowsPerBlock) =
            G.middleRows(blockStart + 1, rowsPerBlock);
        outOffset += rowsPerBlock;
    }
    return newCoeffs;
}

void UnsteadyBBTurb::splitEddyViscositySnapshots() // TODO: (Maybe) Use the averaging/subtractions function from ITHACAutilities
{
    label nSamples = timeSnapshots.size();
    label globalIndex = 0;
    avgNutfield.clear();
    fluctNutfield.clear();

    for (label i = 0; i < nSamples; i++)
    {
        label nSnapshotPerSample = timeSnapshots[i].size();
        Info << "Processing sample " << i + 1 << " with " << nSnapshotPerSample << " snapshots." << endl;
        M_Assert(nSnapshotPerSample > 0,
            "Each parameter sample must have at least one snapshot");
        autoPtr<volScalarField> avgPtr(
            new volScalarField(
                IOobject(
                    "avgNut",
                    Nutfield[globalIndex].time().timeName(),
                    Nutfield[globalIndex].mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                Nutfield[globalIndex]));

        for (label j = 1; j < nSnapshotPerSample; j++)
        {
            avgPtr() += Nutfield[globalIndex + j];
        }

        avgPtr() /= scalar(nSnapshotPerSample);
        avgNutfield.append(avgPtr);
        globalIndex += nSnapshotPerSample;
    }

    label totalSnapshots = Nutfield.size();
    fluctNutfield.setSize(totalSnapshots);

    globalIndex = 0;
    label flatIndex = 0;

    for (label i = 0; i < nSamples; i++)
    {
        label nSnap = timeSnapshots[i].size();
        const volScalarField& avgField = avgNutfield[i];

        for (label j = 0; j < nSnap; j++)
        {
            tmp<volScalarField> tempFluct = Nutfield[globalIndex + j] - avgField;
            tempFluct.ref().rename("fluctNut");
            fluctNutfield.set(flatIndex, tempFluct.ptr());
            flatIndex++;
        }
        globalIndex += nSnap;
    }

    if (DEBUG_MODE)
    {
        ITHACAutilities::createSymLink("./ITHACAoutput/Debug");
        ITHACAstream::exportFields(avgNutfield, "./ITHACAoutput/Debug/", "avgNutfield");
        ITHACAstream::exportFields(fluctNutfield, "./ITHACAoutput/Debug/", "fluctNutfield");
    }
}

// * * * * * * * * * * * * * * Matrices Methods * * * * * * * * * * * * * * //

Eigen::MatrixXd UnsteadyBBTurb::pressure_gradient_term(label NUmodes,
    label NPrghmodes, label NSUPmodes)
{
    label K1size = NUmodes + NSUPmodes + liftfield.size();
    label K2size = NPrghmodes;
    Eigen::MatrixXd K_matrix(K1size, K2size);
    dimensionedVector g = _g();

    // Project everything
    for (label i = 0; i < K1size; i++)
    {
        for (label j = 0; j < K2size; j++)
        {
            K_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] &
                fvc::reconstruct(fvc::snGrad(P_rghmodes[j]) *
                    P_rghmodes[j].mesh().magSf()))
                                 .value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(K_matrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        // Export the matrix
        ITHACAstream::SaveDenseMatrix(K_matrix, "./ITHACAoutput/Matrices/",
            "K_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPrghmodes));
    }

    return K_matrix;
}

Eigen::MatrixXd UnsteadyBBTurb::diffusive_term_temperature(label NUmodes,
    label NTmodes, label NSUPmodes)
{
    label Ysize = NTmodes + liftfieldT.size();
    Eigen::MatrixXd Y_matrix(Ysize, Ysize);

    for (label i = 0; i < Ysize; i++)
    {
        for (label j = 0; j < Ysize; j++)
        {
            Y_matrix(i, j) = fvc::domainIntegrate(L_Tmodes[i] * fvc::laplacian(dimensionedScalar("1", dimless, 1), L_Tmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(Y_matrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        // Export the matrix
        ITHACAstream::SaveDenseMatrix(Y_matrix, "./ITHACAoutput/Matrices/",
            "Y_" + name(liftfieldT.size()) + "_" + name(NTmodes));
    }

    return Y_matrix;
}

Eigen::MatrixXd UnsteadyBBTurb::divergence_term(label NUmodes, label NPrghmodes,
    label NSUPmodes)
{
    label P1size = NPrghmodes;
    label P2size = NUmodes + NSUPmodes + liftfield.size();
    Eigen::MatrixXd P_matrix(P1size, P2size);

    // Project everything
    for (label i = 0; i < P1size; i++)
    {
        for (label j = 0; j < P2size; j++)
        {
            P_matrix(i, j) = fvc::domainIntegrate(P_rghmodes[i] * fvc::div(L_U_SUPmodes[j])).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(P_matrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        // Export the matrix
        ITHACAstream::SaveDenseMatrix(P_matrix, "./ITHACAoutput/Matrices/",
            "P_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPrghmodes));
    }

    return P_matrix;
}

Eigen::MatrixXd UnsteadyBBTurb::buoyant_term(label NUmodes, label NTmodes,
    label NSUPmodes)
{
    label H1size = NUmodes + liftfield.size() + NSUPmodes;
    label H2size = NTmodes + liftfieldT.size();
    Eigen::MatrixXd H_matrix(H1size, H2size);
    dimensionedScalar beta = _beta();
    dimensionedScalar TRef = _TRef();
    dimensionedVector g = _g();
    // volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();

    // Project everything
    for (label i = 0; i < H1size; i++)
    {
        for (label j = 0; j < H2size; j++)
        {
            H_matrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::reconstruct(ghf * fvc::snGrad(1.0 - (beta * (L_Tmodes[j] - TRef))) * L_Tmodes[j].mesh().magSf())).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(H_matrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        // Export the matrix
        ITHACAstream::SaveDenseMatrix(H_matrix, "./ITHACAoutput/Matrices/",
            "H_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(liftfieldT.size()) + "_" + name(NTmodes));
    }

    return H_matrix;
}

Eigen::MatrixXd UnsteadyBBTurb::mass_term_temperature(label NUmodes, label NTmodes,
    label NSUPmodes)
{
    label Wsize = NTmodes + liftfieldT.size();
    Eigen::MatrixXd W_matrix(Wsize, Wsize);

    for (label i = 0; i < Wsize; i++)
    {
        for (label j = 0; j < Wsize; j++)
        {
            W_matrix(i, j) = fvc::domainIntegrate(L_Tmodes[i] * L_Tmodes[j]).value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(W_matrix, sumOp<Eigen::MatrixXd>());
    }
    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(W_matrix, "./ITHACAoutput/Matrices/",
            "W_" + name(liftfieldT.size()) + "_" + name(NTmodes));
    }

    return W_matrix;
}

Eigen::MatrixXd UnsteadyBBTurb::BTturbulence(label NU, label NSUP)
{
    label btSize = NU + NSUP + liftfield.size();
    Eigen::MatrixXd btMatrix(btSize, btSize);
    btMatrix = btMatrix * 0;

    for (label i = 0; i < btSize; i++)
    {
        for (label j = 0; j < btSize; j++)
        {
            btMatrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] & (fvc::div(dev((T(fvc::grad(L_U_SUPmodes[j]))))))).value();
        }
    }
    if (Pstream::parRun())
    {
        reduce(btMatrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(btMatrix, "./ITHACAoutput/Matrices/",
            "BT_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes));
    }
    return btMatrix;
}

Eigen::Tensor<double, 3> UnsteadyBBTurb::temperatureTurbulenceTensor(label NT, label Nnut)
{
    label Stsize = NT + liftfieldT.size();
    Eigen::Tensor<double, 3> YT_tensor;
    YT_tensor.resize(Stsize, Nnut, Stsize);
    for (label i = 0; i < Stsize; i++)
    {
        for (label j = 0; j < Nnut; j++)
        {
            for (label k = 0; k < Stsize; k++)
            {
                YT_tensor(i, j, k) = fvc::domainIntegrate(L_Tmodes[i] * fvc::laplacian(nutmodes[j], L_Tmodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(YT_tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(YT_tensor, "./ITHACAoutput/Matrices/",
            "YT_" + name(liftfield.size()) + "_" + name(NT) + "_" + name(Nnut) + "_t");
    }
    return YT_tensor;
}

Eigen::Tensor<double, 3> UnsteadyBBTurb::turbulenceTensor1(label NU, label NSUP, label Nnut)
{
    label cSize = NU + NSUP + liftfield.size();
    Eigen::Tensor<double, 3> ct1Tensor;
    ct1Tensor.resize(cSize, Nnut, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < Nnut; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(nutmodes[j], L_U_SUPmodes[k])).value();
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(ct1Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ct1Tensor, "./ITHACAoutput/Matrices/",
            "CT1_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(Nnut) + "_t");
    }
    return ct1Tensor;
}

Eigen::Tensor<double, 3> UnsteadyBBTurb::turbulenceAveTensor1(label NU, label NSUP)
{
    label cSize = NU + NSUP + liftfield.size();
    Eigen::Tensor<double, 3> ct1AveTensor;
    label samplesNumber = avgNutfield.size();
    ct1AveTensor.resize(cSize, samplesNumber, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct1AveTensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] & fvc::laplacian(avgNutfield[j], L_U_SUPmodes[k])).value();
            }
        }
    }
    if (Pstream::parRun())
    {
        reduce(ct1AveTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ct1AveTensor, "./ITHACAoutput/Matrices/",
            "CT1Ave_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_t");
    }
    return ct1AveTensor;
}

Eigen::Tensor<double, 3> UnsteadyBBTurb::turbulenceTensor2(label NU, label NSUP, label Nnut)
{
    label cSize = NU + NSUP + liftfield.size();
    Eigen::Tensor<double, 3> ct2Tensor;
    ct2Tensor.resize(cSize, Nnut, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < Nnut; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2Tensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] &
                    (fvc::div(nutmodes[j] * dev((fvc::grad(L_U_SUPmodes[k]))().T()))))
                                         .value();
            }
        }
    }
    if (Pstream::parRun())
    {
        reduce(ct2Tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ct2Tensor, "./ITHACAoutput/Matrices/",
            "CT2_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_" + name(Nnut) + "_t");
    }

    return ct2Tensor;
}

Eigen::Tensor<double, 3> UnsteadyBBTurb::turbulenceAveTensor2(label NU, label NSUP)
{
    label cSize = NU + NSUP + liftfield.size();
    Eigen::Tensor<double, 3> ct2AveTensor;
    label samplesNumber = avgNutfield.size();
    ct2AveTensor.resize(cSize, samplesNumber, cSize);

    for (label i = 0; i < cSize; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < cSize; k++)
            {
                ct2AveTensor(i, j, k) = fvc::domainIntegrate(L_U_SUPmodes[i] &
                    (fvc::div(avgNutfield[j] * dev((fvc::grad(L_U_SUPmodes[k]))().T()))))
                                            .value();
            }
        }
    }
    if (Pstream::parRun())
    {
        reduce(ct2AveTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(ct2AveTensor, "./ITHACAoutput/Matrices/",
            "CT2Ave_" + name(liftfield.size()) + "_" + name(NU) + "_" + name(NSUP) + "_t");
    }
    return ct2AveTensor;
}

Eigen::Tensor<double, 3> UnsteadyBBTurb::turbulenceTemperatureAveTensor(label NT)
{
    label ySize = NT + liftfieldT.size();
    Eigen::Tensor<double, 3> YTAveTensor;
    label samplesNumber = avgNutfield.size();
    YTAveTensor.resize(ySize, samplesNumber, ySize);
    for (label i = 0; i < ySize; i++)
    {
        for (label j = 0; j < samplesNumber; j++)
        {
            for (label k = 0; k < ySize; k++)
            {
                YTAveTensor(i, j, k) = fvc::domainIntegrate(L_Tmodes[i] * fvc::laplacian(avgNutfield[j], L_Tmodes[k])).value();
            }
        }
    }
    if (Pstream::parRun())
    {
        reduce(YTAveTensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(YTAveTensor, "./ITHACAoutput/Matrices/",
            "YT_ave_" + name(liftfieldT.size()) + "_" + name(NT) + "_t");
    }
    return YTAveTensor;
}

// * * * * * * * * * * * * * * Energy Eq. Methods * * * * * * * * * * * * * //
Eigen::Tensor<double, 3> UnsteadyBBTurb::convective_tensor_temperature(label NU,
    label NT, label NSUP)
{
    label Qsize = NU + liftfield.size() + NSUP;
    label Qsizet = NT + liftfieldT.size();
    Eigen::Tensor<double, 3> Q_tensor;
    Q_tensor.resize(Qsizet, Qsize, Qsizet);

    for (label i = 0; i < Qsizet; i++)
    {
        for (label j = 0; j < Qsize; j++)
        {
            for (label k = 0; k < Qsizet; k++)
            {
                Q_tensor(i, j, k) = fvc::domainIntegrate(L_Tmodes[i] * fvc::div(fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(), L_Tmodes[k])).value();
            }
        }
    }
    if (Pstream::parRun())
    {
        reduce(Q_tensor, sumOp<Eigen::Tensor<double, 3>>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseTensor(Q_tensor, "./ITHACAoutput/Matrices/",
            "Q_" + name(liftfield.size()) + "_" + name(NU) + "_" +
                name(NSUP) + "_" + name(liftfieldT.size()) + "_" + name(NT) + "_t");
    }
    return Q_tensor;
}

// * * * * * * * * * * * * * * Penalty term methods * * * * * * * * * * * * * //

List<Eigen::MatrixXd> UnsteadyBBTurb::bcTemperatureVec(label NTmodes)
{
    label BCsize = NTmodes;
    List<Eigen::MatrixXd> bcTempVec(inletIndexT.rows());

    for (label j = 0; j < inletIndexT.rows(); j++)
    {
        bcTempVec[j].resize(BCsize, 1);
    }

    for (label i = 0; i < inletIndexT.rows(); i++)
    {
        label BCind = inletIndexT(i);
        for (label j = 0; j < BCsize; j++)
        {
            bcTempVec[i](j, 0) = gSum(L_Tmodes[j].boundaryField()[BCind] * L_Tmodes[j].mesh().magSf().boundaryField()[BCind]);
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::exportMatrix(bcTempVec, "bcTempVec", "python", "./ITHACAoutput/Matrices/bcTempVec/");
    }

    return bcTempVec;
}

List<Eigen::MatrixXd> UnsteadyBBTurb::bcTemperatureMat(label NTmodes)
// Compute the L2 inner product matrix of the temperature BCs on the inlet patches
{
    label BCsize = NTmodes;
    label BCTsize = inletIndexT.rows();
    List<Eigen::MatrixXd> bcTempMat(BCTsize);

    for (label j = 0; j < BCTsize; j++)
    {
        bcTempMat[j].resize(BCsize, BCsize);
    }

    for (label k = 0; k < BCTsize; k++)
    {
        label BCind = inletIndexT(k);

        for (label i = 0; i < BCsize; i++)
        {
            for (label j = 0; j < BCsize; j++)
            {
                // Corrected to use bcTempMat and L_Tmodes (scalar) instead of velocity modes
                bcTempMat[k](i, j) = gSum(L_Tmodes[i].boundaryField()[BCind] *
                    L_Tmodes[j].boundaryField()[BCind] *
                    L_Tmodes[i].mesh().magSf().boundaryField()[BCind]);
            }
        }
    }

    return bcTempMat;
}

// * * * * * * * * * * * * * * Restart Methods * * * * * * * * * * * * * //

void UnsteadyBBTurb::restart()
{
    Time& runTime = _runTime();
    fvMesh& mesh = _mesh();
    runTime.setTime(0, 0);

    // Read transportProperties dictionary
    IOdictionary transportProperties(
        IOobject(
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE));

    // Update constants
    _nu() = dimensionedScalar("nu", dimViscosity, transportProperties.lookup("nu"));
    _Pr() = dimensionedScalar("Pr", dimless, transportProperties.lookup("Pr"));
    _Prt() = dimensionedScalar("Prt", dimless, transportProperties.lookup("Prt"));
    _beta() = dimensionedScalar("beta", dimless / dimTemperature, transportProperties.lookup("beta"));
    _TRef() = dimensionedScalar("TRef", dimTemperature, transportProperties.lookup("TRef"));

    // Now we reset the fields, reading from disk (OpenFOAM folders)
    volVectorField U_new(IOobject("U", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _U() = U_new;
    // _U().correctBoundaryConditions();

    volScalarField T_new(IOobject("T", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _T() = T_new;
    // _T().correctBoundaryConditions();

    volScalarField p_rgh_new(IOobject("p_rgh", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _p_rgh() = p_rgh_new;
    // _p_rgh().correctBoundaryConditions();

    volScalarField nut_new(IOobject("nut", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _nut() = nut_new;
    // _nut().correctBoundaryConditions();

    volScalarField alphat_new(IOobject("alphat", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _alphat() = alphat_new;
    // _alphat().correctBoundaryConditions();

    volScalarField UliftBC_new(IOobject("UliftBC", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::NO_WRITE), mesh);
    _UliftBC() = UliftBC_new;

    _phi() = linearInterpolate(_U()) & mesh.Sf();
    _rhok() = 1.0 - _beta() * (_T() - _TRef());

    _laminarTransport.clear();
    _laminarTransport = autoPtr<singlePhaseTransportModel>(new singlePhaseTransportModel(_U(), _phi()));
    turbulence.clear();
    turbulence = autoPtr<incompressible::turbulenceModel>(incompressible::turbulenceModel::New(_U(), _phi(), _laminarTransport()));
    turbulence->validate();

    _p() = _p_rgh() + _rhok() * _gh();

    // Fix p reference
    pimpleControl& pimple = _pimple();
    setRefCell(_p(), _p_rgh(), pimple.dict(), pRefCell, pRefValue);
    if (_p_rgh().needReference())
    {
        _p() += dimensionedScalar("p", _p().dimensions(), pRefValue - getRefCellValue(_p(), pRefCell));
    }
    mesh.setFluxRequired(_p_rgh().name());
    Info << "Restart complete." << endl;
}
