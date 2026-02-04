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

    M_Assert(method == "supremizer" || method == "PPE" || method == "VMB",
        "A method must be set to either 'supremizer', 'VMB' or 'PPE' in ITHACAdict");
    M_Assert(bcMethod == "lift" || bcMethod == "penalty",
        "The BC method must be set to lift or penalty in ITHACAdict");
    M_Assert(timedepbcMethod == "yes" || timedepbcMethod == "no",
        "The BC method can be set to yes or no");
    turbulence->validate();
    offline = ITHACAutilities::check_off();
    podex = ITHACAutilities::check_pod();
    supex = ITHACAutilities::check_sup();
    viscDict = ITHACAdict->subDict("viscDict");
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
            nut = turbulence->nut();
            ITHACAstream::exportSolution(U, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(p_rgh, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(T, name(counter), "./ITHACAoutput/Offline/");
            ITHACAstream::exportSolution(nut, name(counter), "./ITHACAoutput/Offline/");
            std::ofstream of("./ITHACAoutput/Offline/" + name(counter) + "/" + runTime.timeName());

            timeSnapshots[nSample](stepCounter) = runTime.value();
            stepCounter++;

            Ufield.append(U.clone());
            Prghfield.append(p_rgh.clone());
            Tfield.append(T.clone());
            Nutfield.append(nut.clone());
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
    Info << "Called the supremizer for UnsteadyBBTurb." << nl << endl;
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
    label NSUP, label Nnut)
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

    offlineRBFInterpolation();
}

void UnsteadyBBTurb::projectPPE(fileName folder, label NU, label NPrgh, label NT,
    label Nnut)
{
    NUmodes = NU;
    NTmodes = NT;
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

    M_matrix = mass_term(NU, NPrgh, 0);
    B_matrix = diffusive_term(NU, NPrgh, 0);
    K_matrix = pressure_gradient_term(NU, NPrgh, 0);
    H_matrix = buoyant_term(NU, NT, 0);
    W_matrix = mass_term_temperature(NU, NT, 0);
    Y_matrix = diffusive_term_temperature(NU, NT, 0);
    P_matrix = divergence_term(NU, NPrgh, 0);
    BT_matrix = BTturbulence(NU, 0);
    C_tensor = convective_term_tens(NU, NPrgh, 0);
    Q_tensor = convective_tensor_temperature(NU, NT, 0);
    D_matrix = laplacian_pressure(NPrghmodes);
    G_tensor = divMomentum(NUmodes, NPrghmodes);
    HP_matrix = buoyant_term_poisson(NPrghmodes, NTmodes);
    BC1_matrix = pressure_BC1(NUmodes, NPrghmodes);
    BC2_tensor = pressure_BC2(NUmodes, NPrghmodes);
    BC3_matrix = pressure_BC3(NUmodes, NPrghmodes);
    CT1_tensor = turbulenceTensor1(NU, 0, Nnut);
    CT2_tensor = turbulenceTensor2(NU, 0, Nnut);
    CT1_ave_tensor = turbulenceAveTensor1(NU, 0);
    CT2_ave_tensor = turbulenceAveTensor2(NU, 0);
    YT_tensor = temperatureTurbulenceTensor(NT, Nnut);
    YT_ave_tensor = turbulenceTemperatureAveTensor(NT);

    if (bcMethod == "penalty")
    {
        Info << "Assembling BC penalty matrices." << endl;
        bcVelVec = bcVelocityVec(NUmodes, 0);
        bcVelMat = bcVelocityMat(NUmodes, 0);
        bcTempVec = bcTemperatureVec(NTmodes);
        bcTempMat = bcTemperatureMat(NTmodes);
    }

    B_total_matrix = B_matrix + BT_matrix;
    label cSize = NU + liftfield.size();
    C_total_tensor.resize(cSize, Nnut, cSize);
    C_total_tensor = CT1_tensor + CT2_tensor;
    C_total_ave_tensor.resize(cSize, avgNutfield.size(), cSize);
    C_total_ave_tensor = CT1_ave_tensor + CT2_ave_tensor;
}

void UnsteadyBBTurb::projectVMB(fileName folder, label NU, label NT, label Nnut)
{
    NUmodes = NU;
    NTmodes = NT;
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
        word W_str = "W_" + name(liftfieldT.size()) + "_" + name(NT);
        word Y_str = "Y_" + name(liftfieldT.size()) + "_" + name(NT);

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + W_str))
        {
            ITHACAstream::ReadDenseMatrix(W_matrix, "./ITHACAoutput/Matrices/", W_str);
        } else
        {
            W_matrix = mass_term_temperature(NU, NT, 0);
        }

        if (ITHACAutilities::check_file("./ITHACAoutput/Matrices/" + Y_str))
        {
            ITHACAstream::ReadDenseMatrix(Y_matrix, "./ITHACAoutput/Matrices/", Y_str);
        } else
        {
            Y_matrix = diffusive_term_temperature(NU, NT, 0);
        }

    } else
    {
        M_matrix = mass_term(NU, NU, 0);
        B_matrix = diffusive_term(NU, NU, 0);
        K_matrix = pressure_gradient_term(NU, NU, 0);
        H_matrix = buoyant_term(NU, NT, 0);
        W_matrix = mass_term_temperature(NU, NT, 0);
        Y_matrix = diffusive_term_temperature(NU, NT, 0);
        BT_matrix = BTturbulence(NU, 0);
        C_tensor = convective_term_tens(NU, NU, 0);
        Q_tensor = convective_tensor_temperature(NU, NT, 0);
        CT1_tensor = turbulenceTensor1(NU, 0, Nnut);
        CT2_tensor = turbulenceTensor2(NU, 0, Nnut);
        CT1_ave_tensor = turbulenceAveTensor1(NU, 0);
        CT2_ave_tensor = turbulenceAveTensor2(NU, 0);
        YT_tensor = temperatureTurbulenceTensor(NT, Nnut);
        YT_ave_tensor = turbulenceTemperatureAveTensor(NT);
    }
    if (bcMethod == "penalty")
    {
        Info << "Assembling BC penalty matrices." << endl;
        bcVelVec = bcVelocityVec(NUmodes, 0);
        bcVelMat = bcVelocityMat(NUmodes, 0);
        bcTempVec = bcTemperatureVec(NTmodes);
        bcTempMat = bcTemperatureMat(NTmodes);
    }

    B_total_matrix = B_matrix + BT_matrix;
    label cSize = NU + liftfield.size();
    C_total_tensor.resize(cSize, Nnut, cSize);
    C_total_tensor = CT1_tensor + CT2_tensor;
    C_total_ave_tensor.resize(cSize, avgNutfield.size(), cSize);
    C_total_ave_tensor = CT1_ave_tensor + CT2_ave_tensor;

    offlineRBFInterpolation();
}

// * * * * * * * * * * * * * * RBF Prep Methods * * * * * * * * * * * * * * //
void UnsteadyBBTurb::offlineRBFInterpolation() // TODO: Make this runnable in parallel
{
    Eigen::MatrixXd weights;
    Eigen::MatrixXd coeffL2nut = ITHACAutilities::getCoeffs(fluctNutfield, nutmodes, Nnutmodes);
    Eigen::MatrixXd coeffL2vel;
    coeffL2vel.resize(0, 0);
    if (bcMethod == "lift")
    {
        coeffL2vel = ITHACAutilities::getCoeffs(Uomfield, Umodes, NUmodes); // Returns a [modes x snapshots]
        skipRBFIndex = liftfield.size();
    } else if (bcMethod == "penalty")
    {
        coeffL2vel = ITHACAutilities::getCoeffs(Ufield, Umodes, NUmodes); // Returns a [modes x snapshots]
    }
    Info << "Shape of the L2 velocity coeff matrix: " << coeffL2vel.rows() << " x " << coeffL2vel.cols() << endl;
    Info << "Shape of the L2 eddy viscosity coeff matrix: " << coeffL2nut.rows() << " x " << coeffL2nut.cols() << endl;

    // Returns a list of two matrices: [0] = velocity derivative coeffs, [1] = eddy viscosity coeffs. Each
    // matrix is [snapshots x coeffs]
    List<Eigen::MatrixXd> velDerCoeff = velDerivativeCoeff(coeffL2vel.transpose(), coeffL2nut.transpose(), timeSnapshots);
    dimA = velDerCoeff[0].cols();

    if (Pstream::master())
    {
        ITHACAutilities::createSymLink("./ITHACAoutput/Debug");
        ITHACAstream::exportMatrix(velDerCoeff[0], "A_RBF", "python", "./ITHACAoutput/Debug/");
        ITHACAstream::exportMatrix(velDerCoeff[1], "G_RBF", "python", "./ITHACAoutput/Debug/");
    }

    if (ITHACAutilities::check_file("./ITHACAoutput/shapeParameters"))
    {
        Info << "### MESSAGE - Reading RBF shape parameters from file is not yet implemented for mathtoolbox usage. Using the one provided in the dictionary." << endl;
    }

    rbfSplines.resize(Nnutmodes);
    for (label i = 0; i < Nnutmodes; i++)
    {
        // Create a RBF interpolator instance
        rbfSplines[i] = std::make_shared<ithacaInterpolator>(viscDict);

        Eigen::MatrixXd x = velDerCoeff[0].transpose();
        Eigen::VectorXd y = velDerCoeff[1].col(i);
        Info << "### DEBUG: Shape of x: " << x.rows() << " x " << x.cols() << endl;
        Info << "### DEBUG: Shape of y: " << y.rows() << " x " << y.cols() << endl;

        rbfSplines[i]->fit(x, y);

        Info << "Fitting ithacaInterpolator for mode " << i + 1 << " completed." << endl;
    }

    Info << "RBF Interpolation completed. Some details of the interpolation are as follows." << endl;
    rbfSplines[0]->printInfo();
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
        M_Assert(nSnapshotPerSample > 0, "Each parameter sample must have at least one snapshot");

        // PtrList<volScalarField> sampleNutFieldPtr;
        // sampleNutFieldPtr.setSize(nSnapshotPerSample);
        // label startIndex = globalIndex;
        // label endIndex = globalIndex + nSnapshotPerSample - 1;
        // for (label j = startIndex; j <= endIndex; j++)
        // {
        //     sampleNutFieldPtr.set(j - startIndex, Nutfield[j]);
        // }
        // ITHACAutilities::computeAverage(sampleNutFieldPtr);

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
    // dimensionedVector g = _g();

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
        ITHACAstream::exportMatrix(K_matrix, "K_matrix", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportMatrix(Y_matrix, "Y_matrix", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportMatrix(P_matrix, "P_matrix", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportMatrix(H_matrix, "H_matrix", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportMatrix(W_matrix, "W_matrix", "python", "./ITHACAoutput/Matrices/python/");
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
            btMatrix(i, j) = fvc::domainIntegrate(L_U_SUPmodes[i] &
                (fvc::div(dev2((T(fvc::grad(L_U_SUPmodes[j])))))))
                                 .value();
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
        ITHACAstream::exportMatrix(btMatrix, "BT_matrix", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportTensor(YT_tensor, "YT_tensor", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportTensor(ct1Tensor, "CT1_tensor", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportTensor(ct1AveTensor, "CT1_ave_tensor", "python", "./ITHACAoutput/Matrices/python/");
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
                    (fvc::div(nutmodes[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T()))))
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
        ITHACAstream::exportTensor(ct2Tensor, "CT2_tensor", "python", "./ITHACAoutput/Matrices/python/");
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
                    (fvc::div(avgNutfield[j] * dev2((fvc::grad(L_U_SUPmodes[k]))().T()))))
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
        ITHACAstream::exportTensor(ct2AveTensor, "CT2_ave_tensor", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportTensor(YTAveTensor, "YT_ave_tensor", "python", "./ITHACAoutput/Matrices/python/");
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
        ITHACAstream::exportTensor(Q_tensor, "Q_tensor", "python", "./ITHACAoutput/Matrices/python/");
    }
    return Q_tensor;
}

Eigen::MatrixXd UnsteadyBBTurb::laplacian_pressure(label NPrghmodes)
{
    label Dsize = NPrghmodes + liftfieldP.size();
    Eigen::MatrixXd D_matrix(Dsize, Dsize);

    // Project everything
    for (label i = 0; i < Dsize; i++)
    {
        for (label j = 0; j < Dsize; j++)
        {
            D_matrix(i, j) = fvc::domainIntegrate(fvc::grad(P_rghmodes[i]) & fvc::grad(P_rghmodes[j])).value();
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(D_matrix, "./ITHACAoutput/Matrices/",
            "D_" + name(NPrghmodes));
    }

    return D_matrix;
}

Eigen::Tensor<double, 3> UnsteadyBBTurb::divMomentum(label NU, label NPrghmodes)
{
    label g1Size = NPrghmodes + liftfieldP.size();
    label g2Size = NU + liftfield.size();
    Eigen::Tensor<double, 3> gTensor;
    gTensor.resize(g1Size, g2Size, g2Size);

    for (label i = 0; i < g1Size; i++)
    {
        for (label j = 0; j < g2Size; j++)
        {
            for (label k = 0; k < g2Size; k++)
            {
                gTensor(i, j, k) = fvc::domainIntegrate(fvc::grad(P_rghmodes[i]) & (fvc::div(fvc::interpolate(L_U_SUPmodes[j]) & L_U_SUPmodes[j].mesh().Sf(), L_U_SUPmodes[k]))).value();
            }
        }
    }

    if (Pstream::master())
    {
        // Export the tensor
        ITHACAstream::SaveDenseTensor(gTensor, "./ITHACAoutput/Matrices/",
            "G_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPrghmodes) + "_t");
    }

    return gTensor;
}

Eigen::MatrixXd UnsteadyBBTurb::buoyant_term_poisson(label NPrghmodes,
    label NTmodes)
{
    label H1size = NPrghmodes;
    label H2size = NTmodes + liftfieldT.size();
    Eigen::MatrixXd HP_matrix(H1size, H2size);
    // Create PTRLIST with lift, velocities and temperatures
    dimensionedScalar beta = _beta();
    dimensionedScalar TRef = _TRef();
    dimensionedVector g = _g();
    volScalarField& gh = _gh();
    surfaceScalarField& ghf = _ghf();

    // Project everything
    for (label i = 0; i < H1size; i++)
    {
        for (label j = 0; j < H2size; j++)
        {
            HP_matrix(i, j) = fvc::domainIntegrate(fvc::reconstruct(fvc::snGrad(
                                                                        P_rghmodes[i]) *
                                                       P_rghmodes[i].mesh().magSf()) &
                fvc::reconstruct(
                    ghf * fvc::snGrad(-(beta * (L_Tmodes[j]))) * L_Tmodes[j].mesh().magSf()))
                                  .value();
        }
    }

    if (Pstream::parRun())
    {
        reduce(HP_matrix, sumOp<Eigen::MatrixXd>());
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(HP_matrix, "./ITHACAoutput/Matrices/",
            "HP_" + name(NPrghmodes) + "_" + name(liftfieldT.size()) + "_" + name(NTmodes));
    }

    return HP_matrix;
}

Eigen::MatrixXd UnsteadyBBTurb::pressure_BC1(label NUmodes, label NPrghmodes)
{
    label P_BC1size = NPrghmodes;
    label P_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC1_matrix(P_BC1size, P_BC2size);
    const fvMesh& mesh = L_U_SUPmodes[0].mesh();

    for (label i = 0; i < P_BC1size; i++)
    {
        for (label j = 0; j < P_BC2size; j++)
        {
            surfaceScalarField lpl((fvc::interpolate(fvc::laplacian(
                                        L_U_SUPmodes[j])) &
                                       mesh.Sf()) *
                fvc::interpolate(P_rghmodes[i]));
            double s = 0;

            for (label k = 0; k < lpl.boundaryField().size(); k++)
            {
                if (lpl.boundaryField()[k].type() != "processor")
                {
                    s += gSum(lpl.boundaryField()[k]);
                }
            }

            BC1_matrix(i, j) = s;
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(BC1_matrix, "./ITHACAoutput/Matrices/",
            "BC1_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPrghmodes));
    }

    return BC1_matrix;
}

Eigen::Tensor<double, 3> UnsteadyBBTurb::pressure_BC2(label NUmodes, label NPrghmodes)
{
    label pressureBC1Size = NPrghmodes;
    label pressureBC2Size = NUmodes + liftfield.size();
    Eigen::Tensor<double, 3> bc2Tensor;
    const fvMesh& mesh = L_U_SUPmodes[0].mesh();
    bc2Tensor.resize(pressureBC1Size, pressureBC2Size, pressureBC2Size);

    for (label i = 0; i < pressureBC1Size; i++)
    {
        for (label j = 0; j < pressureBC2Size; j++)
        {
            for (label k = 0; k < pressureBC2Size; k++)
            {
                surfaceScalarField div_m(fvc::interpolate(fvc::div(fvc::interpolate(
                                                                       L_U_SUPmodes[j]) &
                                                 mesh.Sf(),
                                             L_U_SUPmodes[k])) &
                    mesh.Sf() * fvc::interpolate(P_rghmodes[i]));
                double s = 0;

                for (label k = 0; k < div_m.boundaryField().size(); k++)
                {
                    if (div_m.boundaryField()[k].type() != "processor")
                    {
                        s += gSum(div_m.boundaryField()[k]);
                    }
                }

                bc2Tensor(i, j, k) = s;
            }
        }
    }

    if (Pstream::master())
    {
        // Export the tensor
        ITHACAstream::SaveDenseTensor(bc2Tensor, "./ITHACAoutput/Matrices/",
            "BC2_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NSUPmodes) + "_" + name(NPrghmodes) + "_t");
    }

    return bc2Tensor;
}

Eigen::MatrixXd UnsteadyBBTurb::pressure_BC3(label NUmodes, label NPrghmodes)
{
    label P3_BC1size = NPrghmodes;
    label P3_BC2size = NUmodes + liftfield.size();
    Eigen::MatrixXd BC3_matrix(P3_BC1size, P3_BC2size);
    const fvMesh& mesh = L_U_SUPmodes[0].mesh();
    surfaceVectorField n(mesh.Sf() / mesh.magSf());

    for (label i = 0; i < P3_BC1size; i++)
    {
        for (label j = 0; j < P3_BC2size; j++)
        {
            surfaceVectorField BC3 = fvc::interpolate(fvc::curl(L_U_SUPmodes[j])).ref();
            surfaceVectorField BC4 = (n ^ fvc::interpolate(fvc::grad(P_rghmodes[i]))).ref();
            surfaceScalarField BC5 = ((BC3 & BC4) * mesh.magSf()).ref();
            double s = 0;

            for (label k = 0; k < BC5.boundaryField().size(); k++)
            {
                if (BC5.boundaryField()[k].type() != "processor")
                {
                    s += gSum(BC5.boundaryField()[k]);
                }
            }

            BC3_matrix(i, j) = s;
        }
    }

    if (Pstream::master())
    {
        ITHACAstream::SaveDenseMatrix(BC3_matrix, "./ITHACAoutput/Matrices/",
            "BC3_" + name(liftfield.size()) + "_" + name(NUmodes) + "_" + name(NPrghmodes));
    }

    return BC3_matrix;
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
        ITHACAstream::exportMatrix(bcTempVec, "bcTempVec", "python", "./ITHACAoutput/Matrices/python/");
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

    if (Pstream::master())
    {
        ITHACAstream::exportMatrix(bcTempMat, "bcTempMat", "python", "./ITHACAoutput/Matrices/python/");
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
