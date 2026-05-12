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

  Description
    In this tutorial we recreate the Example 1 in Mendez et al. "Multi-Scale Proper Orthogonal Decomposition of Complex Fluid Flows".
    The dataset is generated on a square cartesian grid of ns = 256*256 points over a domain [20 x 20], while the time domain spans tk [0, 5.11] with a total nt =256 and sampling frequency fs = 100.
    The dataset is composed of the sum of three modes, being identical gaussians of sigma=5, located in x1 = [10, -10], x2 = [-10, 10] and x3 = [0, 0], and pulsing as:
      T1(tk) = A1*sin(2*pi*f1*tk) * exp(0.5(tk -3)^20)
      T2(tk) = A2*sin(2*pi*f2*tk-pi/3)
      T3(tk) = A3*sin(2*pi*f3*tk-pi/4)*(tk - 2.55)^2
    The amplitudes are chosen so as to equal each contribution, that is:
      Ar = 1/norm(Sr Tr^T)
    where Sr is the r-th spatial mode and Tr is the r-th temporal mode reshaped as a column.
    The three frequencies are f1=15, f2=0.1 f3=7
    To adapt with ITHACA, the data generated is stored in OpenFOAM format.
    The mesh is created using blockMesh, while the data is created in C++/Eigen and then converted to OpenFOAM volScalarField format using the Foam2Eigen ITHACA class.

    SourceFiles
    mPOD.C

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "Foam2Eigen.H"
#include "ITHACAPOD.H"

// Function to get the 2D gaussian as eigen for each timestep
Eigen::MatrixXd gaussian2D(const Eigen::MatrixXd& grid, const Eigen::Vector2d& center, double sigma)
{
    Eigen::MatrixXd mode(grid.rows(), 1);
    for (int i = 0; i < grid.rows(); ++i)
    {
        double x = grid(i, 0);
        double y = grid(i, 1);
        mode(i, 0) = std::exp(-((x - center(0)) * (x - center(0)) + (y - center(1)) * (y - center(1))) / (2 * sigma * sigma));
    }
    return mode;
}

// Function to get spatial modes
Eigen::MatrixXd spatialModes(const Eigen::MatrixXd& cellCentres)
{
    Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(cellCentres.rows(), 3);
    Eigen::Vector2d center1(10, -10);
    Eigen::Vector2d center2(-10, 10);
    Eigen::Vector2d center3(0, 0);
    double sigma = 5.0;
    Eigen::MatrixXd mode1 = gaussian2D(cellCentres, center1, sigma);
    Eigen::MatrixXd mode2 = gaussian2D(cellCentres, center2, sigma);
    Eigen::MatrixXd mode3 = gaussian2D(cellCentres, center3, sigma);
    // Fill the phi matrix with the three modes
    phi.col(0) = mode1;
    phi.col(1) = mode2;
    phi.col(2) = mode3;
    return phi;
}

// Function to get the temporal modes
Eigen::MatrixXd temporalModes(const Eigen::VectorXd& time)
{
    Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(3, time.size());
    double f1 = 15.0;
    double f2 = 0.1;
    double f3 = 7.0;
    for (int i = 0; i < time.size(); ++i)
    {
        double tk = time(i);
        psi(0, i) = std::sin(2 * M_PI * f1 * tk) * std::exp(-0.5 * std::pow((tk - 3), 20));
        psi(1, i) = std::sin(2 * M_PI * f2 * tk - M_PI / 3);
        psi(2, i) = std::sin(2 * M_PI * f3 * tk) * std::pow(tk - 2.55, 2);
    }
    return psi;
}

// Function to get the amplitudes
Eigen::MatrixXd computeAmplitudes(const Eigen::MatrixXd& phi, const Eigen::MatrixXd& psi)
{
    Eigen::MatrixXd sigmaMat = Eigen::MatrixXd::Zero(3, 3);
    for (int i = 0; i < 3; ++i)
    {
        double norm = (phi.col(i) * psi.row(i)).norm();
        sigmaMat(i, i) = 1.0 / norm;
    }
    return sigmaMat;
}
class mPODTest
{
public:
    mPODTest(Foam::fvMesh& mesh, Foam::Time& runTime): mesh_(mesh), runTime_(runTime) { }

    Foam::fvMesh& mesh() { return mesh_; }
    Foam::Time& runTime() { return runTime_; }

    PtrList<volScalarField> snapshots;
    PtrList<volScalarField> modes;

private:
    Foam::fvMesh& mesh_;
    Foam::Time& runTime_;
};

int main(int argc, char* argv[])
{
#include "setRootCase.H"
#include "createTime.H" // creates `runTime`
#include "createMesh.H" // creates `mesh`

    mPODTest example(mesh, runTime);
    ITHACAparameters *para = ITHACAparameters::getInstance(mesh, runTime);
    const int ns = example.mesh().nCells();
    const int nt = 256;
    // Set time start, end and time step
    example.runTime().setTime(0, 0);
    example.runTime().setEndTime(5.11);
    example.runTime().setDeltaT(5.11 / nt);

    Eigen::MatrixXd cellCentres(ns, 3);
    // Get cell centres from the mesh
    for (int i = 0; i < ns; ++i)
    {
        cellCentres(i, 0) = example.mesh().C()[i].x();
        cellCentres(i, 1) = example.mesh().C()[i].y();
        cellCentres(i, 2) = example.mesh().C()[i].z();
    }

    Eigen::VectorXd time = Eigen::VectorXd::LinSpaced(nt, 0, 5.11);

    Eigen::MatrixXd phi = Eigen::MatrixXd::Zero(ns, 3);
    Eigen::MatrixXd psi = Eigen::MatrixXd::Zero(3, nt);
    Eigen::MatrixXd sigmaMat = Eigen::MatrixXd::Zero(3, 3);

    phi = spatialModes(cellCentres);
    psi = temporalModes(time);
    sigmaMat = computeAmplitudes(phi, psi);
    
    // Generate the data as the product of the three matrices
    Eigen::MatrixXd data(ns, nt);
    data = phi * sigmaMat * psi;

    Foam::volScalarField datasetFoam(
        Foam::IOobject(
            "u",
            example.runTime().timeName(),
            example.mesh(),
            Foam::IOobject::NO_READ,
            Foam::IOobject::AUTO_WRITE),
        example.mesh(),
        Foam::dimensionedScalar("u", Foam::dimless, 0.0));

    // We now loop through time and fill the datasetFoam with the data generated, then write it to file
    // The translation from Eigen to Foam is done using the Foam2Eigen::Eigen2field, which takes Eigen::VectorXd and writes into a GeometricField
    for (int i = 0; i < nt; ++i)
    {
        example.runTime().setTime(time(i), i);
        Eigen::VectorXd dataCol = data.col(i);
        datasetFoam = Foam2Eigen::Eigen2field(datasetFoam, dataCol, true);
        
        example.snapshots.append(datasetFoam.clone());
        datasetFoam.write();
    }

    Info << "Computing POD modes..." << endl;
    ITHACAPOD::getModes(example.snapshots, example.modes, "u", false, false, false, 5, false);
}