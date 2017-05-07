/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    diffusionFoam

Description
    Solves a multigroup diffusion equation for nuclear reactors.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readControls.H"
    #include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const dimensionedScalar energyPerFission("energyPerFission", dimEnergy, 200e6 * 1.6021892e-19);
    const dimensionedScalar targetFissions("targetFissions", power / energyPerFission);

    Info<< "\nCalculating neutron distribution\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        dimensionedScalar fissions("fissions", targetFissions.dimensions(), 0);

        volScalarField neutrons(*nuSigmaF[0] * phi[0]->oldTime());

        for (int group = 1; group < nEnergyGroups; group++)
        {
            neutrons += *nuSigmaF[group] * phi[group]->oldTime();
        }


        for (int group = nEnergyGroups-1; group >= 0; group--)
        {
            volScalarField source(*chi[group] * neutrons);

            for (int parent = 0; parent < group; parent++)
            {
                source += *SigmaS[parent][group] * phi[parent]->oldTime();
            }


            solve
            (
                -fvm::laplacian(*D[group], *phi[group])
                + fvm::Sp(*SigmaR[group], *phi[group])
                ==
                source
            );

            phi[group]->correctBoundaryConditions();

            fissions += fvc::domainIntegrate(*SigmaF[group] * *phi[group]);
        }


        kEff = fissions / targetFissions;

        Info << "kEff = " << kEff.value() << endl;

        for (int group = 0; group < nEnergyGroups; group++)
        {
            *phi[group] /= kEff;
        }


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
