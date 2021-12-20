/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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
    fullPotentialFoam

Description
    Full potential flow solver which can be used to generate starting fields
    for full Navier-Stokes codes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "psiThermo.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readControls.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Calculating potential flow" << endl;

    while (runTime.loop())
    {
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "setPhiBoundary.H"

        volScalarField aSqr("aSqr", (thermo.Cp()/thermo.Cv()) / psi);

        phi = linearInterpolate(rho*U) & mesh.Sf();

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix PhiEqn
            (
                fvm::laplacian(Phi)
             ==
                ((U * U) && fvc::grad(U) / aSqr)
            );

            PhiEqn.setReference(PhiRefCell, PhiRefValue);
            PhiEqn.solve();
        }

        scalar alpha = 0.1;

        U = alpha * fvc::grad(Phi) + (1 - alpha) * U;
        U.correctBoundaryConditions();
        fvOptions.correct(U);

        if (T0.value() == 0)
        {
            h = h0 - magSqr(U)/2.;
        }
        else
        {
            h = thermo.Cp() * T0 - magSqr(U)/2.;
        }

        h = pos(h) * h + (1 - pos(h)) * hMin;

        if (T0.value() == 0)
        {
            p = p0 /
                pow(
                    h0 / h,
                    thermo.Cp() / (thermo.Cp() - thermo.Cv())
                   );
        }
        else
        {
            p = p0 /
                pow(
                    thermo.Cp() * T0 / h,
                    thermo.Cp() / (thermo.Cp() - thermo.Cv())
                   );
        }

        h0.correctBoundaryConditions();
        p.correctBoundaryConditions();
        thermo.correct();

        rho = thermo.rho();

        runTime.write();

        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
