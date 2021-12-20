    Info << "Reading thermopysical properties\n" << endl;

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h");

    volScalarField& p = thermo.p();
    volScalarField& h0 = thermo.he();
    const volScalarField& psi = thermo.psi();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "compressibleCreatePhi.H"

    U.correctBoundaryConditions();

    #include "PhiBoundaryTypes.H"

    volScalarField Phi
    (
        IOobject
        (
            "Phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", U.dimensions()*dimLength, 0.0),
        PhiBoundaryTypes
    );

    p.correctBoundaryConditions();
    h0.correctBoundaryConditions();
    thermo.correct();

    volScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        h0 - magSqr(U)/2.
    );

    if (p0.value() == 0)
    {
        if (T0.value() == 0)
        {
            p0 = max(p *
                     pow(
                         h0 / h,
                         thermo.Cp() / (thermo.Cp() - thermo.Cv())
                         )
                     );
        }
        else
        {
            p0 = max(p *
                     pow(
                         thermo.Cp() * T0 / h,
                         thermo.Cp() / (thermo.Cp() - thermo.Cv())
                         )
                     );
        }
    }

    label PhiRefCell = 0;
    scalar PhiRefValue = 0.0;

    const dimensionedScalar hMin("hMin", dimEnergy/dimMass, 1e-6);

    #include "createFvOptions.H"
