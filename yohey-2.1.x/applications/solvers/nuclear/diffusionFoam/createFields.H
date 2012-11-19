    Info<< "Reading field phi\n" << endl;

    volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    dimensionedScalar kEff
    (
        "kEff",
        dimless,
        1.
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    dimensionedScalar D
    (
        "D",
        dimLength,
        transportProperties.lookup("D")
    );

    dimensionedScalar sigmaA
    (
        "sigmaA",
        dimless / dimLength,
        transportProperties.lookup("sigmaA")
    );

    dimensionedScalar sigmaF
    (
        "sigmaF",
        dimless / dimLength,
        transportProperties.lookup("sigmaF")
    );

    dimensionedScalar nuSigmaF
    (
        "nuSigmaF",
        dimless / dimLength,
        transportProperties.lookup("nuSigmaF")
    );


    Info<< "Reading operationProperties\n" << endl;

    IOdictionary operationProperties
    (
        IOobject
        (
            "operationProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar power
    (
        "power",
        dimPower,
        operationProperties.lookup("power")
    );
