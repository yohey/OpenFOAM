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

    dimensionedScalar SigmaA
    (
        "SigmaA",
        dimless / dimLength,
        transportProperties.lookup("SigmaA")
    );

    dimensionedScalar SigmaF
    (
        "SigmaF",
        dimless / dimLength,
        transportProperties.lookup("SigmaF")
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
