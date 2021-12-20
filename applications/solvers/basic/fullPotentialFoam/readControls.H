    const dictionary& fullPotentialFlow =
        mesh.solutionDict().subDict("fullPotentialFlow");

    const int nNonOrthCorr =
        fullPotentialFlow.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const dimensionedScalar T0
    (
        "T0",
        dimTemperature,
        fullPotentialFlow.lookupOrDefault<double>("T0", 0.0)
    );

    dimensionedScalar p0
    (
        "p0",
        dimPressure,
        fullPotentialFlow.lookupOrDefault<double>("p0", 0.0)
    );
