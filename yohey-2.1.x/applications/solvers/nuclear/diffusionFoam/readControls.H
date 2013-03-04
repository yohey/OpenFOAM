    const dictionary& neutronDiffusion =
        mesh.solutionDict().subDict("neutronDiffusion");

    const int nEnergyGroups =
        neutronDiffusion.lookupOrDefault<int>("nEnergyGroups", 1);
