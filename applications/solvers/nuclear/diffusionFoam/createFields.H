    Info<< "Reading field phi\n" << endl;

    volScalarField **phi = new volScalarField*[nEnergyGroups];

    char phiName[9];

    for (int group = 0; group < nEnergyGroups; group++)
    {
      sprintf(phiName, "phi%d", group);

      phi[group] = new volScalarField
                   (
                       IOobject
                       (
                           phiName,
                           runTime.timeName(),
                           mesh,
                           IOobject::MUST_READ,
                           IOobject::AUTO_WRITE
                       ),
                       mesh
                   );
    }


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

    dimensionedScalar **D        = new dimensionedScalar*[nEnergyGroups];
    dimensionedScalar **chi      = new dimensionedScalar*[nEnergyGroups];
    dimensionedScalar **SigmaA   = new dimensionedScalar*[nEnergyGroups];
    dimensionedScalar **SigmaR   = new dimensionedScalar*[nEnergyGroups];
    dimensionedScalar **SigmaF   = new dimensionedScalar*[nEnergyGroups];
    dimensionedScalar **nuSigmaF = new dimensionedScalar*[nEnergyGroups];
    dimensionedScalar ***SigmaS   = new dimensionedScalar**[nEnergyGroups];

    char DName[9];
    char chiName[9];
    char SigmaAName[9];
    char SigmaRName[9];
    char SigmaFName[9];
    char nuSigmaFName[10];
    char SigmaSName[9];

    for (int group = 0; group < nEnergyGroups; group++)
    {
      sprintf(DName,        "D%d",        group);
      sprintf(chiName,      "chi%d",      group);
      sprintf(SigmaAName,   "SigmaA%d",   group);
      sprintf(SigmaRName,   "SigmaR%d",   group);
      sprintf(SigmaFName,   "SigmaF%d",   group);
      sprintf(nuSigmaFName, "nuSigmaF%d", group);

      D[group] = new dimensionedScalar
                 (
                     DName,
                     dimLength,
                     transportProperties.lookup(DName)
                 );

      chi[group] = new dimensionedScalar
                   (
                       chiName,
                       dimless,
                       transportProperties.lookup(chiName)
                   );

      SigmaA[group] = new dimensionedScalar
                      (
                          SigmaAName,
                          dimless / dimLength,
                          transportProperties.lookup(SigmaAName)
                      );

      SigmaR[group] = new dimensionedScalar
                      (
                          SigmaRName,
                          *SigmaA[group]
                      );

      delete SigmaA[group];

      SigmaF[group] = new dimensionedScalar
                      (
                          SigmaFName,
                          dimless / dimLength,
                          transportProperties.lookup(SigmaFName)
                      );

      nuSigmaF[group] = new dimensionedScalar
                        (
                            nuSigmaFName,
                            dimless / dimLength,
                            transportProperties.lookup(nuSigmaFName)
                        );

      SigmaS[group] = new dimensionedScalar*[nEnergyGroups];

      for (int child = 0; child < nEnergyGroups; child++)
      {
          sprintf(SigmaSName, "SigmaS%d%d",   group, child);

          if (group < child)
          {
              SigmaS[group][child] = new dimensionedScalar
                                     (
                                         SigmaSName,
                                         dimless / dimLength,
                                         transportProperties.lookup(SigmaSName)
                                     );

              *SigmaR[group] += *SigmaS[group][child];
          }
          else
          {
              SigmaS[group][child] = new dimensionedScalar
                                     (
                                         SigmaSName,
                                         dimless / dimLength,
                                         0
                                     );
          }
      }
    }


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
