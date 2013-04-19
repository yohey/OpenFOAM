const volVectorField::GeometricBoundaryField& ubf = U.boundaryField();
wordList PhiBoundaryTypes = ubf.types();

forAll(PhiBoundaryTypes, patchi)
{
    if(ubf[patchi].fixesValue())
    {
        PhiBoundaryTypes[patchi] = fixedGradientFvPatchScalarField::typeName;
    }
    else if (PhiBoundaryTypes[patchi] == "zeroGradient")
    {
        PhiBoundaryTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}
