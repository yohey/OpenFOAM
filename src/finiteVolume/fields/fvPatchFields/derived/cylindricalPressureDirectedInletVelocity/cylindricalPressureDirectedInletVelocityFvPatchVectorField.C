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

\*---------------------------------------------------------------------------*/

#include "cylindricalPressureDirectedInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::
cylindricalPressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_("phi"),
    rhoName_("rho"),
    origin_(vector::zero),
    axis_(vector::zero),
    inletDir_(p.size())
{}


Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::
cylindricalPressureDirectedInletVelocityFvPatchVectorField
(
    const cylindricalPressureDirectedInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    inletDir_(ptf.inletDir_, mapper)
{}


Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::
cylindricalPressureDirectedInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    inletDir_("inletDirection", dict, p.size())
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::
cylindricalPressureDirectedInletVelocityFvPatchVectorField
(
    const cylindricalPressureDirectedInletVelocityFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    inletDir_(pivpvf.inletDir_)
{}


Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::
cylindricalPressureDirectedInletVelocityFvPatchVectorField
(
    const cylindricalPressureDirectedInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    phiName_(pivpvf.phiName_),
    rhoName_(pivpvf.rhoName_),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    inletDir_(pivpvf.inletDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    inletDir_.autoMap(m);
}


void Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const cylindricalPressureDirectedInletVelocityFvPatchVectorField& tiptf =
        refCast<const cylindricalPressureDirectedInletVelocityFvPatchVectorField>(ptf);

    inletDir_.rmap(tiptf.inletDir_, addr);
}


void Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    vectorField radialCf(patch().Cf() - (patch().Cf() & axis_) * axis_);
    vectorField e1(radialCf / mag(radialCf));
    vectorField e2((axis_ ^ radialCf) / mag(axis_ ^ radialCf));
    vector e3 = (axis_) / mag(axis_);

    cartesianInletDir =
      (inletDir_.component(vector::X) * e1
       + inletDir_.component(vector::Y) * e2
       + inletDir_.component(vector::Z) * e3)
      / mag(inletDir_);

    tmp<vectorField> n = patch().nf();
    tmp<scalarField> ndmagS = (n & cartesianInletDir) * patch().magSf();

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        operator==(cartesianInletDir * phip/ndmagS);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        operator==(cartesianInletDir * phip/(rhop*ndmagS));
    }
    else
    {
        FatalErrorIn
        (
            "cylindricalPressureDirectedInletVelocityFvPatchVectorField::updateCoeffs()"
        )   << "dimensions of phi are not correct"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    inletDir_.writeEntry("inletDirection", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cylindricalPressureDirectedInletVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& pvf
)
{
    fvPatchField<vector>::operator=(cartesianInletDir * (cartesianInletDir & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        cylindricalPressureDirectedInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
