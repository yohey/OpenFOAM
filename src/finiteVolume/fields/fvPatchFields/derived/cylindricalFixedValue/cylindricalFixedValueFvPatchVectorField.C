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

#include "cylindricalFixedValueFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylindricalFixedValueFvPatchVectorField::
cylindricalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    origin_(vector::zero),
    axis_(vector::zero),
    direction_(p.size()),
    magnitude_(p.size())
{}


Foam::cylindricalFixedValueFvPatchVectorField::
cylindricalFixedValueFvPatchVectorField
(
    const cylindricalFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(p, iF),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    direction_(ptf.direction_, mapper),
    magnitude_(ptf.magnitude_, mapper)
{
    vectorField radialCf(ptf.patch().Cf() - (ptf.patch().Cf() & ptf.axis_) * ptf.axis_);
    vectorField e1(radialCf / mag(radialCf));
    vectorField e2((ptf.axis_ ^ radialCf) / mag(ptf.axis_ ^ radialCf));
    vector e3 = (ptf.axis_) / mag(ptf.axis_);

    vectorField cartDir(ptf.direction_.component(vector::X) * e1
                        + ptf.direction_.component(vector::Y) * e2
                        + ptf.direction_.component(vector::Z) * e3);

    // Note: calculate product only on ptf to avoid multiplication on
    // unset values in reconstructPar.
    fixedValueFvPatchVectorField::operator=
    (
        vectorField
        (
            ptf.magnitude_ * cartDir / mag(cartDir),
            mapper
        )
    );
}


Foam::cylindricalFixedValueFvPatchVectorField::
cylindricalFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    direction_("direction", dict, p.size()),
    magnitude_("magnitude", dict, p.size())
{
    vectorField radialCf(patch().Cf() - (patch().Cf() & axis_) * axis_);
    vectorField e1(radialCf / mag(radialCf));
    vectorField e2((axis_ ^ radialCf) / mag(axis_ ^ radialCf));
    vector e3 = (axis_) / mag(axis_);

    vectorField cartDir(direction_.component(vector::X) * e1
                        + direction_.component(vector::Y) * e2
                        + direction_.component(vector::Z) * e3);

    fvPatchVectorField::operator=(magnitude_ * cartDir / mag(cartDir));
}


Foam::cylindricalFixedValueFvPatchVectorField::
cylindricalFixedValueFvPatchVectorField
(
    const cylindricalFixedValueFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    direction_(pivpvf.direction_),
    magnitude_(pivpvf.magnitude_)
{}


Foam::cylindricalFixedValueFvPatchVectorField::
cylindricalFixedValueFvPatchVectorField
(
    const cylindricalFixedValueFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    direction_(pivpvf.direction_),
    magnitude_(pivpvf.magnitude_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylindricalFixedValueFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);
    direction_.autoMap(m);
    magnitude_.autoMap(m);
}


void Foam::cylindricalFixedValueFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const cylindricalFixedValueFvPatchVectorField& tiptf =
        refCast<const cylindricalFixedValueFvPatchVectorField>(ptf);

    direction_.rmap(tiptf.direction_, addr);
    magnitude_.rmap(tiptf.magnitude_, addr);
}


void Foam::cylindricalFixedValueFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("origin") << origin_ << token::END_STATEMENT << nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    direction_.writeEntry("direction", os);
    magnitude_.writeEntry("magnitude", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        cylindricalFixedValueFvPatchVectorField
    );
}

// ************************************************************************* //
