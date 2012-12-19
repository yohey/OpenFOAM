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

Class
    Foam::cylindricalFixedValueFvPatchVectorField

Description
    Describes a surface normal vector boundary condition by its magnitude.
    Note: The value is positive for outward-pointing vectors

SourceFiles
    cylindricalFixedValueFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef cylindricalFixedValueFvPatchVectorField_H
#define cylindricalFixedValueFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class cylindricalFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class cylindricalFixedValueFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        vector origin_;
        vector axis_;
        vectorField direction_;
        scalarField magnitude_;


public:

    //- Runtime type information
    TypeName("cylindricalFixedValue");


    // Constructors

        //- Construct from patch and internal field
        cylindricalFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        cylindricalFixedValueFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  cylindricalFixedValueFvPatchVectorField
        //  onto a new patch
        cylindricalFixedValueFvPatchVectorField
        (
            const cylindricalFixedValueFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        cylindricalFixedValueFvPatchVectorField
        (
            const cylindricalFixedValueFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new cylindricalFixedValueFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        cylindricalFixedValueFvPatchVectorField
        (
            const cylindricalFixedValueFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new cylindricalFixedValueFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }



    // Member functions

        // Access

            //- Return the origin of the cylindricalCS
            const vector& origin() const
            {
                return origin_;
            }

            //- Return the axis of the cylindricalCS
            const vector& axis() const
            {
                return axis_;
            }

            //- Return non-const access to the origin of the cylindricalCS
            vector& origin()
            {
                return origin_;
            }

            //- Return non-const access to the axis of the cylindricalCS
            vector& axis()
            {
                return axis_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
