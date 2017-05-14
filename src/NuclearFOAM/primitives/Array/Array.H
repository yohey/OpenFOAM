
#ifndef Array_H
#define Array_H

#include "ArraySpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T> class List;

/*---------------------------------------------------------------------------*\
                           Class Array Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Array:
    public ArraySpace<Array<Cmpt>, Cmpt>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef Array<label> labelType;


    // Member constants

        //- Rank of Array is 1
        static const direction rank = 1;


    // Static data members

        static const char* const typeName;
        static const Array zero;
        static const Array one;
        static const Array max;
        static const Array min;
        static const Array rootMax;
        static const Array rootMin;


    // Static member functions

        inline static const Array zeros(int);
        inline static const Array ones(int);
        // inline static const Array max(int);
        // inline static const Array min(int);
        // inline static const Array rootMax(int);
        // inline static const Array rootMin(int);


    // Constructors

        //- Construct null
        inline Array();

        //- Construct null
        inline Array(int);

        //- Construct from an uniform value
        inline Array(int, const Cmpt&);

        //- Construct given ArraySpace of the same rank
        template<class Cmpt2>
        inline Array(const ArraySpace<Array<Cmpt2>, Cmpt2>&);

        //- Construct from Istream
        inline Array(int, Istream&);

        //- Construct from Istream
        inline Array(Istream&);


    // Member Functions

        //- Return *this (used for point which is a typedef to Array<scalar>.
        inline const Array<Cmpt>& centre
        (
            const Foam::List<Array<Cmpt> >&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ArrayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
