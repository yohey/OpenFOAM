
#ifndef Array2_H
#define Array2_H

#include <new>

#include "Array.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Array2 Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Array2:
    public Array<Array<Cmpt>>
{

public:

    //- Equivalent type of labels used for valid component indexing
    typedef Array2<label> labelType;


    // Member constants

        //- Rank of Array2 is 2
        static const direction rank = 2;


    // Static data members

        static const char* const typeName;
        static const Array2 zero;
        static const Array2 one;
        static const Array2 max;
        static const Array2 min;
        static const Array2 rootMax;
        static const Array2 rootMin;


    // Static member functions

        inline static const Array2 zeros(int, int);
        inline static const Array2 ones(int, int);
        // inline static const Array2 max(int, int);
        // inline static const Array2 min(int, int);
        // inline static const Array2 rootMax(int, int);
        // inline static const Array2 rootMin(int, int);


    // Constructors

        //- Construct null
        inline Array2(int, int);

        //- Construct from an uniform value
        inline Array2(int, int, const Cmpt&);

        //- Construct as copy
        inline Array2(const Array2<Cmpt>&);

        //- Construct as copy of another Array2 type
        template<class Cmpt2>
        inline Array2(const Array2<Cmpt2>&);

        //- Construct from Istream
        inline Array2(int, int, Istream&);

        //- Construct from Istream
        inline Array2(Istream&);


    // Member Functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Array2I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
