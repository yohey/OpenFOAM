
#ifndef ArraySpace_H
#define ArraySpace_H

#include "direction.H"
#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Form, class Cmpt> class ArraySpace;

template<class Form, class Cmpt>
Istream& operator>>
(
    Istream&,
    ArraySpace<Form, Cmpt>&
);

template<class Form, class Cmpt>
Ostream& operator<<
(
    Ostream&,
    const ArraySpace<Form, Cmpt>&
);


/*---------------------------------------------------------------------------*\
                           Class ArraySpace Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Cmpt>
class ArraySpace
{

public:

    //- Component type
    typedef Cmpt cmptType;


    // Static data members

        int nCmpt_;

        //- The components of this array space
        Cmpt *v_;


    // Constructors

        //- Construct null
        inline ArraySpace(int);

        //- Construct from an uniform value
        inline ArraySpace(int, const Cmpt&);

        //- Construct from Istream
        ArraySpace(int, Istream&);

        //- Construct from Istream
        ArraySpace(Istream&);

        //- Construct as copy
        inline ArraySpace(const ArraySpace<Form, Cmpt>&);

        //- Construct as copy of another ArraySpace type of the same rank
        template<class Form2, class Cmpt2>
        inline ArraySpace(const ArraySpace<Form2, Cmpt2>&);


    // Member Functions

        //- Return the number of elements in the ArraySpace = nCmpt.
        inline label size() const;

        inline const Cmpt& component(const direction) const;
        inline Cmpt& component(const direction);

        inline void component(Cmpt&, const direction) const;
        inline void replace(const direction, const Cmpt&);


    // Member Operators

        inline const Cmpt& operator[](const direction) const;
        inline Cmpt& operator[](const direction);

        inline void operator=(const ArraySpace<Form, Cmpt>&);
        inline void operator+=(const ArraySpace<Form, Cmpt>&);
        inline void operator-=(const ArraySpace<Form, Cmpt>&);

        inline void operator*=(const scalar);
        inline void operator/=(const scalar);


    // IOstream Operators

        friend Istream& operator>> <Form, Cmpt>
        (
            Istream&,
            ArraySpace<Form, Cmpt>&
        );

        friend Ostream& operator<< <Form, Cmpt>
        (
            Ostream&,
            const ArraySpace<Form, Cmpt>&
        );
};


// * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * * //

//- Return a string representation of a ArraySpace
template<class Form, class Cmpt>
word name(const ArraySpace<Form, Cmpt>&);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ArraySpaceI.H"

#ifdef NoRepository
#   include "ArraySpace.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
