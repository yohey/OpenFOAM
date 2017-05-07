
#include "ArraySpace.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from Istream
template<class Form, class Cmpt>
Foam::ArraySpace<Form, Cmpt>::ArraySpace
(
    int nCmpt,
    Istream& is
)
{
    nCmpt_ = nCmpt;
    v_ = new Cmpt[nCmpt_];

    // Read beginning of ArraySpace<Cmpt>
    is.readBegin("ArraySpace<Form, Cmpt>");

    for (int i = 0; i < nCmpt_; ++i)
    {
        is >> v_[i];
    }

    // Read end of ArraySpace<Cmpt>
    is.readEnd("ArraySpace<Form, Cmpt>");

    // Check state of Istream
    is.check("ArraySpace<Form, Cmpt>::ArraySpace(Istream&)");
}


// Construct from Istream
template<class Form, class Cmpt>
Foam::ArraySpace<Form, Cmpt>::ArraySpace
(
    Istream& is
)
{
    is >> nCmpt_;

    v_ = new Cmpt[nCmpt_];

    // Read beginning of ArraySpace<Cmpt>
    is.readBegin("ArraySpace<Form, Cmpt>");

    for (int i = 0; i < nCmpt_; ++i)
    {
        is >> v_[i];
    }

    // Read end of ArraySpace<Cmpt>
    is.readEnd("ArraySpace<Form, Cmpt>");

    // Check state of Istream
    is.check("ArraySpace<Form, Cmpt>::ArraySpace(Istream&)");
}


// Return a string representation
template<class Form, class Cmpt>
Foam::word
Foam::name
(
    const ArraySpace<Form, Cmpt>& as
)
{
    std::ostringstream buf;

    buf << '(';

    for (int i = 0; i < as.nCmpt_-1; ++i)
    {
        buf << as.v_[i] << ',';
    }

    buf << as.v_[as.nCmpt_-1] << ')';

    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Form, class Cmpt>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    ArraySpace<Form, Cmpt>& as
)
{
    token nextToken(is);
    is.putBack(nextToken);

    if (nextToken.isNumber())
    {
      int nCmpt;
      is >> nCmpt;

      if (nCmpt != as.nCmpt_)
      {
        FatalErrorIn("operator>>(Istream&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
      }
    }

    // Read beginning of ArraySpace<Cmpt>
    is.readBegin("ArraySpace<Form, Cmpt>");

    for (int i = 0; i < as.size(); ++i)
    {
        is >> as.v_[i];
    }

    // Read end of ArraySpace<Cmpt>
    is.readEnd("ArraySpace<Form, Cmpt>");

    // Check state of Istream
    is.check("operator>>(Istream&, ArraySpace<Form, Cmpt>&)");

    return is;
}


template<class Form, class Cmpt>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ArraySpace<Form, Cmpt>& as
)
{
    os << as.nCmpt_ << token::SPACE;

    os << token::BEGIN_LIST;

    for (int i = 0; i < as.size()-1; ++i)
    {
        os << as.v_[i] << token::SPACE;
    }

    os << as.v_[as.size()-1] << token::END_LIST;

    // Check state of Ostream
    os.check("operator<<(Ostream&, const ArraySpace<Form, Cmpt>&)");

    return os;
}


// ************************************************************************* //
