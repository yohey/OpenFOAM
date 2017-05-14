
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Array<Cmpt>::Array():
    ArraySpace<Array<Cmpt>, Cmpt>()
{}


template<class Cmpt>
inline Array<Cmpt>::Array(int nCmpt):
    ArraySpace<Array<Cmpt>, Cmpt>(nCmpt)
{}


template<class Cmpt>
inline Array<Cmpt>::Array(int nCmpt, const Cmpt& v):
    ArraySpace<Array<Cmpt>, Cmpt>(nCmpt, v)
{}


template<class Cmpt>
template<class Cmpt2>
inline Array<Cmpt>::Array(const ArraySpace<Array<Cmpt2>, Cmpt2>& as):
    ArraySpace<Array<Cmpt>, Cmpt>(as)
{}


template<class Cmpt>
inline Array<Cmpt>::Array(int nCmpt, Istream& is):
    ArraySpace<Array<Cmpt>, Cmpt>(nCmpt, is)
{}


template<class Cmpt>
inline Array<Cmpt>::Array(Istream& is):
    ArraySpace<Array<Cmpt>, Cmpt>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Cmpt>
inline const Array<Cmpt>& Array<Cmpt>::centre
(
    const Foam::List<Array<Cmpt> >&
) const
{
    return *this;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline typename innerProduct<Array<Cmpt>, Array<Cmpt> >::type
operator&(const Array<Cmpt>& a1, const Array<Cmpt>& a2)
{
#   ifdef FULLDEBUG
    if (a1.nCmpt_ != a2.nCmpt_)
    {
        FatalErrorIn("operator&(Array<Cmpt>&, Array<Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    Cmpt inProd = a1.v_[0] * a2.v_[0];
    for (int i = 1; i < a1.nCmpt_; ++i)
    {
        inProd += a1.v_[i] * a2.v_[i];
    }
    return inProd;
}


// template<class Cmpt>
// inline Array<Cmpt> operator^(const Array<Cmpt>& v1, const Array<Cmpt>& v2)
// {
//     return Array<Cmpt>
//     (
//         (v1.y()*v2.z() - v1.z()*v2.y()),
//         (v1.z()*v2.x() - v1.x()*v2.z()),
//         (v1.x()*v2.y() - v1.y()*v2.x())
//     );
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
