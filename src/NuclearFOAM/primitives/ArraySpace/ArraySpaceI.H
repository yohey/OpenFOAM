
#include "error.H"
#include "products.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Cmpt>
inline ArraySpace<Form, Cmpt>::ArraySpace()
{}


template<class Form, class Cmpt>
inline ArraySpace<Form, Cmpt>::ArraySpace(int nCmpt)
{
    nCmpt_ = nCmpt;
    v_ = new Cmpt[nCmpt_];
}


template<class Form, class Cmpt>
inline ArraySpace<Form, Cmpt>::ArraySpace
(
    int nCmpt,
    const Cmpt& v
)
{
    nCmpt_ = nCmpt;
    v_ = new Cmpt[nCmpt_];

    for (int i = 0; i < nCmpt_; ++i)
    {
      v_[i] = v;
    }
}


template<class Form, class Cmpt>
inline ArraySpace<Form, Cmpt>::ArraySpace
(
    const ArraySpace<Form, Cmpt>& as
)
{
    nCmpt_ = as.nCmpt_;
    v_ = new Cmpt[nCmpt_];

    for (int i = 0; i < nCmpt_; ++i)
    {
      v_[i] = as.v_[i];
    }
}


template<class Form, class Cmpt>
template<class Form2, class Cmpt2>
inline ArraySpace<Form, Cmpt>::ArraySpace
(
    const ArraySpace<Form2, Cmpt2>& as
)
{
    nCmpt_ = as.nCmpt_;
    v_ = new Cmpt[nCmpt_];

    for (int i = 0; i < nCmpt_; ++i)
    {
      v_[i] = as.v_[i];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt>
inline label ArraySpace<Form, Cmpt>::size() const
{
    return nCmpt_;
}


template<class Form, class Cmpt>
inline const Cmpt& ArraySpace<Form, Cmpt>::component
(
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt_)
    {
        FatalErrorIn
        (
            "ArraySpace<Form, Cmpt>::component(direction) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt>
inline Cmpt& ArraySpace<Form, Cmpt>::component
(
    const direction d
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt_)
    {
        FatalErrorIn("ArraySpace<Form, Cmpt>::component(direction)")
            << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt>
inline void ArraySpace<Form, Cmpt>::component
(
    Cmpt& c,
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt_)
    {
        FatalErrorIn
        (
            "ArraySpace<Form, Cmpt>::component(Cmpt&, direction) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    c = v_[d];
}


template<class Form, class Cmpt>
inline void ArraySpace<Form, Cmpt>::replace
(
    const direction d,
    const Cmpt& c
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt_)
    {
        FatalErrorIn
        (
            "ArraySpace<Form, Cmpt>::"
            "replace(direction, const Cmpt&) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    v_[d] = c;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt>
inline const Cmpt& ArraySpace<Form, Cmpt>::operator[]
(
    const direction d
) const
{
#   ifdef FULLDEBUG
    if (d >= nCmpt_)
    {
        FatalErrorIn
        (
            "ArraySpace<Form, Cmpt>::operator[](direction d) const"
        )   << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt>
inline Cmpt& ArraySpace<Form, Cmpt>::operator[]
(
    const direction d
)
{
#   ifdef FULLDEBUG
    if (d >= nCmpt_)
    {
        FatalErrorIn("ArraySpace<Form, Cmpt>::operator[](direction d)")
            << "index out of range"
            << abort(FatalError);
    }
#   endif

    return v_[d];
}


template<class Form, class Cmpt>
inline void ArraySpace<Form, Cmpt>::operator=
(
    const ArraySpace<Form, Cmpt>& as
)
{
#   ifdef FULLDEBUG
    if (nCmpt_ != as.nCmpt_)
    {
        FatalErrorIn("ArraySpace<Form, Cmpt>::operator=(ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    for (int i = 0; i < nCmpt_; ++i)
    {
      v_[i] = as.v_[i];
    }
}


template<class Form, class Cmpt>
inline void ArraySpace<Form, Cmpt>::operator+=
(
    const ArraySpace<Form, Cmpt>& as
)
{
#   ifdef FULLDEBUG
    if (nCmpt_ != as.nCmpt_)
    {
        FatalErrorIn("ArraySpace<Form, Cmpt>::operator+=(ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    for (int i = 0; i < nCmpt_; ++i)
    {
      v_[i] += as.v_[i];
    }
}


template<class Form, class Cmpt>
inline void ArraySpace<Form, Cmpt>::operator-=
(
    const ArraySpace<Form, Cmpt>& as
)
{
#   ifdef FULLDEBUG
    if (nCmpt_ != as.nCmpt_)
    {
        FatalErrorIn("ArraySpace<Form, Cmpt>::operator-=(ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    for (int i = 0; i < nCmpt_; ++i)
    {
      v_[i] -= as.v_[i];
    }
}


template<class Form, class Cmpt>
inline void ArraySpace<Form, Cmpt>::operator*=
(
    const scalar s
)
{
    for (int i = 0; i < nCmpt_; ++i)
    {
      v_[i] *= s;
    }
}


template<class Form, class Cmpt>
inline void ArraySpace<Form, Cmpt>::operator/=
(
    const scalar s
)
{
    for (int i = 0; i < nCmpt_; ++i)
    {
      v_[i] /= s;
    }
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Cmpt>
inline Cmpt& setComponent
(
    ArraySpace<Form, Cmpt>& as,
    const direction d
)
{
    return as.component(d);
}


template<class Form, class Cmpt>
inline const Cmpt& component
(
    const ArraySpace<Form, Cmpt>& as,
    const direction d
)
{
    return as.component(d);
}


// Powers of a Form
// Equivalent to outer-products between the Form and itself
// Form^0 = 1.0
template<class Form, class Cmpt>
inline typename powProduct<Form, 0>::type pow
(
    const ArraySpace<Form, Cmpt>&,
    typename powProduct<Form, 0>::type
    = pTraits<typename powProduct<Form, 0>::type>::zero
)
{
    return 1.0;
}

// // Form^1 = Form
// template<class Form, class Cmpt>
// inline typename powProduct<Form, 1>::type pow
// (
//     const ArraySpace<Form, Cmpt>& v,
//     typename powProduct<Form, 1>::type
//     = pTraits<typename powProduct<Form, 1>::type>::zero
// )
// {
//     return static_cast<const Form&>(v);
// }


// // Form^2 = sqr(Form)
// template<class Form, class Cmpt>
// inline typename powProduct<Form, 2>::type pow
// (
//     const ArraySpace<Form, Cmpt>& v,
//     typename powProduct<Form, 2>::type
//   = pTraits<typename powProduct<Form, 2>::type>::zero
// )
// {
//     return sqr(static_cast<const Form&>(v));
// }


template<class Form, class Cmpt>
inline scalar magSqr
(
    const ArraySpace<Form, Cmpt>& as
)
{
    scalar ms = magSqr(as.v_[0]);
    for (int i = 1; i < as.nCmpt_; ++i)
    {
        ms += magSqr(as.v_[i]);
    }
    return ms;
}


template<class Form, class Cmpt>
inline scalar mag
(
    const ArraySpace<Form, Cmpt>& as
)
{
    return ::sqrt(magSqr(static_cast<const Form&>(as)));
}


template<class Form, class Cmpt>
inline ArraySpace<Form, Cmpt> cmptMultiply
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
    Form v(as1.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = as1.v_[i] * as2.v_[i];
    }
    return v;
}


template<class Form, class Cmpt>
inline ArraySpace<Form, Cmpt> cmptDivide
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
    Form v(as1.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = as1.v_[i] / as2.v_[i];
    }
    return v;
}


template<class Form, class Cmpt>
inline ArraySpace<Form, Cmpt> stabilise
(
    const ArraySpace<Form, Cmpt>& as,
    const Cmpt& small
)
{
    Form v(as.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = stabilise(as.v_[i], small);
    }
    return v;
}


template<class Form, class Cmpt>
inline Cmpt cmptMax
(
    const ArraySpace<Form, Cmpt>& as
)
{
    Cmpt cMax = as.v_[0];
    for (int i = 1; i < as.nCmpt_; ++i)
    {
      cMax = max(cMax, as.v_[i]);
    }
    return cMax;
}


template<class Form, class Cmpt>
inline Cmpt cmptMin
(
    const ArraySpace<Form, Cmpt>& as
)
{
    Cmpt cMin = as.v_[0];
    for (int i = 1; i < as.nCmpt_; ++i)
    {
      cMin = min(cMin, as.v_[i]);
    }
    return cMin;
}


template<class Form, class Cmpt>
inline Cmpt cmptSum
(
    const ArraySpace<Form, Cmpt>& as
)
{
    Cmpt sum = as.v_[0];
    for (int i = 1; i < as.nCmpt_; ++i)
    {
      sum += as.v_[i];
    }
    return sum;
}


template<class Form, class Cmpt>
inline Cmpt cmptAv
(
    const ArraySpace<Form, Cmpt>& as
)
{
    return cmptSum(as) / as.nCmpt_;
}


template<class Form, class Cmpt>
inline Cmpt cmptProduct
(
    const ArraySpace<Form, Cmpt>& as
)
{
    Cmpt product = as.v_[0];
    for (int i = 1; i < as.nCmpt_; ++i)
    {
      product *= as.v_[i];
    }
    return product;
}


template<class Form, class Cmpt>
inline Form cmptMag
(
    const ArraySpace<Form, Cmpt>& as
)
{
    Form v(as.nCmpt_);
    for (int i = 0; i < as.nCmpt_; ++i)
    {
      v.v_[i] = mag(as.v_[i]);
    }
    return v;
}


template<class Form, class Cmpt>
inline Form max
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
    Form v(as1.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = max(as1.v_[i], as2.v_[i]);
    }
    return v;
}


template<class Form, class Cmpt>
inline Form min
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
    Form v(as1.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = min(as1.v_[i], as2.v_[i]);
    }
    return v;
}


template<class Form, class Cmpt>
inline Form minMod
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
    Form v(as1.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = minMod(as1.v_[i], as2.v_[i]);
    }
    return v;
}


template
<
    class Form1, class Cmpt1,
    class Form2, class Cmpt2
>
inline typename innerProduct<Form1, Form2>::type dot
(
    const ArraySpace<Form1, Cmpt1>& t1,
    const ArraySpace<Form2, Cmpt2>& t2
)
{
    return static_cast<const Form1&>(t1) & static_cast<const Form2&>(t2);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Form, class Cmpt>
inline Form operator-
(
    const ArraySpace<Form, Cmpt>& as
)
{
    Form v(as.nCmpt_);
    for (int i = 0; i < as.nCmpt_; ++i)
    {
      v.v_[i] = -as.v_[i];
    }
    return v;
}


template<class Form, class Cmpt>
inline Form operator+
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
#   ifdef FULLDEBUG
    if (as1.nCmpt_ != as2.nCmpt_)
    {
        FatalErrorIn("operator+(ArraySpace<Form, Cmpt>&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    Form v(as1.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = as1.v_[i] + as2.v_[i];
    }
    return v;
}


template<class Form, class Cmpt>
inline Form operator-
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
#   ifdef FULLDEBUG
    if (as1.nCmpt_ != as2.nCmpt_)
    {
        FatalErrorIn("operator-(ArraySpace<Form, Cmpt>&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    Form v(as1.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = as1.v_[i] - as2.v_[i];
    }
    return v;
}


template<class Form, class Cmpt>
inline Form operator*
(
    scalar s,
    const ArraySpace<Form, Cmpt>& as
)
{
    Form v(as.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = s * as.v_[i];
    }
    return v;
}


template<class Form, class Cmpt>
inline Form operator*
(
    const ArraySpace<Form, Cmpt>& as,
    scalar s
)
{
    Form v(as.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = as.v_[i] * s;
    }
    return v;
}


template<class Form, class Cmpt>
inline Form operator/
(
    const ArraySpace<Form, Cmpt>& as,
    scalar s
)
{
    Form v(as.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = as.v_[i] / s;
    }
    return v;
}


template<class Form, class Cmpt>
inline Form operator/
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
#   ifdef FULLDEBUG
    if (as1.nCmpt_ != as2.nCmpt_)
    {
        FatalErrorIn("operator/(ArraySpace<Form, Cmpt>&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    Form v(as1.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = as1.v_[i] / as2.v_[i];
    }
    return v;
}


template<class Form, class Cmpt>
inline Form operator/
(
    scalar s,
    const ArraySpace<Form, Cmpt>& as
)
{
    Form v(as.nCmpt_);
    for (int i = 0; i < v.nCmpt_; ++i)
    {
      v.v_[i] = s / as.v_[i];
    }
    return v;
}


template<class Form, class Cmpt>
inline Cmpt operator&&
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
#   ifdef FULLDEBUG
    if (as1.nCmpt_ != as2.nCmpt_)
    {
        FatalErrorIn("operator&&(ArraySpace<Form, Cmpt>&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    Cmpt ddProd = as1.v_[0] * as2.v_[0];
    for (int i = 1; i < as1.nCmpt_; ++i)
    {
        ddProd += as1.v_[i] * as2.v_[i];
    }
    return ddProd;
}


template<class Form, class Cmpt>
inline bool operator==
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
    bool eq = (as1.nCmpt_ == as2.nCmpt_);
    for (int i = 0; i < as1.nCmpt_; ++i)
    {
        if (!(eq &= (equal(as1.v_[i], as2.v_[i])))) break;
    }
    return eq;
}


template<class Form, class Cmpt>
inline bool operator!=
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
    return !(as1 == as2);
}


template<class Form, class Cmpt>
inline bool operator>
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
#   ifdef FULLDEBUG
    if (as1.nCmpt_ != as2.nCmpt_)
    {
        FatalErrorIn("operator>(ArraySpace<Form, Cmpt>&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    bool gt = true;
    for (int i = 0; i < as1.nCmpt_; ++i)
    {
        if (!(gt &= as1.v_[i] > as2.v_[i])) break;
    }
    return gt;
}


template<class Form, class Cmpt>
inline bool operator<
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
#   ifdef FULLDEBUG
    if (as1.nCmpt_ != as2.nCmpt_)
    {
        FatalErrorIn("operator<(ArraySpace<Form, Cmpt>&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    bool lt = true;
    for (int i = 0; i < as1.nCmpt_; ++i)
    {
        if (!(lt &= as1.v_[i] < as2.v_[i])) break;
    }
    return lt;
}


template<class Form, class Cmpt>
inline bool operator>=
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
#   ifdef FULLDEBUG
    if (as1.nCmpt_ != as2.nCmpt_)
    {
        FatalErrorIn("operator>=(ArraySpace<Form, Cmpt>&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    return !(as1 < as2);
}


template<class Form, class Cmpt>
inline bool operator<=
(
    const ArraySpace<Form, Cmpt>& as1,
    const ArraySpace<Form, Cmpt>& as2
)
{
#   ifdef FULLDEBUG
    if (as1.nCmpt_ != as2.nCmpt_)
    {
        FatalErrorIn("operator<=(ArraySpace<Form, Cmpt>&, ArraySpace<Form, Cmpt>&)")
            << "ranges not match"
            << abort(FatalError);
    }
#   endif

    return !(as1 > as2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
