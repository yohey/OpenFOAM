
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Cmpt>
inline Array2<Cmpt>::Array2(int nRows, int nCols):
    Array<Array<Cmpt>>(nRows)
{
    for (int i = 0; i < this->nCmpt_; ++i)
    {
        new(&this->v_[i]) Array<Cmpt>(nCols);
    }
}


template<class Cmpt>
inline Array2<Cmpt>::Array2(int nRows, int nCols, const Cmpt& v):
    Array<Array<Cmpt>>(nRows)
{
    for (int i = 0; i < this->nCmpt_; ++i)
    {
        new(&this->v_[i]) Array<Cmpt>(nCols, v);
    }
}


template<class Cmpt>
inline Array2<Cmpt>::Array2(const Array2<Cmpt>& as):
    Array<Array<Cmpt>>(as.nCmpt_)
{
    for (int i = 0; i < this->nCmpt_; ++i)
    {
        new(&this->v_[i]) Array<Cmpt>(as[i]);
    }
}


template<class Cmpt>
template<class Cmpt2>
inline Array2<Cmpt>::Array2(const Array2<Cmpt2>& as):
    Array<Array<Cmpt>>(as.nCmpt_)
{
    for (int i = 0; i < this->nCmpt_; ++i)
    {
        new(&this->v_[i]) Array<Cmpt>(as[i]);
    }
}


template<class Cmpt>
inline Array2<Cmpt>::Array2(int nRows, int nCols, Istream& is):
    Array<Array<Cmpt>>(nRows)
{
    for (int i = 0; i < this->nCmpt_; ++i)
    {
        new(&this->v_[i]) Array<Cmpt>(nCols, is);
    }
}


template<class Cmpt>
inline Array2<Cmpt>::Array2(Istream& is):
    Array<Array<Cmpt>>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Cmpt>
inline Array2<Cmpt> operator/
(
    const Array2<Cmpt>& as,
    scalar s
)
{
    Array2<Cmpt> a(as);
    for (int i = 0; i < a.nCmpt_; ++i)
    {
        a.v_[i] = as.v_[i] / s;
    }
    return a;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
