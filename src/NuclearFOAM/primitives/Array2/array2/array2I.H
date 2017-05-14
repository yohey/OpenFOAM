
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const array2 array2::zero(1, 1, 0);

template<>
const array2 array2::one(1, 1, 1);

template<>
const array2 array2::max(1, 1, VGREAT);

template<>
const array2 array2::min(1, 1, -VGREAT);

template<>
const array2 array2::rootMax(1, 1, ROOTVGREAT);

template<>
const array2 array2::rootMin(1, 1, -ROOTVGREAT);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
inline const array2 array2::zeros(int nRows, int nCols) { return array2(nRows, nCols, 0); };

template<>
inline const array2 array2::ones(int nRows, int nCols) { return array2(nRows, nCols, 1); };

// template<>
// inline const array2 array2::max(int nRows, int nCols) { return array2(nRows, nCols, VGREAT); };

// template<>
// inline const array2 array2::min(int nRows, int nCols) { return array2(nRows, nCols, -VGREAT); };

// template<>
// inline const array2 array2::rootMax(int nRows, int nCols) { return array2(nRows, nCols, ROOTVGREAT); };

// template<>
// inline const array2 array2::rootMin(int nRows, int nCols) { return array2(nRows, nCols, -ROOTVGREAT); };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
