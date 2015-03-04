
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const array array::zero(1, 0);

template<>
const array array::one(1, 1);

template<>
const array array::max(1, VGREAT);

template<>
const array array::min(1, -VGREAT);

template<>
const array array::rootMax(1, ROOTVGREAT);

template<>
const array array::rootMin(1, -ROOTVGREAT);


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
inline const array array::zeros(int nCmpt) { return array(nCmpt, 0); };

template<>
inline const array array::ones(int nCmpt) { return array(nCmpt, 1); };

// template<>
// inline const array array::max(int nCmpt) { return array(nCmpt, VGREAT); };

// template<>
// inline const array array::min(int nCmpt) { return array(nCmpt, -VGREAT); };

// template<>
// inline const array array::rootMax(int nCmpt) { return array(nCmpt, ROOTVGREAT); };

// template<>
// inline const array array::rootMin(int nCmpt) { return array(nCmpt, -ROOTVGREAT); };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
