/*---------------------------------------------------------------------------*\

Class
    Foam::mVONOSLimiter

Description
    Class with limiter function which returns the limiter for the
    modified VONOS differencing scheme based on r obtained from the
    LimiterFunc class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    mVONOS.C

\*---------------------------------------------------------------------------*/

#ifndef mVONOS_H
#define mVONOS_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class mVONOSLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class mVONOSLimiter
:
    public LimiterFunc
{

public:

    mVONOSLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        scalar r = LimiterFunc::r
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );

        return max(min(max(min(18.0*r, (3.0 + r)/4.0), r), 1.8), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
