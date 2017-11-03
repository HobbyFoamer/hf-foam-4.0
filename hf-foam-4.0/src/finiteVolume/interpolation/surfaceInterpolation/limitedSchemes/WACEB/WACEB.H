/*---------------------------------------------------------------------------*\

Class
    Foam::WACEBLimiter

Description
    Class with limiter function which returns the limiter for the
    WACEB differencing scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    WACEB.C

\*---------------------------------------------------------------------------*/

#ifndef WACEB_H
#define WACEB_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class WACEBLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class WACEBLimiter
:
    public LimiterFunc
{

public:

    WACEBLimiter(Istream&)
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

        return max(min(min(2.0*r, (3.0 + r)/4.0), 2.0), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
