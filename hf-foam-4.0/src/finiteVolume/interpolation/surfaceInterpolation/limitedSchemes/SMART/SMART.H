/*---------------------------------------------------------------------------*\

Class
    Foam::SMARTLimiter

Description
    Class with limiter function which returns the limiter for the
    SMART differencing scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    SMART.C

\*---------------------------------------------------------------------------*/

#ifndef SMART_H
#define SMART_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class SMARTLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class SMARTLimiter
:
    public LimiterFunc
{

public:

    SMARTLimiter(Istream&)
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

        return max(min(min(18.0*r, (3.0 + r)/4.0), 2.0), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
