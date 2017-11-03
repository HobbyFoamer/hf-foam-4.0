/*---------------------------------------------------------------------------*\

Class
    Foam::VONOSLimiter

Description
    Class with limiter function which returns the limiter for the
    VONOS differencing scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    VONOS.C

\*---------------------------------------------------------------------------*/

#ifndef VONOS_H
#define VONOS_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class VONOSLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class VONOSLimiter
:
    public LimiterFunc
{

public:

    VONOSLimiter(Istream&)
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

        return max(min(max(min(18.0*r, (3.0 + r)/4.0), r), 2.0), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
