/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "LimitedScheme.H"
#include "Limited01.H"
#include "mSMART.H"

#include "DeferredCorrectionLimitedScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLimitedSurfaceInterpolationScheme(mSMART, mSMARTLimiter)
    makeLimitedVSurfaceInterpolationScheme(mSMARTV, mSMARTLimiter)

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        mSMART01,
        Limited01Limiter,
        mSMARTLimiter,
        NVDTVD,
        magSqr,
        scalar
    )


    // Deferred correction schemes
    makeDeferredSurfaceInterpolationScheme(mSMARTDC, mSMARTLimiter)
    makeDeferredVSurfaceInterpolationScheme(mSMARTVDC, mSMARTLimiter)

    makeLDeferredSurfaceInterpolationTypeScheme
    (
        mSMART01DC,
        Limited01Limiter,
        mSMARTLimiter,
        NVDTVD,
        magSqr,
        scalar
    )
}

// ************************************************************************* //
