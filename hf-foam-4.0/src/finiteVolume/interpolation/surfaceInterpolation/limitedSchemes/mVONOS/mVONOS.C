/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "LimitedScheme.H"
#include "Limited01.H"
#include "mVONOS.H"

#include "DeferredCorrectionLimitedScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLimitedSurfaceInterpolationScheme(mVONOS, mVONOSLimiter)
    makeLimitedVSurfaceInterpolationScheme(mVONOSV, mVONOSLimiter)

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        mVONOS01,
        Limited01Limiter,
        mVONOSLimiter,
        NVDTVD,
        magSqr,
        scalar
    )


    // Deferred correction schemes
    makeDeferredSurfaceInterpolationScheme(mVONOSDC, mVONOSLimiter)
    makeDeferredVSurfaceInterpolationScheme(mVONOSVDC, mVONOSLimiter)

    makeLDeferredSurfaceInterpolationTypeScheme
    (
        mVONOS01DC,
        Limited01Limiter,
        mVONOSLimiter,
        NVDTVD,
        magSqr,
        scalar
    )
}

// ************************************************************************* //
