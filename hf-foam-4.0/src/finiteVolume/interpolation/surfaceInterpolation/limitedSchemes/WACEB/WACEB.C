/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "LimitedScheme.H"
#include "Limited01.H"
#include "WACEB.H"

#include "DeferredCorrectionLimitedScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLimitedSurfaceInterpolationScheme(WACEB, WACEBLimiter)
    makeLimitedVSurfaceInterpolationScheme(WACEBV, WACEBLimiter)

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        WACEB01,
        Limited01Limiter,
        WACEBLimiter,
        NVDTVD,
        magSqr,
        scalar
    )


    // Deferred correction schemes
    makeDeferredSurfaceInterpolationScheme(WACEBDC, WACEBLimiter)
    makeDeferredVSurfaceInterpolationScheme(WACEBVDC, WACEBLimiter)

    makeLDeferredSurfaceInterpolationTypeScheme
    (
        WACEB01DC,
        Limited01Limiter,
        WACEBLimiter,
        NVDTVD,
        magSqr,
        scalar
    )
}

// ************************************************************************* //
