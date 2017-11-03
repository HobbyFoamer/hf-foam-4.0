/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "LimitedScheme.H"
#include "Limited01.H"
#include "VONOS.H"

#include "DeferredCorrectionLimitedScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLimitedSurfaceInterpolationScheme(VONOS, VONOSLimiter)
    makeLimitedVSurfaceInterpolationScheme(VONOSV, VONOSLimiter)

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        VONOS01,
        Limited01Limiter,
        VONOSLimiter,
        NVDTVD,
        magSqr,
        scalar
    )


    // Deferred correction schemes
    makeDeferredSurfaceInterpolationScheme(VONOSDC, VONOSLimiter)
    makeDeferredVSurfaceInterpolationScheme(VONOSVDC, VONOSLimiter)

    makeLDeferredSurfaceInterpolationTypeScheme
    (
        VONOS01DC,
        Limited01Limiter,
        VONOSLimiter,
        NVDTVD,
        magSqr,
        scalar
    )
}

// ************************************************************************* //
