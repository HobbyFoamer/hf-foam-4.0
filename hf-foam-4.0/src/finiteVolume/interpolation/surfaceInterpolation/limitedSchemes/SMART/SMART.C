/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "LimitedScheme.H"
#include "Limited01.H"
#include "SMART.H"

#include "DeferredCorrectionLimitedScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeLimitedSurfaceInterpolationScheme(SMART, SMARTLimiter)
    makeLimitedVSurfaceInterpolationScheme(SMARTV, SMARTLimiter)

    makeLLimitedSurfaceInterpolationTypeScheme
    (
        SMART01,
        Limited01Limiter,
        SMARTLimiter,
        NVDTVD,
        magSqr,
        scalar
    )


    // Deferred correction schemes
    makeDeferredSurfaceInterpolationScheme(SMARTDC, SMARTLimiter)
    makeDeferredVSurfaceInterpolationScheme(SMARTVDC, SMARTLimiter)

    makeLDeferredSurfaceInterpolationTypeScheme
    (
        SMART01DC,
        Limited01Limiter,
        SMARTLimiter,
        NVDTVD,
        magSqr,
        scalar
    )
}

// ************************************************************************* //
