{
if (pressureImplicitPorosity)
{
    U = trTU()&UEqn().H();
}
else
{
    U = trAU()*UEqn().H();
}

//UEqn.clear();

phi = fvc::interpolate(rho*U) & mesh.Sf();
bool closedVolume = adjustPhi(phi, U, p);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tpEqn = (fvm::laplacian(rho*trTU(), p) == fvc::div(phi));
    }
    else
    {
        tpEqn = (fvm::laplacian(rho*trAU(), p) == fvc::div(phi));
    }

    tpEqn().setReference(pRefCell, pRefValue);
    // retain the residual from the first iteration
    if (nonOrth == 0)
    {
        eqnResidual = tpEqn().solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);
    }
    else
    {
        tpEqn().solve();
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= tpEqn().flux();
    }
}

#include "incompressible/continuityErrs.H"

string rFacs[2] = {"p0", "p1"};

// Explicitly relax pressure for momentum corrector
// except the final loop
if (nCorr != nCorrMax - 1)
{
    p.relax(mesh.solutionDict().relaxationFactor(rFacs[nCorr]));
    p.correctBoundaryConditions();
}

if (pressureImplicitPorosity)
{
    U -= trTU()&fvc::grad(p);
}
else
{
    U -= trAU()*fvc::grad(p);
}

U.correctBoundaryConditions();

// Relax p AFTER U is corrected
// in case of the final loop
if (nCorr == nCorrMax - 1)
{
    p.relax(mesh.solutionDict().relaxationFactor(rFacs[nCorr]));
    p.correctBoundaryConditions();
}

bound(p, pMin);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

if (nCorr == nCorrMax - 1)
{
    rho = thermo.rho();
    //rho = max(rho, rhoMin);
    //rho = min(rho, rhoMax);
    rho.relax();
    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
}
}
