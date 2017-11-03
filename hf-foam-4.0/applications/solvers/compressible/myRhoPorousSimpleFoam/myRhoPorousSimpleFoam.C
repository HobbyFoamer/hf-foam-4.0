/*---------------------------------------------------------------------------*\

Application
    myRhoPorousSimpleFoam

Description
    Steady-state solver for turbulent flow of compressible fluids with
    RANS turbulence modelling, and implicit or explicit porosity treatment

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "RASModel.H"
#include "porousZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readSIMPLEControls.H"
        #include "initConvergenceCheck.H"

        p.storePrevIter();
        rho.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "hEqn.H"
            #include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        #include "convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
