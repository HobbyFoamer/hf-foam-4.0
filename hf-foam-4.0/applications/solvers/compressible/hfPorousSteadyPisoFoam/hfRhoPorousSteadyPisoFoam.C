/*---------------------------------------------------------------------------*\

Application
    myRhoPorousSteadyPisoFoam

Description
    Steady-state solver for turbulent flow of compressible fluids with
    RANS turbulence modelling, and implicit or explicit porosity treatment

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "basicPsiThermo.H"
#include "RASModel.H"
#include "porousZones.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
			
            int nCorrMax = 2;
            for(int nCorr = 0; nCorr < nCorrMax; nCorr++)
            {
                #include "pEqn.H"
            }
			
            #include "hEqn.H"
			
			UEqn.clear();
        }

        turbulence->correct();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
