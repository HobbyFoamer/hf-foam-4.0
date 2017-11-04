/*---------------------------------------------------------------------------*\

Class
    Foam::compressible::RASModels::SpalartAllmarasHF

Description
    Negative Spalart-Allmaras one-eqn mixing-length model for compressible
    external flows.
    Rotation/Curvature correction, Compressibility correction and modified
    expression of Stilda are included.    
    
    \verbatim
        "Turbulence Modeling in Rotating and Curved Channels: Assessing the
        Spalart-Shur Correction"
        M.L. Shur,
        M.K. Strelets,
        A.K. Travin,
        P.R. Spalart,
        AIAA Journal Vol. 38, No. 5, 2000, pp. 784-792.
    \endverbatim
    using the optional flag \c spalartShurCorrection,

    \verbatim
        "Trends in Turbulence Treatments"
        P.R. Spalart,
        AIAA Paper No. 2000-2603, 2000.
    \endverbatim
    using the optional flag \c compressibilityCorrection,

    \verbatim
        "Modifications and Clarifications for the Implementation of the
        Spalart-Allmaras Turbulence Model"
        S.R. Allmaras,
        F.T. Johnson,
        P.R. Spalart,
        ICCFD7-1902, 7th International Conference on
        Computational Fluid Dynamics, 2012.
    \endverbatim
    using the optional flag \c StildaModification,
    and \c negativeNuTilda

    The default model coefficients correspond to the following:
    @verbatim
        SpalartAllmarasHFCoeffs
        {
            Cb1                         0.1355;
            Cb2                         0.622;
            Cw2                         0.3;
            Cw3                         2.0;
            Cv1                         7.1;
            Cv2                         5.0;
            sigmaNut                    0.66666;
            Prt                         1.0;    // only for compressible
            kappa                       0.41;
            
            // coefficients for rotation/curvature correction
            Cr1                         1.0;
            Cr2                         12.0;
            Cr3                         1.0;
            
            // coefficients for compressibility correction
            C5                          3.5;
            
            // coefficients for negative SA and new Stilda expression
            Cn1                         16.0;
            Cn2                         0.7;
            Cn3                         0.9;
            Ct3                         1.2;
            Ct4                         0.5;
            
            // switches
            compressibilityCorrection   false;
            spalartShurCorrection       false;
            StildaModification          false;
            negativeNutilda             false;
        }
    @endverbatim

SourceFiles
    SpalartAllmarasHF.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleSpalartAllmarasHF_H
#define compressibleSpalartAllmarasHF_H

#include "RASModel.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class SpalartAllmarasHF Declaration
\*---------------------------------------------------------------------------*/

class SpalartAllmarasHF
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;
            dimensionedScalar Prt_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cv2_;
            dimensionedScalar Cr1_;
            dimensionedScalar Cr2_;
            dimensionedScalar Cr3_;
            dimensionedScalar C5_;
            dimensionedScalar Cn1_;
            dimensionedScalar Cn2_;
            dimensionedScalar Cn3_;
            dimensionedScalar Ct3_;
            dimensionedScalar Ct4_;


        //- Optional flag to activate the compressibility correction
        Switch compressibilityCorrection_;

        //- Optional flag to activate the Rotation/Curvature correction
        Switch spalartShurCorrection_;

        //- Optional flag to activate the new STilda correction
        Switch StildaModification_;

        //- Optional flag to activate the negative nuTilda model
        Switch negativeNuTilda_;


        // Fields

            volScalarField nuTilda_;
            volScalarField mut_;
            volScalarField alphat_;
            
            volScalarField Sxx_;
            volScalarField Sxy_;
            volScalarField Sxz_;
            volScalarField Syx_;
            volScalarField Syy_;
            volScalarField Syz_;
            volScalarField Szx_;
            volScalarField Szy_;
            volScalarField Szz_;


        //- Wall distance
        wallDist d_;


    // Private member functions

        tmp<volScalarField> chi() const;
        tmp<volScalarField> fv1(const volScalarField& chi) const;
        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fv3
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;
        tmp<volScalarField> fw(const volScalarField& Stilda) const;

        tmp<volScalarField> fr1(const volSymmTensorField& S, const volTensorField& W) const;

        tmp<volScalarField> ft2(const volScalarField& chi) const;


public:

    //- Runtime type information
    TypeName("SpalartAllmarasHF");


    // Constructors

        //- Construct from components
        SpalartAllmarasHF
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const basicThermo& thermophysicalModel
        );


    //- Destructor
    virtual ~SpalartAllmarasHF()
    {}


    // Member Functions

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff(const volScalarField& chi) const;
/*
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DnuTildaEff",
                    rho_*nuTilda_/sigmaNut_ + mu()
                )
            );
        }
*/
        //- Return the turbulence viscosity
        virtual tmp<volScalarField> mut() const
        {
            return mut_;
        }

        //- Return the effective turbulent thermal diffusivity
        virtual tmp<volScalarField> alphaEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("alphaEff", alphat_ + alpha())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "k",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 2, -2, 0, 0), 0)
                )
            );
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "epslion",
                        runTime_.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0), 0)
                )
            );
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devRhoReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
