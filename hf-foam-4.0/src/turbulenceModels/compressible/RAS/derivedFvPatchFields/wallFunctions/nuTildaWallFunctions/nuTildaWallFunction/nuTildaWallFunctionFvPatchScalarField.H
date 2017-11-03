/*---------------------------------------------------------------------------*\

Class
    Foam::incompressible::RASModels::
    nuTildaWallFunctionFvPatchScalarField

Description
    Wall function boundary condition for walls.
    Calculated nuTilda from mut on the same face.

SourceFiles
    nuTildaWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nuTildaWallFunctionFvPatchScalarField_H
#define nuTildaWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
           Class nuTildaWallFunctionFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nuTildaWallFunctionFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
protected:

    // Protected member functions

        //- Check the type of the patch
        virtual void checkType();

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNuTilda() const;


public:

    //- Runtime type information
    TypeName("compressible::nuTildaWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nuTildaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nuTildaWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nuTildaWallFunctionFvPatchScalarField
        //  onto a new patch
        nuTildaWallFunctionFvPatchScalarField
        (
            const nuTildaWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nuTildaWallFunctionFvPatchScalarField
        (
            const nuTildaWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nuTildaWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nuTildaWallFunctionFvPatchScalarField
        (
            const nuTildaWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nuTildaWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

        //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
