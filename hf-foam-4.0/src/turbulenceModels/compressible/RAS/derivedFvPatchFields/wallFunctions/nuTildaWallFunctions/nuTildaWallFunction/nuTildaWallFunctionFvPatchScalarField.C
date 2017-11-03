/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "nuTildaWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void nuTildaWallFunctionFvPatchScalarField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorIn("nuTildaWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}

tmp<scalarField> nuTildaWallFunctionFvPatchScalarField::calcNuTilda() const
{
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const scalarField& muw = rasModel.mu().boundaryField()[patch().index()];
    const scalarField& rhow = rasModel.rho().boundaryField()[patch().index()];
    const volScalarField& mut = db().lookupObject<volScalarField>("mut");
    const fvPatchScalarField& mutw = mut.boundaryField()[patch().index()];

    tmp<scalarField> tnuTildaw(new scalarField(patch().size(), 0.0));
    scalarField& nuTildaw = tnuTildaw();

    forAll(mutw, faceI)
    {
        scalar nu = muw[faceI]/rhow[faceI];
        scalar nut = mutw[faceI]/rhow[faceI];
        scalar nT = nut;
        
        if (nut > SMALL)
        {
            int iter = 0;
            scalar nuTildaLast = 0.0;

            do
            {
                nuTildaLast = nT;
                nT = pow(nut*(pow3(nT) + pow3(nu*7.1)), 0.25);
            } while (mag((nT - nuTildaLast)/nT) > 1e-8 && ++iter < 1000);
        
            nuTildaw[faceI] = nT;
        } else {
            nuTildaw[faceI] = 0.0;
        }
    }

    return tnuTildaw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nuTildaWallFunctionFvPatchScalarField::nuTildaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


nuTildaWallFunctionFvPatchScalarField::nuTildaWallFunctionFvPatchScalarField
(
    const nuTildaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


nuTildaWallFunctionFvPatchScalarField::nuTildaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


nuTildaWallFunctionFvPatchScalarField::nuTildaWallFunctionFvPatchScalarField
(
    const nuTildaWallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf)
{}


nuTildaWallFunctionFvPatchScalarField::nuTildaWallFunctionFvPatchScalarField
(
    const nuTildaWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nuTildaWallFunctionFvPatchScalarField::updateCoeffs()
{
    operator==(calcNuTilda());

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nuTildaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
