/*---------------------------------------------------------------------------*\

\*---------------------------------------------------------------------------*/

#include "SpalartAllmarasHF.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDdt.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SpalartAllmarasHF, 0);
addToRunTimeSelectionTable(RASModel, SpalartAllmarasHF, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmarasHF::chi() const
{
    return rho_*nuTilda_/mu();
}


tmp<volScalarField> SpalartAllmarasHF::fv1(const volScalarField& chi) const
{
    volScalarField chi3 = pow3(chi);
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> SpalartAllmarasHF::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    // Model adjustment, Eric Paterson: boundedness of
    // turbulent viscosity.  HJ, 24/Aug/2010
    return 1.0 - chi/(1.0 + chi*fv1);
    //return 1.0/pow3(scalar(1) + chi/Cv2);
}


tmp<volScalarField> SpalartAllmarasHF::fv3
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    volScalarField chiByCv2 = (1/Cv2_)*chi;

    return
        (scalar(1) + chi*fv1)
       *(1/Cv2_)
       *(3*(scalar(1) + chiByCv2) + sqr(chiByCv2))
       /pow3(scalar(1) + chiByCv2);
}


tmp<volScalarField> SpalartAllmarasHF::fw(const volScalarField& Stilda) const
{
    volScalarField r = min
    (
        mag(nuTilda_
       /(
           max(Stilda, dimensionedScalar("SMALL", Stilda.dimensions(), SMALL))
           *sqr(kappa_*d_)
        )),
        scalar(10.0)
    );
    r.boundaryField() == 0.0;

    volScalarField g = r + Cw2_*(pow6(r) - r);

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}

tmp<volScalarField> SpalartAllmarasHF::fr1(const volSymmTensorField& S, const volTensorField& W) const
{
    if (spalartShurCorrection_)
    {
        volScalarField sqrS = 2.0*magSqr(S);
        volScalarField sqrW = 2.0*magSqr(W);
        volScalarField sqrD = 0.5*(sqrS + sqrW);

        volScalarField rStar = sqrt(sqrS/max(sqrW, dimensionedScalar("SMALL", sqrW.dimensions(), SMALL)));

        volScalarField Wxx = W.component(tensor::XX);
        volScalarField Wxy = W.component(tensor::XY);
        volScalarField Wxz = W.component(tensor::XZ);
        volScalarField Wyx = W.component(tensor::YX);
        volScalarField Wyy = W.component(tensor::YY);
        volScalarField Wyz = W.component(tensor::YZ);
        volScalarField Wzx = W.component(tensor::ZX);
        volScalarField Wzy = W.component(tensor::ZY);
        volScalarField Wzz = W.component(tensor::ZZ);
        
        volScalarField rTilda =
            2.0/(rho_*sqr(max(sqrD, dimensionedScalar("SMALL", sqrD.dimensions(), SMALL))))
               *(
                   (Wxx*Sxx_ + Wxy*Sxy_ + Wxz*Sxz_)
                      *(fvc::ddt(rho_, Sxx_) + fvc::div(phi_, Sxx_)) //i,j=1,1
                 + (Wxx*Syx_ + Wxy*Syy_ + Wxz*Syz_)
                      *(fvc::ddt(rho_, Sxy_) + fvc::div(phi_, Sxy_)) //i,j=1,2
                 + (Wxx*Szx_ + Wxy*Szy_ + Wxz*Szz_)
                      *(fvc::ddt(rho_, Sxz_) + fvc::div(phi_, Sxz_)) //i,j=1,3
                 + (Wyx*Sxx_ + Wyy*Sxy_ + Wyz*Sxz_)
                      *(fvc::ddt(rho_, Syx_) + fvc::div(phi_, Syx_)) //i,j=2,1
                 + (Wyx*Syx_ + Wyy*Syy_ + Wyz*Syz_)
                      *(fvc::ddt(rho_, Syy_) + fvc::div(phi_, Syy_)) //i,j=2,2
                 + (Wyx*Szx_ + Wyy*Szy_ + Wyz*Szz_)
                      *(fvc::ddt(rho_, Syz_) + fvc::div(phi_, Syz_)) //i,j=2,3
                 + (Wzx*Sxx_ + Wzy*Sxy_ + Wzz*Sxz_)
                      *(fvc::ddt(rho_, Szx_) + fvc::div(phi_, Szx_)) //i,j=3,1
                 + (Wzx*Syx_ + Wzy*Syy_ + Wzz*Syz_)
                      *(fvc::ddt(rho_, Szy_) + fvc::div(phi_, Szy_)) //i,j=3,2
                 + (Wzx*Szx_ + Wzy*Szy_ + Wzz*Szz_)
                      *(fvc::ddt(rho_, Szz_) + fvc::div(phi_, Szz_)) //i,j=3,3
                );

        return 
            (1 + Cr1_)*2.0*rStar/max((1.0 + rStar),SMALL)*(1.0 - Cr3_*atan(Cr2_*rTilda))
             - Cr1_;
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "fr1",
                    mesh_.time().timeName(),
                    U_.db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("fr1", dimless, 1.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}

/*
 * calculation of ft2
 * returns 0.0 when negativeNuTilda = false
 */
tmp<volScalarField> SpalartAllmarasHF::ft2(const volScalarField& chi) const
{
    if(negativeNuTilda_)
    {
        return Ct3_*exp(-Ct4_*sqr(chi));
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "ft2",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("ft2", dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

SpalartAllmarasHF::SpalartAllmarasHF
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const basicThermo& thermophysicalModel
)
:
    RASModel(typeName, rho, U, phi, thermophysicalModel),

    sigmaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Prt_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Prt",
            coeffDict_,
            1.0
        )
    ),

    Cb1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            coeffDict_,
            7.1
        )
    ),
    Cv2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv2",
            coeffDict_,
            5.0
        )
    ),
    Cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr1",
            coeffDict_,
            1.0
        )
    ),
    Cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr2",
            coeffDict_,
            12.0
        )
    ),
    Cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr3",
            coeffDict_,
            1.0
        )
    ),
    C5_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5",
            coeffDict_,
            3.5
        )
    ),
    Cn1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn1",
            coeffDict_,
            16.0
        )
    ),
    Cn2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn2",
            coeffDict_,
            0.7
        )
    ),
    Cn3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cn3",
            coeffDict_,
            0.9
        )
    ),
    Ct3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct3",
            coeffDict_,
            1.2
        )
    ),
    Ct4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ct4",
            coeffDict_,
            0.5
        )
    ),

    compressibilityCorrection_
    (
        coeffDict_.lookupOrDefault("compressibilityCorrection", false)
    ),

    spalartShurCorrection_
    (
        coeffDict_.lookupOrDefault("spalartShurCorrection", false)
    ),

    StildaModification_
    (
        coeffDict_.lookupOrDefault("StildaModification", false)
    ),

    negativeNuTilda_
    (
        coeffDict_.lookupOrDefault("negativeNuTilda", false)
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    mut_
    (
        IOobject
        (
            "mut",
            runTime_.timeName(),
            U_.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    alphat_
    (
        IOobject
        (
            "alphat",
            runTime_.timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateAlphat("alphat", mesh_)
    ),

    Sxx_
    (
        IOobject
        (
            "Sxx",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Sxx", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Sxy_
    (
        IOobject
        (
            "Sxy",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Sxy", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Sxz_
    (
        IOobject
        (
            "Sxz",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Sxz", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Syx_
    (
        IOobject
        (
            "Syx",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Syx", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Syy_
    (
        IOobject
        (
            "Syy",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Syy", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Syz_
    (
        IOobject
        (
            "Syz",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Syz", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Szx_
    (
        IOobject
        (
            "Szx",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Szx", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Szy_
    (
        IOobject
        (
            "Szy",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Szy", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    Szz_
    (
        IOobject
        (
            "Szz",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Szz", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    d_(mesh_)
{
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();

    printCoeffs();

    if (compressibilityCorrection_)
    {
        Info<< "    Employing compressibility correction" << endl;
    }
    if (spalartShurCorrection_)
    {
        Info<< "    Employing Rotation/Curvature correction" << endl;
    }
    if (StildaModification_)
    {
        Info<< "    Enabling new Stilda modification" << endl;
    }
    if (negativeNuTilda_)
    {
        Info<< "    Enabling negative nuTilda" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> SpalartAllmarasHF::DnuTildaEff(const volScalarField& chi) const
{
    volScalarField pow3chi = pow(chi, 3);
    volScalarField fn =
        pos(chi) + neg(chi)*(Cn1_ + pow3chi)/(Cn1_ - pow3chi);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "DnuTildaEff", 
            ((rho_*fn*nuTilda_ + mu())/sigmaNut_)
        )
    );
}


tmp<volSymmTensorField> SpalartAllmarasHF::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k() - (mut_/rho_)*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<volSymmTensorField> SpalartAllmarasHF::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                U_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -muEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> SpalartAllmarasHF::divDevRhoReff() const
{
    volScalarField muEff_ = muEff();

    return
    (
      - fvm::laplacian(muEff_, U_)
      - fvc::div(muEff_*dev2(T(fvc::grad(U_))))
    );
}


bool SpalartAllmarasHF::read()
{
    if (RASModel::read())
    {
        sigmaNut_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Prt_.readIfPresent(coeffDict());

        Cb1_.readIfPresent(coeffDict());
        Cb2_.readIfPresent(coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        Cv1_.readIfPresent(coeffDict());
        Cv2_.readIfPresent(coeffDict());
        Cr1_.readIfPresent(coeffDict());
        Cr2_.readIfPresent(coeffDict());
        Cr3_.readIfPresent(coeffDict());
        C5_.readIfPresent(coeffDict());
        Cn1_.readIfPresent(coeffDict());
        Cn2_.readIfPresent(coeffDict());
        Cn3_.readIfPresent(coeffDict());
        Ct3_.readIfPresent(coeffDict());
        Ct4_.readIfPresent(coeffDict());

        compressibilityCorrection_.readIfPresent
        (
            "compressibilityCorrection", coeffDict()
        );
        spalartShurCorrection_.readIfPresent
        (
            "spalartShurCorrection", coeffDict()
        );
        StildaModification_.readIfPresent
        (
            "StildaModification", coeffDict()
        );
        negativeNuTilda_.readIfPresent
        (
            "negativeNuTilda", coeffDict()
        );

        return true;
    }
    else
    {
        return false;
    }
}


void SpalartAllmarasHF::correct()
{
    // Bound in case of topological change
    // HJ, 22/Aug/2007
    if (mesh_.changing())
    {
        bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    }

    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*nuTilda_*fv1(chi());
        mut_.correctBoundaryConditions();

        // Re-calculate thermal diffusivity
        alphat_ = mut_/Prt_;
        alphat_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    if (mesh_.changing())
    {
        d_.correct();
    }

    volScalarField chi = this->chi();
    volScalarField fv1 = this->fv1(chi);

    volTensorField gradU = fvc::grad(U_);
    volSymmTensorField S = symm(gradU);
    volTensorField W = skew(gradU);
    volScalarField Omega = sqrt(2.0)*mag(W);

    volScalarField Stilda =
//        fv3(chi, fv1)*::sqrt(2.0)*mag(skew(fvc::grad(U_)))
        sqrt(2.0)*mag(W)
      + fv2(chi, fv1)*nuTilda_/sqr(kappa_*d_);

    // new Stilda modification used when StildaModification = true
    if (StildaModification_)
    {
        volScalarField Sbar
        (
            fv2(chi, fv1)*nuTilda_/sqr(kappa_*d_)
        );
        Stilda = Omega
               + pos(Cn2_*Omega + Sbar)*Sbar
               + neg(Cn2_*Omega + Sbar)
               *(Omega*(sqr(Cn2_)*Omega + Cn3_*Sbar))
               /max(((Cn3_ - 2.0*Cn2_)*Omega - Sbar),
                    dimensionedScalar("SMALL", Omega.dimensions(), SMALL));
    }
    
    Sxx_ = S.component(tensor::XX);
    Sxy_ = S.component(tensor::XY);
    Sxz_ = S.component(tensor::XZ);
    Syx_ = S.component(tensor::YX);
    Syy_ = S.component(tensor::YY);
    Syz_ = S.component(tensor::YZ);
    Szx_ = S.component(tensor::ZX);
    Szy_ = S.component(tensor::ZY);
    Szz_ = S.component(tensor::ZZ);
    
    // calculation of compressibility correction
    volScalarField compCorr
    (
        IOobject
        (
            "compCorr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("compCorr", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
    );
    if(compressibilityCorrection_)
    {
        volScalarField Cp = thermo().Cp();
        volScalarField Cv = thermo().Cv();
        volScalarField a = sqrt((Cp/Cv)*(Cp - Cv)*thermo().T());
        compCorr = -C5_*rho_*sqr(nuTilda_/a)*magSqr(gradU);
    }

    volScalarField fr1 = this->fr1(S, W);
    volScalarField ft2 = this->ft2(chi);

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(rho_, nuTilda_)
      + fvm::div(phi_, nuTilda_)
      - fvm::laplacian(DnuTildaEff(chi), nuTilda_)
      - Cb2_/sigmaNut_*rho_*magSqr(fvc::grad(nuTilda_))
     ==
        pos(nuTilda_)
       *(
         // RHS for nuTilda >= 0, rotation/curvature correction factor is multiplied
            Cb1_*(1.0 - ft2)*fr1*Stilda*rho_*nuTilda_
          - fvm::Sp((Cw1_*fw(Stilda) - Cb1_/sqr(kappa_)*ft2)*rho_*nuTilda_/sqr(d_), nuTilda_)
        )
      + neg(nuTilda_)
       *(
         // RHS for nuTilda < 0
            Cb1_*(1.0 - Ct3_)*Omega*rho_*nuTilda_
          + fvm::Sp(Cw1_*rho_*nuTilda_/sqr(d_), nuTilda_)
        )
    );

    nuTildaEqn().relax();
    solve(nuTildaEqn);
    if(!negativeNuTilda_)
    {
        // bound nuTilda when negativeNuTilda = false
        bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    }
    //nuTilda_.correctBoundaryConditions();
    chi = this->chi();
    fv1 = this->fv1(chi);

    // Re-calculate viscosity
    mut_.internalField() = mag(fv1)*nuTilda_.internalField()*rho_.internalField();
    if(negativeNuTilda_)
    {
        // bound mut when negativeNuTilda = true
        bound(mut_, dimensionedScalar("0", mut_.dimensions(), 0.0));
    }
    mut_.correctBoundaryConditions();
    nuTilda_.correctBoundaryConditions();

    // Re-calculate thermal diffusivity
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
