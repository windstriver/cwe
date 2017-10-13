/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kOngWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void kOngWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}



// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void kOngWallFunctionFvPatchScalarField::updateCoeffs() 
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];


    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar Cmu25 = pow025(Cmu_);

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));  
    const scalarField magGradU(mag(Uw.snGrad()));
    scalarField& kw = *this;

    // Set k wall values
    forAll(kw, facei)
    {
       label faceCelli = patch().faceCells()[facei];

        scalar uTau = Cmu25*sqrt(k[faceCelli]);

        scalar yPlus = uTau*y[facei]/nuw[facei]; 
       tmp<scalarField> tuTau = calcUTau(magGradU);
    scalarField& uts = tuTau.ref();

     tmp<scalarField> tuTau2 = calcUTau2(magGradU);
    scalarField& uts2 = tuTau2.ref();

         if (yPlus <= 5)
        {
            
            kw[facei] = magUp[facei]*nuw[facei]/mag(y[facei]);
        } 
        else if (yPlus > 5 && yPlus < 30)
      {
       
       kw[facei] = uts2[facei]*uts2[facei];

        }
        else
        {
       
       kw[facei] = uts[facei]*uts[facei];

        }

        kw[facei] /= sqrt(Cmu_);
    }

    // Limit kw to avoid failure of the turbulence model due to division by kw
    kw = max(kw, SMALL);

  fixedValueFvPatchField<scalar>::updateCoeffs();
//kLowReWallFunctionFvPatchScalarField ::updateCoeffs();


    // TODO: perform averaging for cells sharing more than one boundary face
}


tmp<scalarField> kOngWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau.ref();

    forAll(uTau, facei)
    {
        scalar ut = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = max(kappa_*magUp[facei]/ut, 13.86);
                scalar fkUu = exp(kUu);

                scalar f =
                    - ut*y[facei]/nuw[facei]+ 1/E_*fkUu;

                scalar df =
                    y[facei]/nuw[facei] + 1/E_*kUu*fkUu/ut;

                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > ROOTVSMALL && err > 0.01 && ++iter < 10);

            uTau[facei] = max(0.0, ut);
        }
    }

    return tuTau;
}

tmp<scalarField> kOngWallFunctionFvPatchScalarField::calcUTau2
(
    const scalarField& magGradU
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    tmp<scalarField> tuTau2(new scalarField(patch().size(), 0.0));
    scalarField& uTau2 = tuTau2.ref();

    forAll(uTau2, facei)
    {
        scalar ut2 = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);

        if (ut2 > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                 
            scalar lg = log(E_*y[facei]*ut2/nuw[facei]);
            scalar yp = min(y[facei]*ut2/nuw[facei], 30);
                   yp = max(5,yp);
                  
                scalar f =
                    -ut2/magUp[facei]+ (kappa_*(yp-5)*yp-(30-yp)*lg)/(25*yp*lg);

                scalar df1=
                   ((-2*kappa_*yp+5*kappa_-lg)*y[facei]/nuw[facei]+(30-yp)/ut2)/(25*yp*lg);
               scalar df2=
                   -(-kappa_*(yp-5)*yp+(30-yp)*lg)*(yp/ut2+lg*y[facei]/nuw[facei]) /(25*sqr(yp*lg));
                scalar df =
                  1/magUp[facei]+df1+df2;

                scalar uTauNew = ut2 + f/df;
                err = mag((ut2 - uTauNew)/ut2);
                ut2 = uTauNew;

            } while (ut2 > ROOTVSMALL && err > 0.01 && ++iter < 10);

            uTau2[facei] = max(0.0, ut2);
        }
    }

    return tuTau2;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOngWallFunctionFvPatchScalarField::kOngWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    //kLowReWallFunctionFvPatchScalarField(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    Ceps2_(1.9)
{
    checkType();
}


kOngWallFunctionFvPatchScalarField::kOngWallFunctionFvPatchScalarField
(
    const kOngWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    //kLowReWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    Ceps2_(ptf.Ceps2_)
{
    checkType();
}


kOngWallFunctionFvPatchScalarField::kOngWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    //kLowReWallFunctionFvPatchScalarField(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    Ceps2_(dict.lookupOrDefault<scalar>("Ceps2", 1.9))
{
    checkType();
}


kOngWallFunctionFvPatchScalarField::kOngWallFunctionFvPatchScalarField
(
    const kOngWallFunctionFvPatchScalarField& kwfpsf
)
:
   fixedValueFvPatchField<scalar>(kwfpsf),
    //kLowReWallFunctionFvPatchScalarField(kwfpsf),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_)
{
    checkType();
}


kOngWallFunctionFvPatchScalarField::kOngWallFunctionFvPatchScalarField
(
    const kOngWallFunctionFvPatchScalarField& kwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
   fixedValueFvPatchField<scalar>(kwfpsf, iF),
    //kLowReWallFunctionFvPatchScalarField(kwfpsf, iF),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_)
{
    checkType();
}




// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void kOngWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);

}


void kOngWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ceps2") << Ceps2_ << token::END_STATEMENT << nl;
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    kOngWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
