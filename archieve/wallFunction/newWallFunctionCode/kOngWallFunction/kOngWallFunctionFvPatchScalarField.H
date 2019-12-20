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

Class
    Foam::kOngWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on velocity.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            kOngWallFunction;
    }
    \endverbatim

See also
    Foam::kLowReWallFunctionFvPatchScalarField

SourceFiles
    kOngWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef kOngWallFunctionFvPatchScalarField_H
#define kOngWallFunctionFvPatchScalarField_H


#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class kOngWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class kOngWallFunctionFvPatchScalarField
:
    //public kLowReWallFunctionFvPatchScalarField
public fixedValueFvPatchField<scalar>
    
{ 
protected:
 
//- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;

        //- Ceps2 coefficient
        scalar Ceps2_;

        
        //- Check the type of the patch
        virtual void checkType();
       // Protected Member Functions                
        virtual void updateCoeffs();
        //- Calculate yPLus
       virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;
       virtual tmp<scalarField> calcUTau2(const scalarField& magGradU) const;


public:

    //- Runtime type information
    TypeName("kOngWallFunction");


    // Constructors

        //- Construct from patch and internal field
        kOngWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        kOngWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  kOngWallFunctionFvPatchScalarField
        //  onto a new patch
        kOngWallFunctionFvPatchScalarField
        (
            const kOngWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        kOngWallFunctionFvPatchScalarField
        (
            const kOngWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new kOngWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        kOngWallFunctionFvPatchScalarField
        (
            const kOngWallFunctionFvPatchScalarField&,
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
                new kOngWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions
//- Evaluate the patchField
            virtual void evaluate(const Pstream::commsTypes);
              //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
