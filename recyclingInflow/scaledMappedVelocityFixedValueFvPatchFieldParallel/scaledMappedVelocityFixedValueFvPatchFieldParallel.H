/*---------------------------------------------------------------------------*\
Class
    Foam::scaledMappedVelocityFixedValueFvPatchField

Group
    grpInletBoundaryConditions grpCoupledBoundaryConditions

Description
    This boundary condition maps the velocity and flux from a neighbour patch
    to this patch.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            scaledMappedVelocity;
        phi             phi;
        value           uniform 0;  // place holder
    }
    \endverbatim

    The underlying sample mode should be set to \c nearestPatchFace or
    \c nearestFace

Note
    This boundary condition can only be applied to patches that are of
    the \c mappedPolyPatch type.

SeeAlso
    Foam::mappedPatchBase
    Foam::mappedPolyPatch
    Foam::mappedFvPatch
    Foam::fixedValueFvPatchVectorField

SourceFiles
    scaledMappedVelocityFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef scaledMappedVelocityFixedValueFvPatchFieldParallel_H
#define scaledMappedVelocityFixedValueFvPatchFieldParallel_H

#include "fixedValueFvPatchFields.H"
#include "mappedFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class scaledMappedVelocityFixedValueFvPatch Declaration
\*---------------------------------------------------------------------------*/

class scaledMappedVelocityFixedValueFvPatchField
:
    public fixedValueFvPatchVectorField
{
    // Private data
        //scalar deltaInlet_; // boundary layer displacement thickness at the inlet
        //scalar thetaInlet_; // boundary layer momentumt thickness at the inlet
        //scalar nu_;         // kinematic viscosity
        //scalar Ue_;         // freestream velocity
        //scalar t_;          // initial averaging time
        //vectorField UMeanSpanTime_; // Averaged velocity at the recycled plane

public:

    //- Runtime type information
    TypeName("scaledMappedVelocity");

    // Constructors

        //- Construct from patch and internal field
        scaledMappedVelocityFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        scaledMappedVelocityFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  scaledMappedVelocityFixedValueFvPatchField
        //  onto a new patch
        scaledMappedVelocityFixedValueFvPatchField
        (
            const scaledMappedVelocityFixedValueFvPatchField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        scaledMappedVelocityFixedValueFvPatchField
        (
            const scaledMappedVelocityFixedValueFvPatchField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new scaledMappedVelocityFixedValueFvPatchField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        scaledMappedVelocityFixedValueFvPatchField
        (
            const scaledMappedVelocityFixedValueFvPatchField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new scaledMappedVelocityFixedValueFvPatchField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 

#endif 

// ************************************************************************* //
