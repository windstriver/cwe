#include "scaledMappedVelocityFixedValueFvPatchField.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scaledMappedVelocityFixedValueFvPatchField::
scaledMappedVelocityFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    deltaInlet_(0.02835),
    thetaInlet_(0.002667),
    nu_(1.4612e-5),
    Ue_(6.355),
    t_(1e-3),
    UMeanSpanTime_(patch().size(), vector::zero)
{}


Foam::scaledMappedVelocityFixedValueFvPatchField::
scaledMappedVelocityFixedValueFvPatchField
(
    const scaledMappedVelocityFixedValueFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    deltaInlet_(ptf.deltaInlet_),
    thetaInlet_(ptf.thetaInlet_),
    nu_(ptf.nu_),
    Ue_(ptf.Ue_),
    t_(ptf.t_),
    UMeanSpanTime_(ptf.UMeanSpanTime_)
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "scaledMappedVelocityFixedValueFvPatchField::"
            "scaledMappedVelocityFixedValueFvPatchField"
            "("
                "const scaledMappedVelocityFixedValueFvPatchField&, "
                "const fvPatch&, "
                "const DimensionedField<vector, volMesh>&, "
                "const fvPatchFieldMapper&"
            ")"
        )   << "Patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
}


Foam::scaledMappedVelocityFixedValueFvPatchField::
scaledMappedVelocityFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict),
    deltaInlet_(readScalar(dict.lookup("deltaInlet"))),
    thetaInlet_(readScalar(dict.lookup("thetaInlet"))),
    nu_(readScalar(dict.lookup("nu"))),
    Ue_(readScalar(dict.lookup("Ue"))),
    t_(readScalar(dict.lookup("t"))),
    UMeanSpanTime_("UMeanSpanTime", dict, patch().size())
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "scaledMappedVelocityFixedValueFvPatchField::"
            "scaledMappedVelocityFixedValueFvPatchField"
            "("
                "const fvPatch&, "
                "const DimensionedField<vector, volMesh>&, "
                "const dictionary&"
            ")"
        )   << "Patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << " for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        this->patch().patch()
    );
    if (mpp.mode() == mappedPolyPatch::NEARESTCELL)
    {
        FatalErrorIn
        (
            "scaledMappedVelocityFixedValueFvPatchField::"
            "scaledMappedVelocityFixedValueFvPatchField"
            "("
                "const fvPatch&, "
                "const DimensionedField<vector, volMesh>&, "
                "const dictionary&"
            ")"
        )   << "Patch " << p.name()
            << " of type '" << p.type()
            << "' can not be used in 'nearestCell' mode"
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }
}


Foam::scaledMappedVelocityFixedValueFvPatchField::
scaledMappedVelocityFixedValueFvPatchField
(
    const scaledMappedVelocityFixedValueFvPatchField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    deltaInlet_(ptf.deltaInlet_),
    thetaInlet_(ptf.thetaInlet_),
    nu_(ptf.nu_),
    Ue_(ptf.Ue_),
    t_(ptf.t_),
    UMeanSpanTime_(ptf.UMeanSpanTime_)
{}


Foam::scaledMappedVelocityFixedValueFvPatchField::
scaledMappedVelocityFixedValueFvPatchField
(
    const scaledMappedVelocityFixedValueFvPatchField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    deltaInlet_(ptf.deltaInlet_),
    thetaInlet_(ptf.thetaInlet_),
    nu_(ptf.nu_),
    Ue_(ptf.Ue_),
    t_(ptf.t_),
    UMeanSpanTime_(ptf.UMeanSpanTime_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scaledMappedVelocityFixedValueFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    // Get the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        scaledMappedVelocityFixedValueFvPatchField::patch().patch()
    );
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());
    const word& fieldName = internalField().name();
    const volVectorField& UField = nbrMesh.lookupObject<volVectorField>(fieldName);
    const label nbrPatchID = nbrMesh.boundaryMesh().findPatchID(mpp.samplePatch());
    vectorField URecycled = UField.boundaryField()[nbrPatchID];

    // Read the predicted utau generated from wall model
    const volScalarField& uTauPredictedField = nbrMesh.lookupObject<volScalarField>("uTauPredicted");
    const label groundPatchID = nbrMesh.boundaryMesh().findPatchID("groundLeft"); 
    const scalarField uTauPredictedGround = uTauPredictedField.boundaryField()[groundPatchID];
    // Read ground patch face centres
    const vectorField& groundPatchFaceCentres = nbrMesh.Cf().boundaryField()[groundPatchID];
    //Info << "uTauPredictedGround: " << uTauPredictedGround << endl;
    //Info << "groundPatchFaceCentres: " << groundPatchFaceCentres << endl; 

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;
    //mpp.distribute(URecycled);

    // Get the number of cells in the wall-normal direction and the coordinate (random order)
    label Nnorm = 0;
    scalar newVal = 1.0;
    scalarField coord(patch().size(), 0);
    forAll(patch(), patchI)
    {
        newVal = 1.0;
        forAll(patch(), patchII)
        {
            if ( Nnorm == 0 ) { break; }
            if ( mag(patch().Cf()[patchI].component(1) - coord[patchII]) < 10e-9 )
            {
                newVal = 0;
                break;
            }
            if ( patchII == Nnorm - 1 ) { break; }
        }
        if ( newVal == 1.0 ) { coord[Nnorm] = patch().Cf()[patchI].component(1); Nnorm++; }
    }

    // Get the number of cells in the spanwise direction and the coordinate (random order)
    label Nspan = 0;
    scalarField coord2(patch().size(), 0);
    forAll(patch(), patchI)
    {
        newVal = 1.0;
        forAll(patch(), patchII)
        {
            if ( Nspan == 0 ) { break; }
            if ( mag(patch().Cf()[patchI].component(2) - coord2[patchII]) < 10e-9 )
            {
               newVal = 0;
               break;
            }
            if ( patchII == Nspan - 1 ) { break; }
        }
        if ( newVal == 1.0 ) { coord2[Nspan] = patch().Cf()[patchI].component(2); Nspan++; }
    }

    // Get the number of cells in the streamwise direction and the coordinate (random order)
    label Nstream = 0;
    scalarField coordX(groundPatchFaceCentres.size(), 0);
    forAll(groundPatchFaceCentres, patchI)
    {
        newVal = 1.0;
        forAll(groundPatchFaceCentres, patchII)
        {
            if ( Nstream == 0 ) { break; }
            if ( mag(groundPatchFaceCentres[patchI].component(0) - coordX[patchII]) < 10e-9 )
            {
               newVal = 0;
               break;
            }
            if ( patchII == Nstream - 1 ) { break; }
        }
        if ( newVal == 1.0 )
        {
            coordX[Nstream] = groundPatchFaceCentres[patchI].component(0); Nstream++;
        }
    }

    //Info << "Nstream: " << Nstream << endl;
    //Info << "coordX: " << coordX << endl;

    // Sort the normal coordinates from lowest to highest
    scalarField normCoord(Nnorm, 10e9);
    forAll(patch(), patchI)
    {
        forAll(patch(), patchII)
        {
            if ( coord[patchI] < normCoord[patchII] )    
            {
                forAll(patch(), patchIII)
                {
                    normCoord[Nnorm - 1 - patchIII] = normCoord[Nnorm - 2 - patchIII];
                    if (Nnorm - 1 - patchIII == patchII) { break; }
                }
                normCoord[patchII] = coord[patchI];
                break;
            }
        }
        if (patchI == Nnorm - 1) { break; }
    }

    // Sort the spanwise coordinates from lowest to highest
    scalarField spanCoord(Nspan, 10e9);
    forAll(patch(), patchI)
    {
        forAll(patch(), patchII)
        {
            if ( coord2[patchI] < spanCoord[patchII] )    
            {
                forAll(patch(), patchIII)
                {
                    spanCoord[Nspan - 1 - patchIII] = spanCoord[Nspan - 2 - patchIII];
                    if (Nspan - 1 - patchIII == patchII) { break; }
                }
                spanCoord[patchII] = coord2[patchI];
                break;
            }
        }
        if (patchI == Nspan - 1) { break; }
    }

    // Sort the streamwise coordinates from lowest to highest
    scalarField streamCoord(Nstream, 10e9);
    forAll(groundPatchFaceCentres, patchI)
    {
        forAll(groundPatchFaceCentres, patchII)
        {
            if ( coordX[patchI] < streamCoord[patchII] )    
            {
                forAll(patch(), patchIII)
                {
                    streamCoord[Nstream - 1 - patchIII] = streamCoord[Nstream - 2 - patchIII];
                    if (Nstream - 1 - patchIII == patchII) { break; }
                }
                streamCoord[patchII] = coordX[patchI];
                break;
            }
        }
        if (patchI == Nstream - 1) { break; }
    }

    // calculation of the friction velocity at the recycled plane from wall model
    scalar UTauRecycle = 0.2;

    forAll(groundPatchFaceCentres, patchI)
    {
       if (mag(groundPatchFaceCentres[patchI].component(0) - streamCoord[Nstream-1]) < 10e-9)
       {
           UTauRecycle += uTauPredictedGround[patchI];
       }
    }

    UTauRecycle /= Nspan;

    scalar UTauInlet = 0.2;

    forAll(groundPatchFaceCentres, patchI)
    {
       if (mag(groundPatchFaceCentres[patchI].component(0) - streamCoord[0]) < 10e-9)
       {
           UTauInlet += uTauPredictedGround[patchI];
       }
    }

    UTauInlet /= Nspan;


    //Info << "streamCoord: " << streamCoord << endl;
    //Info << "xMax: " << streamCoord[Nstream-1] << endl;
    //Info << "UTauRecycle: " << UTauRecycle << endl;

    // Get the normal label position
    labelField normLab(patch().size(), 0);
    forAll(patch(), patchI)
    {
        forAll(patch(), patchII)
        {
            if ( mag(patch().Cf()[patchI].component(1) - normCoord[patchII]) < 1e-9 )
            {
                normLab[patchI] = patchII;
                break;
            }
        }
    }

    // Get the spanwise label position
    labelField spanLab(patch().size(), 0);
    forAll(patch(), patchI)
    {
        forAll(patch(), patchII)
        {
            if ( mag(patch().Cf()[patchI].component(2) - spanCoord[patchII]) < 1e-9 )
            {
                spanLab[patchI] = patchII;
                break;
            }
        }
    }

    // average the instantaneous velocity in the recycle plane in the spanwise direction
    vectorField UMeanSpan(Nnorm, vector::zero);
    forAll(patch(), patchI)
    {
        UMeanSpan[normLab[patchI]] += URecycled[patchI];
    }
    UMeanSpan /= Nspan;

    // Read the time step and averaging time in controlDict
    const dictionary& controlDict = this->db().time().controlDict();
    scalar deltaT_(readScalar(controlDict.lookup("deltaT")));
    scalar T_(readScalar(controlDict.lookup("averagingTime")));

    // average the spanwise averaged velocity field in time
    // NB: many solvers call this function several times per time step. Hence, t_ is different to the actual simulation time!
    /*
    forAll(patch(), patchI)
    {
        if (t_/deltaT_<=1.0) {UMeanSpanTime_[patchI] = UMeanSpan[normLab[patchI]];} 
        else
        {
            UMeanSpanTime_[patchI] = (1 - (deltaT_/t_))*UMeanSpanTime_[patchI]  + (deltaT_/t_)*UMeanSpan[normLab[patchI]];
        }
    }
    t_ = t_ + deltaT_; // The averaging time increases progressively to remove the initial transient solution
    Info<< "Averaging time: " << t_ << endl;
    */

    // average the spanwise averaged velocity field in time
    t_ = this->db().time().value();
    if (t_ > T_) { t_ = T_;}
    forAll(patch(), patchI)
    {
        if (t_/deltaT_<=1.0) {UMeanSpanTime_[patchI] = UMeanSpan[normLab[patchI]];} 
        else
        {
            UMeanSpanTime_[patchI] = (1 - (deltaT_/t_))*UMeanSpanTime_[patchI]  + (deltaT_/t_)*UMeanSpan[normLab[patchI]];
        }
    }

    // average the spanwise averaged velocity field in time
    // Reynolds decomposition of the velocity field
    vectorField UPrime(patch().size(), vector::zero);
    forAll(patch(), patchI)
    {
        UPrime[patchI] = URecycled[patchI] - UMeanSpanTime_[patchI];
    }

    // Decompose the velocity vector (x_ is the flow direction in this case)
    vector x_(1,0,0);
    scalarField UxMeanSpanTime = (UMeanSpanTime_ & x_);
    vectorField UyzMeanSpanTime = UMeanSpanTime_ - (UxMeanSpanTime*x_);

    // Use UxMeanSpanTime with size Nnorm instead of patch().size(). It will be used to compute the boundary layer thickness, etc.
    scalarField UxMeanSpanTime_(Nnorm, 0);
    forAll(patch(), patchI)
    {
        UxMeanSpanTime_[normLab[patchI]] = UxMeanSpanTime[patchI];
    }
    //Info << UxMeanSpanTime_ << endl;

    // calculation of the friction velocity at the recycled plane
    // scalar UTauRecycle = 0;
    // UTauRecycle = sqrt(nu_*UxMeanSpanTime_[0]/normCoord[0]);

    // get the boundary layer thickness at the recycled plane
    scalar deltaRecycle = deltaInlet_;
    forAll(UxMeanSpanTime_, patchI)
    {
        if (UxMeanSpanTime_[patchI] > 0.99*Ue_)
        {
            deltaRecycle = (normCoord[patchI] - normCoord[patchI-1])/(UxMeanSpanTime_[patchI]-UxMeanSpanTime_[patchI-1])
                *(0.99*Ue_ - UxMeanSpanTime_[patchI-1]) + normCoord[patchI-1];
            break;
        }
    }
    //if (deltaRecycle<deltaInlet_) {deltaRecycle=1.1*deltaInlet_;}
 
    // trapezoidal integration for momentum thickness
    scalar thetaRecycle = 0;
    forAll(patch(), patchI)
    {
        if (patchI == 0)
        {
            thetaRecycle = 0.5*(UxMeanSpanTime_[patchI]/Ue_)*(1.0 - 0.5*(UxMeanSpanTime_[patchI]/Ue_))*(normCoord[patchI]);
        }
        else if (patchI < Nnorm)
        {
            thetaRecycle += 0.5*((UxMeanSpanTime_[patchI] + UxMeanSpanTime_[patchI-1])/Ue_)
                *(1.0 - 0.5*((UxMeanSpanTime_[patchI] + UxMeanSpanTime_[patchI-1])/Ue_))*(normCoord[patchI] - normCoord[patchI-1]);
        }
        else { break; }
    }

    // Info << "UxMeanSpanTime_ " << UxMeanSpanTime_ << endl;
    // power law rule relation for inlet friction velocity
    // scalar UTauInlet_ = UTauRecycle*pow((thetaRecycle/thetaInlet_), 0.125); // inlet friction velocity
    scalar UTauInlet_ = UTauInlet;  // from wall model
    Info << "deltaRecycle = " << deltaRecycle << "\t deltaInlet = " << deltaInlet_  << endl;
    Info << "thetaRecycle = " << thetaRecycle << "\t thetaInlet = " << thetaInlet_  << endl;
    Info << "UTauRecycle = " << UTauRecycle << "\t UTauInlet = " << UTauInlet_ << endl;

    // Find y target in the recycle plane (corresponding to etaInlet and yPlusInlet)
    scalarField yTargetInner(Nnorm, 0);
    scalarField yTargetOuter(Nnorm, 0);
    forAll(patch(), patchI)
    {
        if ( patchI < Nnorm)
        {
            yTargetInner[patchI] = (UTauInlet_/UTauRecycle)*normCoord[patchI];
            if ( yTargetInner[patchI] > normCoord[Nnorm-1] ) { yTargetInner[patchI] = normCoord[Nnorm-1]; }

            yTargetOuter[patchI] = (deltaRecycle/deltaInlet_)*normCoord[patchI];
            if ( yTargetOuter[patchI] > normCoord[Nnorm-1] ) { yTargetOuter[patchI] = normCoord[Nnorm-1]; }
        }
        else { break; }
    }

    // Find the lower and higher bounds of the y target
    labelField yLowInner(Nnorm, 0);
    labelField yLowOuter(Nnorm, 0);
    labelField yHighInner(Nnorm, 1);
    labelField yHighOuter(Nnorm, 1);
    scalarField fractionInner(Nnorm, 0);
    scalarField fractionOuter(Nnorm, 0);
    forAll(patch(), patchI)
    {
        if (patchI == Nnorm) { break;}

        forAll(patch(), patchII)
        {
            if (normCoord[patchII] <= yTargetInner[patchI])
            {
                yLowInner[patchI] = patchII;
                if (patchII == Nnorm - 1){ yHighInner[patchI] = patchII; break; }
                else { yHighInner[patchI] = patchII + 1; }
            }
            else { break; }
        }

        forAll(patch(), patchII)
        {
            if (normCoord[patchII] <= yTargetOuter[patchI])
            {
                yLowOuter[patchI] = patchII;
                if (patchII == Nnorm - 1){ yHighOuter[patchI] = patchII; break; }
                else { yHighOuter[patchI] = patchII + 1; }
            }
	    else { break; }
        }

        // Also compute the fraction within the interval
        if( yLowInner[patchI] == Nnorm - 1 ) { fractionInner[patchI] = 0.0; }
        else
        {
            fractionInner[patchI] = (yTargetInner[patchI] - normCoord[yLowInner[patchI]])
                                /(normCoord[yHighInner[patchI]] - normCoord[yLowInner[patchI]]);
        }

        if( yLowOuter[patchI] == Nnorm - 1 ) { fractionOuter[patchI] = 0.0; }
        else
        {
            fractionOuter[patchI] = (yTargetOuter[patchI] - normCoord[yLowOuter[patchI]])
                                /(normCoord[yHighOuter[patchI]] - normCoord[yLowOuter[patchI]]);
        }
    }

    // Find neighbors above and below yTargetInnet (in the normal direction)
    labelField lowInnerLab(patch().size(), 0);
    labelField highInnerLab(patch().size(), 0);
    labelField lowOuterLab(patch().size(), 0);
    labelField highOuterLab(patch().size(), 0);
    forAll(patch(), patchI)
    {
        forAll(patch(), patchII)
        {
            if ( spanLab[patchI] == spanLab[patchII] && normLab[patchII] == yLowInner[normLab[patchI]] )
            {
                lowInnerLab[patchI] = patchII;
            }
            if ( spanLab[patchI] == spanLab[patchII] && normLab[patchII] == yHighInner[normLab[patchI]] )
            {
                highInnerLab[patchI] = patchII;
            }
            if ( spanLab[patchI] == spanLab[patchII] && normLab[patchII] == yLowOuter[normLab[patchI]] )
            {
                lowOuterLab[patchI] = patchII;
            }
            if ( spanLab[patchI] == spanLab[patchII] && normLab[patchII] == yHighOuter[normLab[patchI]] )
            {
                highOuterLab[patchI] = patchII;
            }
        }
    }

    // Interpolate velocity components
    vectorField UPrimeInner = UPrime;
    vectorField UPrimeOuter = UPrime;
    scalarField UxMeanInner = UxMeanSpanTime;
    scalarField UxMeanOuter = UxMeanSpanTime;
    vectorField UyzMeanInner = UyzMeanSpanTime;
    vectorField UyzMeanOuter = UyzMeanSpanTime;
    forAll(patch(), patchI)
    {
        UPrimeInner[patchI] = UPrime[lowInnerLab[patchI]] + (UPrime[highInnerLab[patchI]] - UPrime[lowInnerLab[patchI]])*fractionInner[normLab[patchI]];
        UPrimeOuter[patchI] = UPrime[lowOuterLab[patchI]] + (UPrime[highOuterLab[patchI]] - UPrime[lowOuterLab[patchI]])*fractionOuter[normLab[patchI]];
        UxMeanInner[patchI] = UxMeanSpanTime[lowInnerLab[patchI]] + (UxMeanSpanTime[highInnerLab[patchI]] - UxMeanSpanTime[lowInnerLab[patchI]])*fractionInner[normLab[patchI]];
        UxMeanOuter[patchI] = UxMeanSpanTime[lowOuterLab[patchI]] + (UxMeanSpanTime[highOuterLab[patchI]] - UxMeanSpanTime[lowOuterLab[patchI]])*fractionOuter[normLab[patchI]];
        UyzMeanInner[patchI] = UyzMeanSpanTime[lowInnerLab[patchI]] + (UyzMeanSpanTime[highInnerLab[patchI]] - UyzMeanSpanTime[lowInnerLab[patchI]])*fractionInner[normLab[patchI]];
        UyzMeanOuter[patchI] = UyzMeanSpanTime[lowOuterLab[patchI]] + (UyzMeanSpanTime[highOuterLab[patchI]] - UyzMeanSpanTime[lowOuterLab[patchI]])*fractionOuter[normLab[patchI]];
    }

    // Scale velocity components
    forAll(patch(), patchI)
    {
        UPrimeInner[patchI] = UPrimeInner[patchI]*(UTauInlet_/UTauRecycle);
        UPrimeOuter[patchI] = UPrimeOuter[patchI]*(UTauInlet_/UTauRecycle);
        UxMeanInner[patchI] = UxMeanInner[patchI]*(UTauInlet_/UTauRecycle);
        UxMeanOuter[patchI] = (UTauInlet_/UTauRecycle)*UxMeanOuter[patchI] + (1.0 - (UTauInlet_/UTauRecycle))*Ue_;
    }

    // Define the Inner-Outer blending function W
    scalarField etaInlet(patch().size(), 0);
    forAll(patch(), patchI)
    {
        etaInlet[patchI] = patch().Cf()[patchI].component(1)/deltaInlet_;
    }
    scalarField W = 0.5*(1.0+tanh(4.0*(etaInlet-0.2)/((1.0-2*0.2)*etaInlet+0.2))/tanh(4.0)); // Blending function

    forAll(patch(), patchI)
    {
        if( W[patchI] > 1.0 ){ W[patchI] = 1.0; }
    }

    // Blend the inner and outer velocity profiles
    scalarField UxMean(patch().size(), 0);
    vectorField UyzMean(patch().size(), vector::zero);
    forAll(patch(), patchI)
    {
        UPrime[patchI] = (1.0-W[patchI])*UPrimeInner[patchI] + W[patchI]*UPrimeOuter[patchI];
        UxMean[patchI] = (1.0-W[patchI])*UxMeanInner[patchI] + W[patchI]*UxMeanOuter[patchI];
        UyzMean[patchI] = (1.0-W[patchI])*UyzMeanInner[patchI] + W[patchI]*UyzMeanOuter[patchI];

        if ( UxMean[patchI] > Ue_ ) { UxMean[patchI] = Ue_; }
    }

    // reconstruct the scaled velocity field
    URecycled = UyzMean + (UxMean*x_) + UPrime;
    forAll(patch(), patchI)
    {
        if (normCoord[normLab[patchI]] > 2*deltaInlet_) // Damp perturbations in the freestream
        {
             URecycled[patchI] += (Ue_ - (URecycled[patchI] & x_))*x_;
        }
    }

    // construct UxMean of size Nnorm instead of patch().size() and sort it.
    scalarField UxMean_(Nnorm, 0);
    forAll(patch(), patchI)
    {
        UxMean_[normLab[patchI]] = UxMean[patchI];
    }

    // trapezoidal integration for inlet momentum thickness
    forAll(patch(), patchI)
    {
        if (patchI == 0)
        {
            thetaInlet_ = 0.5*(UxMean_[patchI]/Ue_)*(1.0 - 0.5*(UxMean_[patchI]/Ue_))*(normCoord[patchI]);
        }
        else if (patchI < Nnorm)
        {
            thetaInlet_ += 0.5*((UxMean_[patchI] + UxMean_[patchI-1])/Ue_)
                *(1.0 - 0.5*((UxMean_[patchI] + UxMean_[patchI-1])/Ue_))*(normCoord[patchI] - normCoord[patchI-1]);
        }
        else { break; }
    }


    operator==(URecycled);

    fixedValueFvPatchVectorField::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}


void Foam::scaledMappedVelocityFixedValueFvPatchField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "value", *this);
    os.writeKeyword("deltaInlet") << deltaInlet_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaInlet") << thetaInlet_ << token::END_STATEMENT << nl;
    os.writeKeyword("nu") << nu_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ue") << Ue_ << token::END_STATEMENT << nl;
    os.writeKeyword("t") << t_ << token::END_STATEMENT << nl;
    writeEntry(os, "UMeanSpanTime", UMeanSpanTime_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        scaledMappedVelocityFixedValueFvPatchField
    );
}


// ************************************************************************* //
