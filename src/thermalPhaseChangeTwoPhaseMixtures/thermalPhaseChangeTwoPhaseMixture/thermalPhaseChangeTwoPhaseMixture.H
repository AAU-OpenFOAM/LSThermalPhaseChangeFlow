/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::thermalPhaseChangeTwoPhaseMixture

Description

SourceFiles
    thermalPhaseChangeTwoPhaseMixture.C
    newPhaseChangeModel.C
    
    HArdt and Lee models here and in this solver is based on the work the solver phaseChangeHeatFoam by N. Samkhaniani
    
    [N. Samkhaniani, M.R. Ansari, "The evaluation of the diffuse interface method for phase change simulations using OpenFOAM", Heat Transfer–Asian Research. 2017;00:1–31. doi:10.1002/htj.21268](http://onlinelibrary.wiley.com/doi/10.1002/htj.21268/full)

\*---------------------------------------------------------------------------*/

#ifndef thermalPhaseChangeTwoPhaseMixture_H
#define thermalPhaseChangeTwoPhaseMixture_H

#include "thermalIncompressibleTwoPhaseMixture.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class thermalPhaseChangeTwoPhaseMixture Declaration
\*---------------------------------------------------------------------------*/

class thermalPhaseChangeTwoPhaseMixture
:
    public thermalIncompressibleTwoPhaseMixture
{

protected:

    // Protected data

        dictionary thermalPhaseChangeTwoPhaseMixtureCoeffs_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;
        //- Saturation vapour Temperature
        dimensionedScalar TSat_; //add
        //- true if Saturation vapour Temperature based on pressure
        bool TSatLocal_; //add
        
        //- gas specific constant
        dimensionedScalar R_;//add
        
        dimensionedScalar Hfg_;//add



    // Private Member Functions

        //- No copy construct
        thermalPhaseChangeTwoPhaseMixture(const thermalPhaseChangeTwoPhaseMixture&) = delete;

        //- No copy assignment
        void operator=(const thermalPhaseChangeTwoPhaseMixture&) = delete;


public:

    //- Runtime type information
    TypeName("thermalPhaseChangeTwoPhaseMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            thermalPhaseChangeTwoPhaseMixture,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<thermalPhaseChangeTwoPhaseMixture> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        thermalPhaseChangeTwoPhaseMixture
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~thermalPhaseChangeTwoPhaseMixture() = default;


    // Member Functions

        //- Return const-access to the saturation vapour Temperature
        const dimensionedScalar& TSat() const //add
        {
        	return	TSat_;
        }
        
        
        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }
        
        //- Return const-access to the saturation vapour Temperature
        volScalarField TSatLocal() const;
        
        //- Return const-access to the latent heat
        const dimensionedScalar& Hfg() const //add
        {
        	return Hfg_;
        }

        

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const = 0;
        
        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (T - TSat)        
        virtual Pair<tmp<volScalarField>> mDotT() const = 0;

        //- Return the volumetric condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        Pair<tmp<volScalarField>> vDotAlphal() const;

        //- Return the volumetric condensation and vaporisation rates as
        //  coefficients to multiply (p - pSat)
        Pair<tmp<volScalarField>> vDotP() const;
        
        Pair<tmp<volScalarField> > vDotT() const;
        

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
