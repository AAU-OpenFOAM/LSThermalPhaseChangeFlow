/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::interfacePropertiesLS

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfacePropertiesLS.C

\*---------------------------------------------------------------------------*/

#ifndef interfacePropertiesLS_H
#define interfacePropertiesLS_H

#include "IOdictionary.H"
#include "surfaceTensionModel.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class interfacePropertiesLS Declaration
\*---------------------------------------------------------------------------*/

class interfacePropertiesLS
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;

        //- Compression coefficient
        scalar cAlpha_;

        //- Surface tension
        autoPtr<surfaceTensionModel> sigmaPtr_;

        //- Parameter for S-CLSVOF method (deltaX,gamma,epsilon,deltaTau)
        dimensionedScalar deltaX_;
        dimensionedScalar gamma_;
        dimensionedScalar epsilon_;
        dimensionedScalar deltaTau_;
        dimensionedScalar dimChange_;


        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& psi_;
        const volScalarField& alpha1_;
        const volVectorField& U_;
        surfaceScalarField nVecf_;
        surfaceVectorField nVecfv_;
        volScalarField C_;
        volScalarField psi0_;
        volScalarField delta_;
        volScalarField H_;
        volScalarField Hscale_;
        volScalarField deltaScale_;


    // Private Member Functions

        //- No copy construct
        interfacePropertiesLS(const interfacePropertiesLS&) = delete;

        //- No copy assignment
        void operator=(const interfacePropertiesLS&) = delete;

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::Boundary& nVec,
            const surfaceVectorField::Boundary& gradPsif
        ) const;

        //- Re-calculate the interface curvature
        void calculateC();

        //- Re-calculate the delta function and Heaviside function
        void calculateDelta();
        void calculateH();
        void calculateHscale();
        void calculateDeltaScale();

        
        //- Re-calculate the interface curvature
        void calculatePsi0();


public:

    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfacePropertiesLS
        (
            const volScalarField& psi,
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions
        
        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const surfaceScalarField& nVecf() const
        {
            return nVecf_;
        }
        
       const surfaceVectorField& nVecfv() const
        {
            return nVecfv_;
        }

        const volScalarField& C() const
        {
            return C_;
        }

          const dimensionedScalar& epsilon() const
        {
            return epsilon_;
        }

        const dimensionedScalar& deltaTau() const
        {
            return deltaTau_;
        }
        
        const dimensionedScalar& dimChange() const
        {
            return dimChange_;
        }

        const volScalarField& psi0() const
        {
            return psi0_;
        }

        const volScalarField& delta() const
        {
            return delta_;
        }

        const volScalarField& H() const
        {
            return H_;
        }

        const volScalarField& Hscale() const
        {
            return Hscale_;
        }

        const volScalarField& deltaScale() const
        {
            return deltaScale_;
        }

        tmp<volScalarField> sigmaC() const;

        tmp<surfaceScalarField> surfaceTensionForce() const;
        
        tmp<surfaceScalarField> surfaceTensionForceScale() const;


        void correct();
        
        void correctPsi0()
        {
            calculatePsi0();
        }

        void correctDeltaH()
        {
            calculateDelta();
            calculateH();
            delta_.correctBoundaryConditions();
            H_.correctBoundaryConditions();
        }

        void correctDeltaHscale()
        {
            calculateDelta();
            calculateH();
            calculateHscale();
            calculateDeltaScale();
            delta_.correctBoundaryConditions();
            H_.correctBoundaryConditions();
            deltaScale_.correctBoundaryConditions();
            Hscale_.correctBoundaryConditions();
        }

        //- Read transportProperties dictionary
        bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
