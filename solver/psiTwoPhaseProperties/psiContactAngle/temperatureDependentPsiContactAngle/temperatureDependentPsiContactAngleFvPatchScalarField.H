/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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
    Foam::temperatureDependentPsiContactAngleFvPatchScalarField

Description
    Temperature-dependent psiContactAngle scalar boundary condition.

    The contact angle is evaluated from the specified Foam::Function1 for the
    temperature field looked-up from the mesh database the name of which
    may optionally be provided.

Usage
    \table
        Property     | Description             | Required    | Default value
        T            | Temperature field name  | no          | T
        theta0       | Contact angle function  | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            temperatureDependentPsiContactAngle;
        theta0          constant 60;
    }
    \endverbatim

See also
    Foam::psiContactAngleTwoPhaseFvPatchScalarField
    Foam::constantPsiContactAngleFvPatchScalarField
    Foam::Function1

SourceFiles
    temperatureDependentPsiContactAngleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureDependentPsiContactAngleFvPatchScalarField_H
#define temperatureDependentPsiContactAngleFvPatchScalarField_H

#include "psiContactAngleTwoPhaseFvPatchScalarField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
       Class temperatureDependentPsiContactAngleFvPatch Declaration
\*---------------------------------------------------------------------------*/

class temperatureDependentPsiContactAngleFvPatchScalarField
:
    public psiContactAngleTwoPhaseFvPatchScalarField
{
    // Private data

        //- Name of temperature field, default = "T"
        word TName_;

        //- Equilibrium contact angle function
        autoPtr<Function1<scalar>> theta0_;


public:

    //- Runtime type information
    TypeName("temperatureDependentPsiContactAngle");


    // Constructors

        //- Construct from patch and internal field
        temperatureDependentPsiContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperatureDependentPsiContactAngleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  temperatureDependentPsiContactAngleFvPatchScalarField
        //  onto a new patch
        temperatureDependentPsiContactAngleFvPatchScalarField
        (
            const temperatureDependentPsiContactAngleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        temperatureDependentPsiContactAngleFvPatchScalarField
        (
            const temperatureDependentPsiContactAngleFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new temperatureDependentPsiContactAngleFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        temperatureDependentPsiContactAngleFvPatchScalarField
        (
            const temperatureDependentPsiContactAngleFvPatchScalarField&,
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
                new temperatureDependentPsiContactAngleFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Return the equilibrium contact-angle
        virtual tmp<scalarField> theta
        (
            const fvPatchVectorField& Up,
            const fvsPatchVectorField& nVec
        ) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
