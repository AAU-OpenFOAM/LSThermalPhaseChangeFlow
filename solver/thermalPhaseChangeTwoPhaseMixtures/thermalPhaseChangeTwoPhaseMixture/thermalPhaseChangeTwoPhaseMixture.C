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

\*---------------------------------------------------------------------------*/

#include "thermalPhaseChangeTwoPhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalPhaseChangeTwoPhaseMixture, 0);
    defineRunTimeSelectionTable(thermalPhaseChangeTwoPhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalPhaseChangeTwoPhaseMixture::thermalPhaseChangeTwoPhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    thermalIncompressibleTwoPhaseMixture(U, phi),
    thermalPhaseChangeTwoPhaseMixtureCoeffs_(optionalSubDict(type + "Coeffs")),
    pSat_("pSat", dimPressure, *this),
    TSat_("TSat", dimTemperature, *this),
    TSatLocal_(readBool(lookup("TSatLocal"))),
    R_("R", pow2(dimLength)/pow2(dimTime)/dimTemperature, *this),
    Hfg_("Hfg",Hf2() - Hf1())


{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::thermalPhaseChangeTwoPhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff(1.0/rho1() - alpha1_*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::thermalPhaseChangeTwoPhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho1() - 1.0/rho2());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}


Foam::Pair<Foam::tmp<Foam::volScalarField> >
Foam::thermalPhaseChangeTwoPhaseMixture::vDotT() const
{
//         dimensionedScalar Hfg = Hf2() - Hf1();
        volScalarField rhoCp =  rho1()*Cp1()*alpha1() + rho2()*Cp2()*(1.0-alpha1());
        //volScalarField TCoeff = (Hfg()+(Cp1()-Cp2())*TSat_)/rhoCp;
        volScalarField TCoeff = Hfg_/rhoCp;
        //volScalarField TCoeff = Hfg()*(alpha1_/Cp1() + (1.0 - alpha1_)/Cp2());
        Pair<tmp<volScalarField> > mDotT = this->mDotT();

            return Pair<tmp<volScalarField> >
            (
                TCoeff*mDotT[0],
                TCoeff*mDotT[1]
            );
}

Foam::volScalarField
Foam::thermalPhaseChangeTwoPhaseMixture::TSatLocal() const
{
    if (TSatLocal_)
    {
//         dimensionedScalar Hfg = Hf2() - Hf1();
        //Info <<"TSatlocal" << endl;
        const volScalarField& p = alpha1().db().lookupObject<volScalarField>("p");
        volScalarField oneByT = 1.0/TSat() - R_/Hfg_*log(max(p/pSat_,1E-08));

        return volScalarField
        (
            1.0/oneByT
        );
    }
    else
    {
        //Info <<"TSat" << endl;
        volScalarField one
        (
            IOobject
            (
                "one",
                alpha1().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            alpha1().mesh(),
            dimensionedScalar ("one",dimless, 1.0)
        );

        return volScalarField
        (
            one*TSat_
        );
    }
}

bool Foam::thermalPhaseChangeTwoPhaseMixture::read()
{
    if (thermalIncompressibleTwoPhaseMixture::read())
    {
        thermalPhaseChangeTwoPhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");
        readEntry("pSat", pSat_);
        readEntry("TSat", TSat_);
        readEntry("TSatLocal", TSatLocal_);
        readEntry("R", R_);


        return true;
    }

    return false;
}


// ************************************************************************* //
