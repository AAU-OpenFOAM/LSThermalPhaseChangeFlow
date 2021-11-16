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

#include "interfacePropertiesLSH.H"
#include "alphaContactAngleTwoPhaseFvPatchScalarField.H"
#include "psiContactAngleTwoPhaseFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nVec on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfacePropertiesLS::correctContactAngle
(
    surfaceVectorField::Boundary& nVecb,
    const surfaceVectorField::Boundary& gradPsif
) const
{
    const fvMesh& mesh = psi_.mesh();
    const volScalarField::Boundary& abf = psi_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<psiContactAngleTwoPhaseFvPatchScalarField>(abf[patchi]))
        {
            psiContactAngleTwoPhaseFvPatchScalarField& acap =
                const_cast<psiContactAngleTwoPhaseFvPatchScalarField&>
                (
                    refCast<const psiContactAngleTwoPhaseFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nVecp = nVecb[patchi];
            const scalarField theta
            (
                degToRad() * acap.theta(U_.boundaryField()[patchi], nVecp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nVecp to correspond to the contact angle

            const scalarField a12(nVecp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nVecp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nVecp = a*nf + b*nVecp;
            nVecp /= (mag(nVecp) + deltaN_.value());

            acap.gradient() = (nf & nVecp)*mag(gradPsif[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfacePropertiesLS::calculateC()
{
    const fvMesh& mesh = psi_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of psi
    const volVectorField gradPsi(fvc::grad(psi_, "nVec"));

    // Interpolated face-gradient of psi
    surfaceVectorField gradPsif(fvc::interpolate(gradPsi));

    // Face unit interface normal
    surfaceVectorField nVecfv(gradPsif/(mag(gradPsif) + /*deltaN_+*/scalar(1.0e-6)/dimChange_));
    
    correctContactAngle(nVecfv.boundaryFieldRef(), gradPsif.boundaryField());  

    // Face unit interface normal flux
    nVecf_ = nVecfv & Sf;

    // Simple expression for curvature
    C_ = -fvc::div(nVecf_);

}


void Foam::interfacePropertiesLS::calculatePsi0()   // construct of psi0 from alpha
{
    psi0_ == (double(2.0)*alpha1_ - double(1.0))*gamma_;
}

void Foam::interfacePropertiesLS::calculateDelta()   // calculate Dirac function
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(mag(psi_[celli]) > epsilon_.value())
          delta_[celli] = double(0.0);
       else
          delta_[celli] = double(1.0)/(double(2.0)*epsilon_.value())*(double(1.0)+cos(M_PI*psi_[celli]/epsilon_.value()));
    }
}

void Foam::interfacePropertiesLS::calculateH()  // calculate the heaviside function 
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(psi_[celli] < -epsilon_.value())
          H_[celli] = double(0.0);
       else if(epsilon_.value() < psi_[celli])
          H_[celli] = double(1.0);
       else
          H_[celli] = double(1.0)/double(2.0)*(double(1.0)+psi_[celli]/epsilon_.value()+sin(M_PI*psi_[celli]/epsilon_.value())/M_PI);
    }
}

void Foam::interfacePropertiesLS::calculateHscale()
{
    forAll(psi_.mesh().cells(),celli)
    {
       if(psi_[celli] < -epsilon_.value())
          Hscale_[celli] = double(0.0);
       else if(epsilon_.value() < psi_[celli])
          Hscale_[celli] = double(1.0);
       else
          Hscale_[celli] = double(1.0)/double(2.0)*(double(1.0)/double(2.0)+psi_[celli]/epsilon_.value()+psi_[celli]*psi_[celli]/(double(2.0)*epsilon_.value()*epsilon_.value())-(cos(double(2.0)*M_PI*psi_[celli]/epsilon_.value())-double(1.0))/(double(4.0)*M_PI*M_PI)+sin(M_PI*psi_[celli]/epsilon_.value())*(epsilon_.value()+psi_[celli])/(M_PI*epsilon_.value()));
    }
}

void Foam::interfacePropertiesLS::calculateDeltaScale()
{
    deltaScale_ == double(2.0)*H_*delta_;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacePropertiesLS::interfacePropertiesLS
(
    const volScalarField& psi,
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        psi.mesh().solverDict(psi.name()).get<scalar>("cAlpha")
    ),

    sigmaPtr_(surfaceTensionModel::New(dict, alpha1.mesh())),
    
    deltaX_("deltaX", dimensionSet(0, 0, 0, 0, 0), dict),  // deltaX which is used in construct of psi0 from alpha
    
    gamma_
    (
        "gamma",    // used in construct of psi0 from alpha
        deltaX_*double(0.75)
    ),

    epsilon_
    (
        "epsilon",  // used in both dirac function and heaviside function
        deltaX_*double(1.5)
    ),

    deltaTau_
    (
        "deltaTau",  // this is used in re-initialization of psi
        deltaX_*double(0.1)
    ),
    
    dimChange_
    (
        dimensionedScalar("dimChange",dimLength, 1.0)
    ),


    deltaN_
    (
        "deltaN",
        1e-8/cbrt(average(psi.mesh().V()))
    ),
    
    psi_(psi),

    alpha1_(alpha1),
    
    U_(U),

    nVecf_
    (
        IOobject
        (
            "nVecf",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar(dimArea, Zero)
    ),
    
     nVecfv_
    (
        IOobject
        ( 
            "nVecfv", // normal vector of interface
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedVector(dimless, vector::zero)
    ),

    C_
    (
        IOobject
        (
            "interfacePropertiesLS:C",
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar(dimless/dimLength, Zero)
    ),
    
     psi0_
    (
        IOobject
        (
            "psi0",  // Initial level set
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("psi0", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    delta_
    (
        IOobject
        (
            "delta",  // Dirac function
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("delta", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    ),

    H_
    (
        IOobject
        (
            "H",   //Heaviside function
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("H", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    Hscale_
    (
        IOobject
        (
            "Hscale", // magnitude of the heaviside function
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("Hscale", dimless, 0.0),
        psi_.boundaryField().types()
    ),

    deltaScale_
    (
        IOobject
        (
            "deltaScale", // magnitude of dirac function
            psi_.time().timeName(),
            psi_.mesh()
        ),
        psi_.mesh(),
        dimensionedScalar("deltaScale", dimless, 0.0),
        calculatedFvPatchScalarField::typeName
    )
{
    calculateC();
    calculatePsi0();
    calculateDelta();
    calculateH();
    calculateHscale();
    calculateDeltaScale();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::interfacePropertiesLS::sigmaC() const
{
    return sigmaPtr_->sigma()*C_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfacePropertiesLS::surfaceTensionForce() const
{
    return fvc::snGrad(psi_)*fvc::interpolate(sigmaC())*fvc::interpolate(delta());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfacePropertiesLS::surfaceTensionForceScale() const
{
    return fvc::snGrad(psi_)*fvc::interpolate(sigmaC())*fvc::interpolate(deltaScale());
}

void Foam::interfacePropertiesLS::correct()
{
    calculateC();
}


bool Foam::interfacePropertiesLS::read()
{
    psi_.mesh().solverDict(psi_.name()).readEntry("cAlpha", cAlpha_);
    sigmaPtr_->readDict(transportPropertiesDict_);

    return true;
}


// ************************************************************************* //
