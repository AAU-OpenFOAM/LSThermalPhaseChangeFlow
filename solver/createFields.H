Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));

volScalarField mDotT
(
    IOobject
    (
        "mDotT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("mDotT",dimDensity/dimTime,0)
);

volScalarField mDotP
(
    IOobject
    (
        "mDotP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("mDotP",dimDensity/dimTime,0)
);

volScalarField mDotAlpha
(
    IOobject
    (
        "mDotAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("mDotAlpha",dimDensity/dimTime,0)
);

#include "createPhi.H"


Info<< "Creating thermalPhaseChangeTwoPhaseMixture\n" << endl;
autoPtr<thermalPhaseChangeTwoPhaseMixture> mixture =
    thermalPhaseChangeTwoPhaseMixture::New(U, phi);

volScalarField& alpha1(mixture->alpha1());
volScalarField& alpha2(mixture->alpha2());

const dimensionedScalar& rho1 = mixture->rho1();
const dimensionedScalar& rho2 = mixture->rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Construct interface from alpha1 distribution
interfaceProperties interface(alpha1, U, mixture());

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, mixture())
);


IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        U.time().constant(),
        U.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);



surfaceForces surfForces(alpha1,phi,U,transportProperties);

const volScalarField& gh = surfForces.acc();


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createFvOptions.H"

IOobject alphaPhi10Header
(
    IOobject::groupName("alphaPhi0", alpha1.group()),
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

// MULES flux from previous time-step
surfaceScalarField alphaPhi10
(
    alphaPhi10Header,
    phi*fvc::interpolate(alpha1)
);

// Looking up which alpha1 solver to use
const dictionary& alpha1Dict(mesh.solverDict(alpha1.name()));

word interfaceMethod
(
    alpha1Dict.lookupOrDefault<word>("interfaceMethod", "MULES")
);

#include "createIsoAdvection.H"

#include "createLSFields.H"
