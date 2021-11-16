wordList psiBCTypes
(
    alpha1.boundaryField().size(),
    zeroGradientFvPatchScalarField::typeName
);

forAll(alpha1.boundaryField(), patchi)
{
    if (alpha1.boundaryField()[patchi].fixesValue())
    {
        psiBCTypes[patchi] = fixedValueFvPatchScalarField::typeName;
    }
}

volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero),
    psiBCTypes
);
    
interfacePropertiesLS interfaceLS(psi, alpha1, U, mixture);



