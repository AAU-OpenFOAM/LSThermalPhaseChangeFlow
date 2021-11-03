Info << "solve TEqn" << endl;
{
    const dimensionedScalar& K1 = mixture->kappa1();
    const dimensionedScalar& K2 = mixture->kappa2();

    const dimensionedScalar& C1 = mixture->Cp1();
    const dimensionedScalar& C2 = mixture->Cp2();    
    const volScalarField& TSatLocal = mixture->TSatLocal();

    Pair<tmp<volScalarField> > vDotT = mixture->vDotT();
        const volScalarField& vDotcT = vDotT[0]();
        const volScalarField& vDotvT = vDotT[1]();

    volScalarField k= K1*alpha1 + (1.0 - alpha1)*K2;
    volScalarField rhoC = rho1*C1*alpha1 + (1.0 - alpha1)*rho2*C2;

    fvScalarMatrix TEqn
    (
          fvm::ddt(T)
        + fvm::div(phi, T)
        - fvm::Sp(fvc::div(phi),T)
        - fvm::laplacian( k/rhoC , T,"laplacian(alphaEff,T)")
    );

    TEqn.relax();
    fvOptions.constrain(TEqn);
        
        
    solve
    (
        TEqn
        ==
        fvm::Sp(vDotcT-vDotvT,T) - (vDotcT-vDotvT)*TSatLocal
    );

    fvOptions.correct(T);
    
    Info<< "TAve = "
        << T.weightedAverage(mesh.V()).value()
        << "  Min(T) = " << min(T).value()
        << "  Max(T) = " << max(T).value()
        << endl;

    Pair<tmp<volScalarField> > mDot = mixture->mDotT();
    const volScalarField& mDotcT = mDot[0]();
    const volScalarField& mDotvT = mDot[1]();
    mDotT = (mDotcT-mDotvT)*(T-TSatLocal);
}

