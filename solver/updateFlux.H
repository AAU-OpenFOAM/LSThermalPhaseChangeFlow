{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // Standard face-flux compression coefficient
    surfaceScalarField phic(interface.cAlpha()*mag(phi/mesh.magSf()));
    volScalarField divU(fvc::div(phi));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (interface.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    surfaceScalarField phir("phir", phic*interfaceLS.nVecf());

    Pair<tmp<volScalarField>> vDotAlphal = mixture->vDotAlphal();
        
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);

    surfaceScalarField phiAlpha
    (
        fvc::flux
        (
            phi,
            alpha1,
            alphaScheme
        )
        + fvc::flux
        (
            -fvc::flux(-phir, alpha2, alpharScheme),
            alpha1,
            alpharScheme
        )
    );

    MULES::explicitSolve
    (
        geometricOneField(),
        alpha1,
        phi,
        phiAlpha,
        vDotvmcAlphal,
        (divU*alpha1 + vDotcAlphal)(),
        oneField(),
        zeroField()
    );
    rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
}
