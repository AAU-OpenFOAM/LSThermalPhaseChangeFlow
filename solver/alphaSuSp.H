volScalarField divU(fvc::div(phi));
volScalarField Su(divU*alpha1 + vDotcAlphal);
volScalarField Sp(vDotvmcAlphal);
