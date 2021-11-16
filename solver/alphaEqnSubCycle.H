{
    if (interfaceMethod == "isoAdvector" && pimple.nCorrPIMPLE() > 1)
    {
        // If nOuterCorrectors > 1 then for all but the first loop the advection
        // of alpha is done using an average, 0.5*phi+0.5*phiNew where phi is
        // the flux at the beginning of the time step and phiNew is the flux
        // estimate at the end of the time step from the previous outer
        // iteration. Similarly we use 0.5*U + 0.5*UNew in later iterations.
        if (pimple.firstIter())
        {
            // To recalculate the alpha1 update in subsequent iterations, we
            // must store its current value before overwriting with the new
            // value
            alpha1.storePrevIter();
            // Storing initial phi and U for use in later outer iterations.
            phi.storePrevIter();
            U.storePrevIter();
        }
        else
        {
            // Resetting alpha1 to value before advection in first PIMPLE
            // iteration.
            alpha1 = alpha1.prevIter();

            // Temporarily setting U and phi with which to advect interface.
            U = 0.5*U.prevIter() + 0.5*U;
            phi = 0.5*phi.prevIter() + 0.5*phi;
        }
    }
    
    if (interfaceMethod == "isoAdvectorPLIC" && pimple.nCorrPIMPLE() > 1)
    {
        // If nOuterCorrectors > 1 then for all but the first loop the advection
        // of alpha is done using an average, 0.5*phi+0.5*phiNew where phi is
        // the flux at the beginning of the time step and phiNew is the flux
        // estimate at the end of the time step from the previous outer
        // iteration. Similarly we use 0.5*U + 0.5*UNew in later iterations.
        if (pimple.firstIter())
        {
            // To recalculate the alpha1 update in subsequent iterations, we
            // must store its current value before overwriting with the new
            // value
            alpha1.storePrevIter();
            // Storing initial phi and U for use in later outer iterations.
            phi.storePrevIter();
            U.storePrevIter();
        }
        else
        {
            // Resetting alpha1 to value before advection in first PIMPLE
            // iteration.
            alpha1 = alpha1.prevIter();

            // Temporarily setting U and phi with which to advect interface.
            U = 0.5*U.prevIter() + 0.5*U;
            phi = 0.5*phi.prevIter() + 0.5*phi;
        }
    }
    
    
    // Standard face-flux compression coefficient
    surfaceScalarField phic(interface.cAlpha()*mag(phi/mesh.magSf()));

    // Add the optional isotropic compression contribution
    if (icAlpha > 0)
    {
        phic *= (1.0 - icAlpha);
        phic += (interface.cAlpha()*icAlpha)*fvc::interpolate(mag(U));
    }

    volScalarField divU(fvc::div(phi));
    
    Pair<tmp<volScalarField> > mDot = mixture->mDotAlphal();
    const volScalarField& mDotcAlpha = mDot[0]();
    const volScalarField& mDotvAlpha = mDot[1]();

    mDotAlpha = mDotvAlpha*alpha1 + mDotcAlpha*(alpha2);

    if (nAlphaSubCycles > 1)
    {
        dimensionedScalar totalDeltaT = runTime.deltaT();
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(rhoPhi.dimensions(), Zero)
        );

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            #include "alphaEqn.H"
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
        }

        rhoPhi = rhoPhiSum;
    }
    else
    {
        #include "alphaEqn.H"
    }
    
    if (interfaceMethod == "isoAdvection" && !pimple.firstIter())
    {

        // Resetting U and phi to value at latest iteration.
        U = 2.0*U - U.prevIter();
        phi = 2.0*phi - phi.prevIter();
    }
    
    if (interfaceMethod == "isoAdvectionPLIC" && !pimple.firstIter())
    {

        // Resetting U and phi to value at latest iteration.
        U = 2.0*U - U.prevIter();
        phi = 2.0*phi - phi.prevIter();
    }
    
    if (interfaceMethod == "levelSet" || interfaceMethod == "levelSetScale")
    {
        #include "solveLSFunction.H"
        interfaceLS.correct();
        #include "updateFlux.H"
        rho == alpha1*rho1 + alpha2*rho2;
    }
    else if (interfaceMethod == "isoAdvector" || interfaceMethod == "MULES" || interfaceMethod == "isoAdvectorPLIC")
    {
        rho == alpha1*rho1 + alpha2*rho2;
    }
}
