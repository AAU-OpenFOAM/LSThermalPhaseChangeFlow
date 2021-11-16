fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
  - fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
  + turbulence->divDevRhoReff(rho, U)
);

UEqn.relax();

if (interfaceMethod == "levelSet")
{
    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
            ==
            fvc::reconstruct
            (
                (
                    interfaceLS.surfaceTensionForce()
                  + surfForces.accelerationForce()  
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
    }
}
else if (interfaceMethod == "isoAdvectorPLIC")
{
    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
            ==
            fvc::reconstruct
            (
                (
                    surfForces.surfaceTensionForce()
                  + surfForces.accelerationForce()
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
    }
}
else if (interfaceMethod == "MULES" || interfaceMethod == "isoAdvector") 
{
    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
            ==
            fvc::reconstruct
            (
                (
                    interface.surfaceTensionForce()
                  + surfForces.accelerationForce()
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
    }
}
else if (interfaceMethod == "levelSetScale")
{
    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
            ==
            fvc::reconstruct
            (
                (
                    interfaceLS.surfaceTensionForceScale()
                  + surfForces.accelerationForce()  
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );
    }
}


