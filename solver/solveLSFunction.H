// solve Level-Set function as the re-initialization equation
Info<< "solve the reinitialization equation"<< nl << endl;
       
interfaceLS.correctPsi0();
psi == interfaceLS.psi0();

for
(
    int corr=0;
    corr<int(interfaceLS.epsilon().value()/interfaceLS.deltaTau().value());
    corr++
)
{
    psi = psi
  + interfaceLS.psi0()/mag(interfaceLS.psi0())*(double(1)
  - mag(fvc::grad(psi)*interfaceLS.dimChange()))*interfaceLS.deltaTau();
    psi.correctBoundaryConditions();
}
   
if (interfaceMethod == "levelSet")
{
    interfaceLS.correctDeltaH();
}
else if (interfaceMethod == "levelSetScale")
{
    interfaceLS.correctDeltaHscale();
}

   
   
