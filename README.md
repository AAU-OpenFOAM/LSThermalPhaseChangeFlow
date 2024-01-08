# LSThermalPhaseChangeFlow Solver

## Introduction
`LSThermalPhaseChangeFlow` is an OpenFOAM solver developed for simulating thermal phase change processes. This solver, compatible with OpenFOAM v2006, integrates both VOF-isoAdvector and coupled Level Set with VOF (CLSVOF) methods, providing versatility in interface capturing.

## Features
- Utilizes the VOF-isoAdvector method for sharp interface capturing and reduced non-physical currents.
- Includes Level Set techniques for precise curvature prediction.
- Capable of simulating a wide range of thermal phase change phenomena.
- Provides examples of various thermal phase change benchmark cases.
- Demonstrates improvements in computational efficiency and accuracy over traditional methods.

## Installation
To install `LSThermalPhaseChangeFlow`, ensure you have OpenFOAM v2006 installed and sourced. Clone this repository into your OpenFOAM user directory and compile the solver using the provided Makefiles.

```bash
cd $FOAM_RUN
git clone https://github.com/AAU-OpenFOAM/LSThermalPhaseChangeFlow.git
cd LSThermalPhaseChangeFlow
./Allwmake
