# SDURail

SDURail extends SoRoSim with a rail model and simulation framework.

## Prerequisites

Before using SDURail, install and set up:

- MATLAB
- SoRoSim

## Installation

### 1. Clone SoRoSim

Clone the SoRoSim repository to a suitable location:

```bash
git clone https://github.com/SoRoSim/SoRoSim.git
```

### 2. Clone SDURail

Clone this repository into the Differentiable_SoRoSim directory of your SoRoSim installation:
```text
SoRoSim/
└── Differentiable_SoRoSim/
    └── SDURail/
```

## Configuration

Rail parameters can be modified in:

```bash
link/RailParameters.m
```

Adjust the values in this file to match the rail geometry and simulation setup you wish to test.

## Running a Simulation

From MATLAB:

Open the project.
Navigate to the SDURail directory.
Run:

```bash
simulation_script.m
```

## Output

The solution is plotted in the plots/ folder, and an animation called Dynamics.mp4 is created.

![SDURail Simulation](RadialRail.gif)



