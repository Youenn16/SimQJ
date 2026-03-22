# SimQJ

> SimQJ is a versatile toolbox designed for the simulation of Quantum Jump Trajectories, with a current focus on efficiently handling systems of free fermions.

![Build Status](https://img.shields.io/badge/build-passing-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue) 
 

## Table of Contents

- [Overview](#overview)
- [Features](#features) 
  - [Current Features](#current-features)
  - [Features to support](#features-to-support)
  - [Supported Models](#supported-models)
- [Getting Started](#getting-started)
  - [Setup](#setup)
  - [Use](#use)  
- [Project Structure](#project-structure)
- [Add a model](#add-a-model)
- [Contributing](#contributing) 
- [License](#license) 


## Overview

This package is developped to provide tools in C++ and Lapack library to simulate quantum jumps trajectories. In particular, it focuses on the case of free fermions Hamiltonian, and jumps preserving the Gaussianity. Moreover, this code allows to generate histograms for entanglement statistics (see https://link.aps.org/pdf/10.1103/PRXQuantum.5.030329). In this context, the simulation is done by default with high order integration scheme of the Stoschastic Schrodinger Equation (see https://arxiv.org/pdf/1405.6694). This creates some overhead, but can still simulate - in an acceptable time - 1D chains of 128-512 sites, without particle conservation depending on statistics wanted. 
 
 
## Features

### Current Features

- ✅ Simulate Ising Chain with density monitoring (QJ)
- ✅ BenchMark Ising Chain Monitoring with exact diagonalisation in low dimension
- ✅ Core features for simulation without particle conservation using the correlation matrix
- ✅ Density jumps 
 
 

### Features to support

- ❌ Losses jumps
- ❌ Hole density jumps
- ❌ Particle conservation support (wavefunction representation and correlation matrix) + associate method like wavefunction overlap
- ❌ Wavefunction representation for non particle conservation
- ❌ First order integration scheme
- ❌ Tests for core features
- ❌ Benchmark of new models
- ❌ PBC 


### Supported Models

- ✅ Ising Chain with density monitoring (QJ)
- ❌ Ising Chain with losses jumps (QJ)
- ❌ SSH Chain with density monitoring
- ❌ Free Fermion Chain with density monitoring
 
 



## 🚀 Getting Started

### Setup

The project is given with a bash script `setup.sh` to install necessary packages. First clone the repo


```bash
    git clone https://github.com/Youenn16/SimQJ.git
```

And run the install script
 

```bash
    bash setup.sh
```

This script installs necessary C++ packages, such as Openblas and Lapack and compile the project through Cmake in the `build` folder. On the python side, it builds helpfull commands. Then do 

```bash
    source toolbox-binary/.venv/bin/activate
```

to activate the python environment. 

---
### Use
 
#### Simulation 

To simulate a model you need to specify in `config/config.param` the configuration of the simulation, the file is structured as follow :

```
ModelParameters:
name=isingChain_density
gamma=0.3
h=0.7
kappa=1.0
J=1.0

SimulationParameters:
L=8
dt=0.005
time_max=10
no_click=False
high_order=True
record_probas=True

DataMonitoring:
time
norm
entanglement(2,4)
density(0,1)
```

The `ModelParameters` specifies the name of the model and the related parameters, `SimulationParameters` concerns some simulation options, in particular :
- `L` is the size of the chain [*mandatory*]
- `dt` is time step used in the simulation [*mandatory*]
- `time_max` is the time until which the trajectory simulation runs [*mandatory*]
- `no_click` runs only the deterministic non-Hermitian evolution [*not mandatory*, default=False]
- `high_order` choose the simulation methods between high_order and first order(not supported yet) [*not mandatory* default=True]
- `record_probas` allows to record the random numbers used to generate the trajectory, useful to compare with other simulations or rerun the simulation without a seed  [*not mandatory* default=False]
- `init_state` choice of the init state, it can be chosen between `Neel` state, `Random` in that case need to specify the initial number or particle by `nb_part`. Or it can be a precise pattern such as `0101011110` matching the full system size [*not mandatory* default=Neel]

`DataMonitoring` specifies the quantities one wants to save during the trajectories, hence 
- `time` save time points (each time step + jumps times in high order)
- `norm` save time evolution 
- `entanglement` save the bi-partite von neumann entanglement for all the arguments $\ell$ between brackets, where this area starts at site 1 and finish at site $\ell$ . 
- `density` save the on site density of the specified sites between brackets
- `entanglement_statistics` save the entanglement statistics for all the arguments $\ell$ between brackets (see https://link.aps.org/pdf/10.1103/PRXQuantum.5.030329)
 

Then the simulation can be run in `build` through 


```bash
    ./simQj
```

and the results will be recorded in `build/res`. Here are the parameters the command can support 

```bash
  --version 
  -v, --verbose                 Enable verbose output
  -o, --output <folder>         Output folder
  -b, --batch_id <n>            Batch Id 
  -l, --loop_number <n>         Number of iterations
  -c, --config_path <file>      Config path
  -h, --help                    Show this help
```

The `<folder>` and `<file>` are specified from project folder, `<n>` are integer arguments and `batch_id` is used to number the batch and `loop_number` to specify the number of trajectory runs. For example on a slurm cluster we can sample many trajectories by doing 

```bash
#!/bin/bash
#SBATCH --job-name=my_array_job     # Job name 
#SBATCH --ntasks=1                    # One task per array job
#SBATCH --cpus-per-task=2             # CPUs per task 
#SBATCH --array=1-100                 # Array range: 100 jobs (indexed 1-100)
 

# Runs 400 trajectories
./simQJ -b ${SLURM_ARRAY_TASK_ID} -l 4 
```

#### Python Scripts 

**a - Benchmark**

By doing 
```bash
simQJ benchmark --model IsingModel
```
(only `IsingModel` is supported for the moment) you can compare the C++ simulation (that can be used for large scale simulation) and an exact diagonalisation simulation in python (up to 8-10 sites). The command plots entanglement evolution for both. 

**b - Plot**

By doing 
```bash
simQJ plot -q quantity -n <number>
```
one can plot a quantity in function of time that has been saved (works only in case of one trajectory generated). For example, if one wants to plot entanglement for a subsystem of size 4, we do `simQJ plot -q entanglement -n 2`.

## Project Structure

The project has the following structure :

```
SimQJ/
├── build/
│   ├── simQJ
│   └── res
├── config/
│   └── config.param
├── include/
│   ├── models/
│   └── *.h
├── src/
│   ├── models/
│   └── *.cpp 
├── tests/
├── include/
│   └── cli.py
├── toolbox/
│   ├── exact_QJ.py
│   └── visualize_traj.py
├── toolbox-binary/
├── setup.sh
├── pyproject.toml
├── CMakeLists.txt
└── README.md
└── version.h.in

```
 
The C++ code lies in `src` (for the `.cpp`) and `include` (for the headers `.h`). 


## Add a model
TODO
 
## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a new branch: `git checkout -b feature/your-feature`
3. Make your changes and commit: `git commit -m 'Add your feature'`
4. Push to your branch: `git push origin feature/your-feature`
5. Open a Pull Request
 

---
 
 
## License

Distributed under the MIT License.   