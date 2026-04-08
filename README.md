# Assignment1_HPSC_DEM
# OpenMP Parallelisation of a 3D DEM Solver

This project implements a three-dimensional **Discrete Element Method (DEM)** solver for spherical particles and parallelises it using OpenMP. The aim is to study performance improvements and understand bottlenecks in particle-based simulations.

---

## Overview

The DEM solver simulates particle motion by computing forces between particles and updating their positions over time. Each particle follows Newton’s laws, and interactions are modelled using a simple spring-dashpot contact model.

The main challenge in DEM is that particle interactions scale as **O(N²)**, which makes simulations slow for large systems. This project focuses on identifying this bottleneck and improving performance using parallel computing.

---

## Features

* 3D particle simulation inside a box
* Spring-dashpot contact model
* Semi-implicit Euler time integration
* Multiple verification tests:

  * Free fall
  * Constant velocity
  * Bouncing particle
* Profiling to identify bottlenecks
* OpenMP-based parallel implementation
* Performance analysis (speedup and runtime scaling)

---

## File Structure

```
.
├── dem_serial.cpp        # Serial implementation
├── dem_parallel.cpp      # OpenMP parallel version
├── plots/                # Generated plots (PDFs)
├── report.tex            # LaTeX report
└── README.md
```

---

## Compilation and Execution

### Serial Version

```bash
g++ dem_serial.cpp -O2 -o dem_serial
./dem_serial
```

### Parallel Version (OpenMP)

```bash
g++ dem_parallel.cpp -O2 -fopenmp -o dem_parallel
export OMP_NUM_THREADS=2
./dem_parallel
```

---

## Parallelisation Strategy

* Loops over particles are parallelised using OpenMP
* The particle contact loop is the main bottleneck
* Race conditions arise because multiple threads update shared forces
* This is handled using `#pragma omp atomic`

---

## Results

* The particle contact loop accounts for **95–99%** of runtime
* Speedup is limited for small systems due to:

  * Atomic operation overhead
  * Thread management overhead
* Performance improves as the number of particles increases

---

## Plots Included

* Free fall verification
* Constant velocity
* Bouncing motion
* Convergence study
* Kinetic energy vs time
* Speedup vs particle count
* Runtime comparison

---

## Limitations

* O(N²) scaling limits performance
* Atomic operations reduce parallel efficiency
* Not suitable for very large particle systems

---

## Possible Improvements

* Use neighbour lists (Verlet lists)
* Use thread-local force accumulation
* Domain decomposition for large-scale simulations

---

## Acknowledgements

Some help from large language models (ChatGPT and Claude) was used for:

* Debugging C++ code
* Fixing LaTeX issues
* Generating plots

All concepts and implementation were understood and verified independently.

---

## Author

**Anushka Chauhan**
IIT Mandi

