# Installation

Flan is built as a modern C++23 project using CMake and supports both serial
and hybrid MPI/OpenMP execution. This page walks you through installing and
building Flan on laptops, workstations, and HPC systems such as Perlmutter.

---

## Prerequisites

Flan has the following dependencies:

- C++23-capable compiler  
- CMake ≥ 3.26
- OpenMP support  
- MPI  
  - MPICH, OpenMPI, Cray-MPICH, etc.
- HDF5
- NetCDF
- zlib
- conda (or mamba)

HDF5, NetCDF and zlib are installed by machine-specific scripts to install
dependencies, found in `machines/`. Other dependencies must be provided by the
computing environment.

---

## Clone the Repository

Cloning the repository is straightforward. This will clone it into `flan/`:

```bash
git clone https://github.com/shawnzamperini/flan.git
```

---

## Machine specific instructions

Select your machine from the following list. Use 
[Linux](install-linux.md) for a general purpose installation,
e.g., a personal laptop.

- [Linux](install-linux.md)
- [Perlmutter](install-perlmutter.md)
