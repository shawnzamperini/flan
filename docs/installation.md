<link rel="stylesheet" href="assets/style.css">

<div class="container">
<div class="sidebar-container">
{% include_relative _sidebar.md %}
</div>

<div class="content-container">

# Installation

FLAN is built as a modern C++20 project using CMake and supports both serial
and hybrid MPI/OpenMP execution. This page walks you through installing and
building FLAN on laptops, workstations, and HPC systems such as Perlmutter.

---

## Prerequisites

Before building FLAN, ensure the following dependencies are available:

### **Required**
- C++20-capable compiler  
  - GCC ≥ 11  
  - Clang ≥ 14  
  - NVHPC ≥ 23  
- CMake ≥ 3.20  
- OpenMP support  
- MPI implementation (for parallel runs)  
  - MPICH, OpenMPI, Cray-MPICH, etc.

### **Optional**
- **HDF5** — for background plasma profiles  
- **OpenADAS rate tables** — for ionization/recombination physics  
- **Python + NumPy** — for post-processing and analysis  

---

## Clone the Repository

</div>
</div>
