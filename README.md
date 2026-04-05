# Advanced Programming Techniques – Octree-Based Simulation Framework in C++23

This repository contains my work for the course **Advanced Programming Techniques**.  
The project was developed as a modular **C++23 library** built with **CMake**, focusing on octree-based spatial data structures, numerical field abstractions, PDE solving, and fluid simulation. The implementation combines reusable header/source components (`.hpp` / `.cpp`) with scientific computing functionality and VTK-based visualization support.

The first milestone established the core octree and grid infrastructure together with a PDE solver, while the second milestone extended the framework with multidimensional field containers and a working **Lattice Boltzmann Method (LBM)** simulation.

---

## Project Overview

The goal of this project was to develop a reusable simulation framework in modern C++ that combines:

- octree-based spatial data structures
- structured cell-grid abstractions
- numerical PDE solving
- multidimensional field storage
- scientific visualization with VTK
- transient fluid simulation using LBM

The work progressed from low-level infrastructure and traversal utilities to full simulation pipelines with exportable results for visualization in ParaView.

---

## Milestone 1 – Octree Data Structures and PDE Solver

The first milestone focused on building the foundation of the framework.

### Main work completed
- Set up a **C++23** project with:
  - CMake
  - CTest
  - clang-format
  - clang-tidy
- Implemented core support classes such as:
  - `Vec` for fixed-size algebraic vectors
  - `MortonIndex` for octree path encoding and traversal
  - `Box` and `OctreeGeometry` for geometric operations
- Built the main `CellOctree` data structure with:
  - node-stream storage
  - descriptor-based tree construction
  - cell views
- Added **VTK HyperTreeGrid export** for octree visualization in ParaView
- Implemented traversal infrastructure:
  - cursors
  - iterators
  - ranges
  - Morton-order depth-first traversal
  - horizontal per-level traversal
- Developed the `CellGrid` abstraction with:
  - cell enumeration
  - adjacency lists
  - periodicity handling
  - range-based access
- Implemented a **Jacobi solver** for the **3D Poisson equation**
- Exported numerical solution, right-hand side, and residual fields to VTK

### Applications
- **create-htgfile**  
  Creates and exports octrees from descriptor strings
- **poisson**  
  Solves a 3D Poisson problem on a uniformly refined octree grid and exports the results

---

## Milestone 2 – Numerical Fields and Fluid Simulation with LBM

The second milestone extended the framework toward transient and field-based numerical simulation.

### Main work completed
- Implemented `GridVector<T, Q>` for scalar and vector-valued numerical fields on a `CellGrid`
- Added proper **rule-of-five** semantics for memory-safe multidimensional field storage
- Developed `mdspan`-based non-owning views for efficient structured access
- Extended the VTK export pipeline to support:
  - scalar fields
  - vector-valued fields
  - multidimensional simulation output
- Implemented the core kernels of the **Lattice Boltzmann Method (LBM)** using the **D3Q19 stencil**
- Added kernels for:
  - PDF initialization
  - macroscopic density computation
  - velocity computation
  - collision
  - streaming
- Built a **Taylor–Green vortex decay** simulation on a fully periodic uniform grid
- Compared the numerical solution against the analytical benchmark
- Exported density and velocity fields for visualization in ParaView

### Applications
- **tgv**  
  Runs a Taylor–Green vortex decay simulation using LBM and exports the final fields

---

## Key Features

- Modern **C++23** design
- Octree-based spatial representation
- Morton-index-based traversal
- Structured cell-grid abstraction
- Custom iterators and ranges
- VTK export for ParaView visualization
- Jacobi solver for elliptic PDEs
- Multidimensional field containers
- `mdspan`-based views
- LBM fluid simulation with D3Q19
- Periodic boundary handling
- Analytical benchmark validation

---

## Skills and Concepts Applied

This project involved practical work in:

- Modern C++23
- Template programming
- Rule of five
- Builder-style design
- Tree-based spatial data structures
- Scientific computing
- PDE discretization
- Iterative solvers
- Fluid simulation
- Lattice Boltzmann Method
- Multidimensional data layouts
- Scientific visualization with VTK

---

## Summary

This project demonstrates the development of a complete **C++23 scientific simulation framework**, starting from octree and grid infrastructure and extending to numerical PDE solving and transient fluid simulation. It combines software engineering, numerical methods, and scientific visualization in a structured and reusable implementation.