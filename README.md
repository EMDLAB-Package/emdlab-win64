
# EMDLAB Finite Element Package

### Electrical Machines Design Laboratory

---

## Overview

**EMDLAB** is an open-source numerical package developed in MATLAB for the design and analysis of electrical machines, including motors, generators, transformers, and actuators.

Designed to support academic research and education, EMDLAB provides a transparent and reproducible computational environment. Its open-source architecture encourages collaboration, facilitates code extension, and makes it particularly suitable for developing and investigating new numerical methods.

The package is structured as a **library of MATLAB objects**, allowing users to select and combine the modules required for their specific applications. Furthermore, EMDLAB provides the flexibility to develop customized, standalone software tailored to individual design and analysis requirements.

---

## Setup Instructions (Windows 64-bit)

EMDLAB can be installed using either of the following methods.

### Method 1: Web-Based Installer (Recommended)

The web-based installer automatically downloads and installs the latest version of EMDLAB.

1. Download **[EMDLAB Setup (Windows 64-bit)](https://comprogexpert.com/emdlab_setup.exe)**.
2. Run `emdlab_setup.exe` and follow the installation instructions.
3. To update EMDLAB, simply run the same installer again. It automatically downloads and installs the latest version.

### Method 2: Manual Installation

1. Download the **[`emdlab-win64.zip`](https://github.com/EMDLAB-Package/emdlab-win64/archive/refs/heads/main.zip)** archive.
2. Extract the ZIP file.
3. Copy the extracted **`emdlab-win64`** folder directly into the `C:\` directory.
4. Make sure the folder name remains unchanged.

---

## MATLAB Configuration

To use EMDLAB in your MATLAB projects, add the following command at the beginning of your MATLAB script (`.m` file):

```matlab
addpath(genpath('C:\emdlab-win64'));
```

This command automatically adds the EMDLAB directory and all its subdirectories to the MATLAB search path.

---

## Tutorials and Academic References

### Installation Tutorial

For a step-by-step guide to installing EMDLAB from GitHub, watch the following video:

**[Watch the Installation Tutorial on YouTube](https://youtu.be/ifwybm4r2_0)**

### Main EMDLAB Publication

For academic use, please cite the following publication:

**[Read the EMDLAB paper](https://www.sciencedirect.com/science/article/pii/S2352711025004121)**

### Second-Order Magnetostatic Solver

The following paper presents the second-order magnetostatic solver developed for EMDLAB:

**[Read the paper on Springer](https://link.springer.com/article/10.1007/s00366-026-02336-y)**

---

## Using Gmsh as an External Mesh Generator

In addition to its built-in mesh generators, EMDLAB supports **Gmsh** as an external mesh-generation tool.

Follow these steps to configure and use Gmsh.

### Step 1. Check Python Availability

In the MATLAB Command Window, execute:

```matlab
pyenv
```

This command displays the Python environment configured for MATLAB. Make sure that Python is accessible.

### Step 2. Install Gmsh

Open the Windows Command Prompt (`cmd`) and execute:

```bash
pip install gmsh
```

### Step 3. Generate Meshes Using Gmsh

When calling the `.generateMesh` method, specify `'gmsh'` as the input argument:

```matlab
.generateMesh('gmsh')
```

EMDLAB will then use Gmsh as the external mesh generator.

---

**EMDLAB | Electrical Machines Design Laboratory**
