### EMDLAB finite element package

**Electrical Machines Design Laboratory**

**Overview:**
EMDLAB is an open-source numerical package developed in the MATLAB environment for the design
and analysis of electrical machines, including motors, generators, transformers, and actuators. 
The package provides a transparent and reproducible tool for academic research and education. 
Its open-source nature encourages collaboration and facilitates easy extension of the code, making it particularly valuable for method development and academic use.

The package is organized as a library of MATLAB objects. Depending on your specific design needs, 
you can select the appropriate modules to obtain your desired results. With EMDLAB, 
you can also develop customized, standalone software tailored to your applications.

### Setup Instructions (Windows 64-bit):

✅ 1) Download the emdlab-win64.zip file.

✅ 2) Extract the zip file and place the "emdlab-win64" folder in the "C:\" directory, without changing the folder name.

✅ 3) To use the EMDLAB package in your MATLAB code, add the following line at the beginning of your mfile:

***--->> addpath(genpath('C:\emdlab-win64'));***

**How to install EMDLAB from GitHub? (Follow the video link):**
https://youtu.be/ifwybm4r2_0

**For academic use, cite this paper:**
https://www.sciencedirect.com/science/article/pii/S2352711025004121

**This paper presents the second-order magnetic-static solver of EMDLAB:**
https://link.springer.com/article/10.1007/s00366-026-02336-y

### How to use Gmsh as a mesh generator in EMDLAB:
In addition to the built-in EMDLAB mesh generators, it is also possible to use Gmsh as an external mesh generator.

✅ 1) **Install Python:** In the MATLAB command window, type `>> pyenv` to check whether Python is accessible.

✅ 2) **Install Gmsh** via cmd window: `pip install gmsh`

✅ 3) When using the .generateMesh method, use 'gmsh' as input: **.generateMesh('gmsh')**
