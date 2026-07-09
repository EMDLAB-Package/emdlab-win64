# run_geo.py
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("Mesh.Smoothing", 10)
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.Algorithm", 8)

# Load .geo file
gmsh.open("C:\\emdlab-win64\\tmp\\emdlab_gmsh_geoFile.geo") 

# Generate 2D mesh
gmsh.model.mesh.generate(2)
gmsh.model.mesh.recombine()

# Write mesh to file
gmsh.write("C:\\emdlab-win64\\tmp\\emdlab_gmsh_mshFile.m")

gmsh.finalize()
