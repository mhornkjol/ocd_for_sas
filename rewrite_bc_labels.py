from fenics import *
import sys
import os

mesh = Mesh()

meshfolder = sys.argv[1]
meshname = sys.argv[2]
outfolder = sys.argv[3]

hdf = HDF5File(mesh.mpi_comm(), os.path.join(meshfolder, meshname), "r")
hdf.read(mesh, "/mesh", False)

boundaries  = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
hdf.read(boundaries, "/boundaries")
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
hdf.read(subdomains, "/subdomains")
hdf.close()

for facet in facets(mesh):
    if boundaries[facet.index()] == 5:
        boundaries[facet.index()] = 8
    if boundaries[facet.index()] == 7:
        boundaries[facet.index()] = 5

for facet in facets(mesh):
    if boundaries[facet.index()] == 8:
        boundaries[facet.index()] = 7

hdf = HDF5File(mesh.mpi_comm(), os.path.join(outfolder, meshname), "w")
hdf.write(mesh, "/mesh")
hdf.write(boundaries, "/boundaries")
hdf.write(subdomains, "/subdomains")

boundary_file = File("%s/boundaries.pvd" % outfolder)
boundary_file << boundaries
