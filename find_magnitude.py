from fenics import *
import numpy as np
import sys

meshfile = sys.argv[1]
velfile = sys.argv[2]

mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
hdf.read(mesh, "/mesh", False)
hdf.close()

V = VectorFunctionSpace(mesh, "CG", 1)
phi = Function(V)

hdf = HDF5File(mesh.mpi_comm(), velfile, "r")
hdf.read(phi, "/function")

u, v, w = phi.split(deepcopy=True)

magnitude = np.sqrt(u.vector()[:]**2 + v.vector()[:]**2 + w.vector()[:]**2)

print("The maximum magnitude is: ", max(magnitude))
