from fenics import *

# Load mesh, subdomains and boundaries
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "PatID-004_cut/mesh/sasflow.h5", "r")
hdf.read(mesh, "/mesh", False)
print ("mesh num vertices ", mesh.num_vertices())
print ("mesh num cells", mesh.num_cells())
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
hdf.read(subdomains, "/subdomains")
boundaries  = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
hdf.read(boundaries, "/boundaries")


V = FunctionSpace(mesh, 'CG', 1)

c = Function(V)
xdmf =  XDMFFile(MPI.comm_world, "PatID-004_cut/results/solution.xdmf")
vtkfile = File("PatID-004_cut/results/solution.pvd")
for i in range(10):
        xdmf.read_checkpoint(c,"concentration",i) 
        vtkfile << (c, i)
xdmf.close()
