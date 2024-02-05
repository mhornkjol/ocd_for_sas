from fenics import *
import sys

patient = sys.argv[1]
num = sys.argv[2]
img = sys.argv[3]

mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "%s/mesh/sasflow%s.h5" % (patient, num), "r")
hdf.read(mesh, "/mesh", False)
hdf.close()

C = FunctionSpace(mesh, "CG", 1)
c0 = Function(C)
c1 = Function(C)
c2 = Function(C)
c3 = Function(C)
c4 = Function(C)

c_data = Function(C)


xdmf =  XDMFFile(MPI.comm_world, "%s/forward_imgs/%s_img_%s/solution.xdmf" % (patient, img, num))

xdmf.read_checkpoint(c0, "concentration", 0)
xdmf.read_checkpoint(c1, "concentration", 42)
xdmf.read_checkpoint(c2, "concentration", 83)
xdmf.read_checkpoint(c3, "concentration", 125)
xdmf.read_checkpoint(c4, "concentration", 165)
xdmf.close()


hdf = HDF5File(mesh.mpi_comm(), "%s/data_concentrations/res_%s/c1.h5" % (patient, num), "r")
hdf.read(c_data, "/function")

c_data_int = assemble(c_data*dx)

c0.vector()[:] = c0.vector()*c_data_int/assemble(c4*dx)
c1.vector()[:] = c1.vector()*c_data_int/assemble(c4*dx)
c2.vector()[:] = c2.vector()*c_data_int/assemble(c4*dx)
c3.vector()[:] = c3.vector()*c_data_int/assemble(c4*dx)
c4.vector()[:] = c4.vector()*c_data_int/assemble(c4*dx)



xdmf = XDMFFile(MPI.comm_world, "%s/forward_imgs_scaled/%s_img_%s/solution.xdmf" % (patient, img, num))

xdmf.write_checkpoint(c0, "concentration", 0, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c1, "concentration", 42, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c2, "concentration", 83, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c3, "concentration", 125, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c4, "concentration", 165, XDMFFile.Encoding.HDF5, True)


save_file = File("%s/forward_imgs_scaled/%s_img_%s/solution.pvd" % (patient, img, num))
save_file << (c0, 0)
save_file << (c1, 42)
save_file << (c2, 83)
save_file << (c3, 125)
save_file << (c4, 165)