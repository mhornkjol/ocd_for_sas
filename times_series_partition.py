from fenics import *
import sys

patient = sys.argv[1]
num = int(sys.argv[2])
stokes = sys.argv[3]
assert stokes in ["true", "false"]
stokes = stokes == "true"
assert isinstance(stokes, bool)

mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "%s/mesh/sasflow%s.h5" %(patient, num), "r")
hdf.read(mesh, "/mesh", False)

C = FunctionSpace(mesh, "CG", 1)
c0 = Function(C)
c1 = Function(C)
c2 = Function(C)
c3 = Function(C)
c4 = Function(C)

if stokes:
    xdmf =  XDMFFile(MPI.comm_world, "%s/forward_imgs/stokes_vel_img_%s/solution.xdmf" % (patient, num))
else:
    xdmf =  XDMFFile(MPI.comm_world, "%s/forward_imgs/z_vel_img_%s/solution.xdmf" % (patient, num))

xdmf.read_checkpoint(c0, "concentration", 0)
xdmf.read_checkpoint(c1, "concentration", 42)
xdmf.read_checkpoint(c2, "concentration", 83)
xdmf.read_checkpoint(c3, "concentration", 125)
xdmf.read_checkpoint(c4, "concentration", 165)

xdmf.close()

if stokes:
    xdmf =  XDMFFile(MPI.comm_world, "%s/forward_imgs/stokes_vel_img_%s_mini/solution.xdmf" % (patient, num))
    save_file = File("%s/forward_imgs/stokes_vel_img_%s_mini/solution.pvd" % (patient, num))
else:
    xdmf =  XDMFFile(MPI.comm_world, "%s/forward_imgs/z_vel_img_%s_mini/solution.xdmf" % (patient, num))
    save_file = File("%s/forward_imgs/z_vel_img_%s_mini/solution.pvd" % (patient, num))

xdmf.write_checkpoint(c0, "concentration", 0, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c1, "concentration", 42, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c2, "concentration", 83, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c3, "concentration", 125, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c4, "concentration", 165, XDMFFile.Encoding.HDF5, True)

save_file << (c0, 0)
save_file << (c1, 42)
save_file << (c2, 83)
save_file << (c3, 125)
save_file << (c4, 165)