from fenics import *
import random
import sys

patient = sys.argv[1]
num = int(sys.argv[2])
stokes = sys.argv[3]
assert stokes in ["true", "false"]
stokes = stokes == "true"
assert isinstance(stokes, bool)
error_rate = int(sys.argv[4])

mesh = Mesh()

meshfile = "%s/mesh/sasflow%s.h5" % (patient, num)
info("Reading mesh from %s" % meshfile)
hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
hdf.read(mesh, "/mesh", False)

print("loading")
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


print("Start loop")
max_c0 = max(c0.vector()[:])
max_c1 = max(c1.vector()[:])
max_c2 = max(c2.vector()[:])
max_c3 = max(c3.vector()[:])
max_c4 = max(c4.vector()[:])

for i, value in enumerate(c0.vector()[:]):
    c0.vector()[i] = c0.vector()[i] + random.uniform(-1,1)*error_rate*0.01*max_c0
    c1.vector()[i] = c1.vector()[i] + random.uniform(-1,1)*error_rate*0.01*max_c1
    c2.vector()[i] = c2.vector()[i] + random.uniform(-1,1)*error_rate*0.01*max_c2
    c3.vector()[i] = c3.vector()[i] + random.uniform(-1,1)*error_rate*0.01*max_c3
    c4.vector()[i] = c4.vector()[i] + random.uniform(-1,1)*error_rate*0.01*max_c4

print("Start save")
if stokes:
    save_file = File("%s/forward_imgs/stokes_vel_noise_%s_img_%s/solution.pvd" % (patient, error_rate, num))
else:
    save_file = File("%s/forward_imgs/z_vel_noise_%s_img_%s/solution.pvd" % (patient, error_rate, num))
save_file << (c0, 0)
save_file << (c1, 42)
save_file << (c2, 83)
save_file << (c3, 125)
save_file << (c4, 165)

if stokes:
    xdmf = XDMFFile(MPI.comm_world, "%s/forward_imgs/stokes_vel_noise_%s_img_%s/solution.xdmf" % (patient, error_rate, num))
else:
    xdmf = XDMFFile(MPI.comm_world, "%s/forward_imgs/z_vel_noise_%s_img_%s/solution.xdmf" % (patient, error_rate, num))
xdmf.write_checkpoint(c0, "concentration", 0, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c1, "concentration", 42, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c2, "concentration", 83, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c3, "concentration", 125, XDMFFile.Encoding.HDF5, True)
xdmf.write_checkpoint(c4, "concentration", 165, XDMFFile.Encoding.HDF5, True)