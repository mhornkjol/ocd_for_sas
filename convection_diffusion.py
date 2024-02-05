from fenics import *
import numpy as np
import sys

patient = sys.argv[1]
num = int(sys.argv[2])
stokes = sys.argv[3]
assert stokes in ["true", "false"]
stokes = stokes == "true"
assert isinstance(stokes, bool)

D = 3.8e-4*3600
T = 6 #7*24*3600 # final time
num_steps = 6000 # number of time steps
tau = T/num_steps

mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), "%s/mesh/sasflow%s.h5" % (patient, num), "r")
hdf.read(mesh, "/mesh", False)
print(f'mesh num vertices {mesh.num_vertices()}')
print(f'mesh num cells {mesh.num_cells()}')
print(f'h max {mesh.hmax()}')
print(f'h min {mesh.hmin()}')
subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
hdf.read(subdomains, "/subdomains")
boundaries  = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
hdf.read(boundaries, "/boundaries")

# Build vectorfunction and functionspace
print("Building function space")
V = FunctionSpace(mesh, 'CG', 1)
O = VectorFunctionSpace(mesh, 'CG', 1)

# Define initial value
c_ = Function(V)

# Define variational problem
c = TrialFunction(V)
d = TestFunction(V)

dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
phi = Function(O)

if stokes:
    xdmf = XDMFFile(MPI.comm_world, "%s/stokes_vel_fields/stokes_vel_field_%s/velocity.xdmf" % (patient, num))
    xdmf.read_checkpoint(phi, "velocity", 0)

    # Scale stokes to 20 mm/h max velocity magnitude
    u, v, w = phi.split(deepcopy=True)
    magnitude = np.sqrt(u.vector()[:]**2 + v.vector()[:]**2 + w.vector()[:]**2)
    max_vel = 20
    factor = max_vel/max(magnitude)
    phi.vector()[:] = phi.vector()[:]*factor
else:
    z_vector = Expression((0,0,'20'), element=O.ufl_element())
    phi.interpolate(z_vector)


F = (1.0/tau*(c - c_)*d - inner(c*phi,grad(d)) + inner(D*grad(c), grad(d)))*dx()
a, L = system(F)
bc = DirichletBC(V, 1, boundaries, 7)


boundary_value_file = File("boundary_values.pvd")
m = Function(V)
bc.apply(m.vector())
boundary_value_file << m

if stokes:
    save_file = File("%s/forward_imgs/stokes_vel_img_%s/solution.pvd" % (patient, num))
    xdmf = XDMFFile(MPI.comm_world, "%s/forward_imgs/stokes_vel_img_%s/solution.xdmf" % (patient, num))
else:
    save_file = File("%s/forward_imgs/z_vel_img_%s/solution.pvd" % (patient, num))
    xdmf = XDMFFile(MPI.comm_world, "%s/forward_imgs/z_vel_img_%s/solution.xdmf" % (patient, num))

c = Function(V)
t = 0
times = []
for i in range(num_steps):
    print(f'{i}')
    solve(a == L, c, bc, solver_parameters={"linear_solver": "mumps"})

    c_.assign(c)

    # if t>=5 and t <= 6:
    if i == 5000 or i == 5250 or i == 5500 or i == 5750 or i == 5999:
        save_file << (c, t)
        xdmf.write_checkpoint(c, "concentration", len(times), XDMFFile.Encoding.HDF5, True)
        times.append(t)

    t += tau

xdmf.close()
