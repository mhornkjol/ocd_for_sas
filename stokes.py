from fenics import *
import sys

patient = sys.argv[1]
num = int(sys.argv[2])

parameters["krylov_solver"]["monitor_convergence"] = True

mu = 7e-4

meshfile = "%s/mesh/sasflow%s.h5" % (patient, num) # sys.argv[1]
outfolder = "%s/stokes_vel_field/stokes_vel_field_%s" % (patient, num) # sys.argv[2]

print("Loading mesh")
mesh = Mesh()
hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
hdf.read(mesh, "/mesh", False)

subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
hdf.read(subdomains, "/subdomains")
boundaries  = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
hdf.read(boundaries, "/boundaries")


# Build function space
print("Building function space")
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 1)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)

for facet in facets(mesh):
    if boundaries[facet.index()] == 6:
        for vertex in vertices(facet):
            if vertex.point().array()[2] > -17:
                boundaries[facet.index()] = 8

boundary_file = File('%s/boundaries.pvd' % outfolder)
boundary_file << boundaries


# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)

n = FacetNormal(mesh)

# Boundary conditions
print("Setting boundary conditions")
no_slip = Constant((0.0, 0.0, 0.0)) # Noslip

bc1 = DirichletBC(W.sub(0), no_slip, boundaries, 5)
bc2 = DirichletBC(W.sub(0), no_slip, boundaries, 6)

bcs = [bc1, bc2]

dx = Measure("dx", domain=mesh, subdomain_data=subdomains)
ds = Measure("ds", domain=mesh, subdomain_data=boundaries)

print("Defining equations")
f = Constant((0.0, 0.0, 0.0))
pres = Constant(8e-8)
h = mesh.hmin()
beta = 0.01
a = mu*inner(grad(u), grad(v))*dx - div(v)*p*dx - q*div(u)*dx #- inner(grad(p),grad(q))*h**2*beta*dxF
p = mu*inner(grad(u), grad(v))*dx + (1.0/mu)*p*q*dx #+ inner(grad(p),grad(q))*h**2*beta*dxF

L = inner(f, v)*dx - inner(pres*n, v)*ds(7)

print("Assembeling")
A = assemble(a, keep_diagonal=True)
P = assemble(p, keep_diagonal=True)
b = assemble(L)
[bc.apply(P,b) for bc in bcs]
[bc.apply(A,b) for bc in bcs]

print("Solving system")
U = Function(W)
solver = KrylovSolver("minres", "amg")
solver.set_operators(A, P)
it = solver.solve(U.vector(), b)
u, p = U.split(deepcopy=True)

print("Saving solution")
ufile_pvd = File("%s/velocity.pvd" % outfolder)
pfile_pvd = File("%s/pressure.pvd" % outfolder)
ufile_pvd << u
pfile_pvd << p

# Save velocity
# Save velocity
save_file_velocity = "%s/velocity.xdmf" % outfolder 
save_file_pressure = "%s/pressure.xdmf" % outfolder 
with XDMFFile(MPI.comm_world, save_file_velocity) as xdmf:
    xdmf.write_checkpoint(u, "velocity", 0)

with XDMFFile(MPI.comm_world, save_file_pressure) as xdmf:
    xdmf.write_checkpoint(p, "pressure", 0)