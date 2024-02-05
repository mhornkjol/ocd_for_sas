import sys
import os
import glob
import shutil
import resource
import numpy
import math
import argparse
import time
import datetime
import csv

from mpi4py import MPI

from dolfin import *

from compute_ocdr import diffusion
from postprocess_phi import csvwrite, compute_magnitude_field

import matplotlib.pyplot as plt


from dolfin_adjoint import *

def bar(phi, volume=None):
    phi_mag = compute_magnitude_field(phi)
    if not volume:
        mesh = phi.function_space().mesh()
        volume = assemble(1.0*dx(domain=mesh))
    avg_phi_mag = assemble(phi_mag*dx())/volume
    return avg_phi_mag
    
def compute_ocd_reduced(c0, c1, tau, D, alpha, n, results_dir, boundaries, space="CG", reg="H1"):
    # c0:
    # c1:
    # tau:   time step
    # D:     diffusion coefficient
    # alpha: regularization parameter
    # space: finite element space for the velocity field ("CG" | "RT" | "BDM")
    # reg:   if "H1", use H1-regularization, else use H(div)
    
    info("Computing OCD via reduced approach")

    # Define mesh and function space for the concentration
    mesh = c0.function_space().mesh()
    C = FunctionSpace(mesh, "CG", 1)

    # Space for the convective velocity field phi
    if space == "CG":
        Q = VectorFunctionSpace(mesh, "CG", 1)
    else:
        Q = FunctionSpace(mesh, space, 1)
    phi = Function(Q, name="Control")
    # phi_val = Expression((0, 0, '1'), element=Q.ufl_element())
    # phi.interpolate(phi_val)

    # Regularization term
    def R(phi, alpha, mesh):
        if reg == "H1":
            form = 0.5*alpha*(inner(grad(phi), grad(phi)))*dx(domain=mesh) # 0.5*alpha*(inner(phi, phi) + inner(grad(phi), grad(phi)))*dx(domain=mesh)
        else:
            form = 0.5*alpha*(inner(div(phi), div(phi)))*dx(domain=mesh) # 0.5*alpha*(inner(phi, phi) + inner(div(phi), div(phi)))*dx(domain=mesh)
        return form

    # Define previous solution
    c_ = Function(C)
    c_.assign(c0) # Hack to make dolfin-adjoint happy, maybe just start tape here?

    c2 = Function(C)
    c2.assign(c1)
    
    # Define variational problem
    c = TrialFunction(C)
    d = TestFunction(C)

    F = (1.0/tau*(c - c_)*d + div(c*phi)*d)*dx() + diffusion(D, c, d)
    a, L = system(F)
    
    bc = DirichletBC(C, c2, boundaries, 7)


    name = lambda s: os.path.join(results_dir, "pvd", s)
    boundary_file = File(name("boundary_file.pvd"))
    boundary_file << boundaries
    boundary_value_file = File(name("boundary_values.pvd"))
    m = Function(C)
    bc.apply(m.vector())
    boundary_value_file << m

    # ... and solve it once
    c = Function(C, name="State")
    solve(a == L, c, bc, solver_parameters={"linear_solver": "mumps"})

    # Output max values of target and current solution for progress
    # purposes
    info("\max c_0 = %f" % c_.vector().max())
    info("\max c_1 = %f" % c2.vector().max())
    info("\max c = %f" % c.vector().max())

    # Define the objective functional
    j = 0.5*(c - c2)**2*dx(domain=mesh) + R(phi, alpha, mesh)
    J = assemble(j)
    j0 = 0.5*(c - c2)**2*dx(domain=mesh)
    jr = R(phi, alpha, mesh)
    J0 = assemble(j0)
    Jr = assemble(jr)
    info("J (initial) = %f" % J)
    info("J0 (initial) = %f" % J0)
    info("Jr (initial) = %f" % Jr)
    
    # Define control field
    m = Control(phi)

    # Define call-back for output at each iteration of the optimization algorithm
    name = lambda s: os.path.join(results_dir, s)   
    dirname = results_dir

    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    header = ("j", "\max \phi")
    js = []

    def eval_cb(j, phi):
        jr = assemble(R(phi, alpha, mesh))
        j0 = j-jr
        values = (j, jr, j0, phi.vector().max())
        mem = resource.getrusage(resource.RUSAGE_SELF)[2]
        info("Current memory usage: %g (MB)" % (mem/1024))
        info("\tj = %f, \max phi = %f (mm/h)" % (values[0],values[3]))
        csvwrite(name("optimization_values.csv"), values, header, mode="a")

        # Read the optimization counter file, update counter, and
        # write it back, geez. 
        counter = 0
        with open(os.path.join(results_dir, "counter.csv"), "r") as f:
            reader = csv.reader(f)
            for row in reader:
                counter = int(row[0])
        counter += 1

        with open(os.path.join(results_dir, "counter.csv"), "w") as f:
            info("Updating counter file, counter is now %d " % counter)
            writer = csv.writer(f)
            writer.writerow((counter,))

        # Write current control variable to file in HDF5 and PVD formats
        # file = HDF5File(mesh.mpi_comm(), name("opt_phi_%d.h5" % counter), "w")
        # file.write(phi, "/function", 0)
        # file.close()

    # Define reduced functional in terms of J and m
    Jhat = ReducedFunctional(J, m, eval_cb_post=eval_cb)

    # Minimize functional
    tol = 1.0e-8
    phi_opt = minimize(Jhat,
                       tol=tol, 
                       options={"gtol": tol, "maxiter": 100, "disp": True})
    pause_annotation()

    # Update phi, and do a final solve to compute c
    phi.assign(phi_opt)
    solve(a == L, c, bc, solver_parameters={"linear_solver": "mumps"})

    J = assemble(j)
    j0 = 0.5*(c - c2)**2*dx(domain=mesh)
    jr = R(phi, alpha, mesh)
    J0 = assemble(j0)
    Jr = assemble(jr)
    info("J  = %f" % J)
    info("J0 = %f" % J0)
    info("Jr = %f" % Jr)
    
    return (c, phi) 
        
def map_patient(patient_dir, n, alpha, key, results_dir, normalized, imfolder, dispersion):
    """
    patient_dir: Absolute patient directory e.g. "../../241"
    n:           Mesh resolution e.g. 16 (depending on available meshes)
    alpha:       Regularization parameter
    key:         0h, 6h or 24h (starting time)
    results_dir: Name for results directory (subdirectories will be created)
    """

    # Read FEniCS mesh 
    mesh = Mesh()
    
    meshfile = os.path.join(patient_dir, "mesh", "sasflow%d.h5" % n)
    
    if not os.path.isfile(meshfile):
        raise Exception("Missing meshfile %s" % meshfile)
    
    
    info("Reading mesh from %s" % meshfile)
    hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
    hdf.read(mesh, "/mesh", False)

    info("Reading subdomains (SAS 1) from %s" % meshfile)
    subdomains = MeshFunction("size_t", mesh, mesh.topology().dim())
    hdf.read(subdomains, "/subdomains")

    info("Reading boundaries (SAS 1) from %s" % meshfile)
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    hdf.read(boundaries, "/boundaries")
    
    hdf.close()
    info("... with %d vertices and %d cells" % (mesh.num_vertices(), mesh.num_cells()))

    D = 3.8e-4*3600*dispersion

    C = FunctionSpace(mesh, "CG", 1)
    c0 = Function(C)
    c1 = Function(C)

    # Read forward images
    t0, t1 = key.split("_")
    t0 = int(t0)
    t1 = int(t1)
    info("Reading constructed solution data from %s" % imfolder)
    info(os.path.join(patient_dir, imfolder, "solution.xdmf"))
    xdmf =  XDMFFile(MPI.comm_world, os.path.join(patient_dir, imfolder, "solution.xdmf"))
    
    xdmf.read_checkpoint(c0, "concentration", t0)
    xdmf.read_checkpoint(c1, "concentration", t1)
    tau = (t1-t0)*0.006
    info("The time difference is %f" % tau)

    if normalized:
        info("Normalizing")
        c0_norm = assemble(c0*dx)
        c1_norm = assemble(c1*dx)
        c0.vector()[:] = numpy.divide(c0.vector()[:], c0_norm/50000)
        c1.vector()[:] = numpy.divide(c1.vector()[:], c1_norm/50000)

    # Store images and solution to .h5 format for postprocessing 
    name = lambda s: os.path.join(results_dir, "hdf5", s)
    info("Storing images and data to %s" % name("..."))
    file = HDF5File(mesh.mpi_comm(), name("c0.h5"), "w")
    file.write(c0, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("c1.h5"), "w")
    file.write(c1, "/function", 0)
    file.close()

    # Make an optimization counter file
    # if rank == 0:
    with open(os.path.join(results_dir, "counter.csv"), "w") as f:
        writer = csv.writer(f)
        writer.writerow((0,))     
            
    # Compute OCD approximation via reduced method
    space = "CG"
    reg = "H1"
    (c, phi) = compute_ocd_reduced(c0, c1, tau, D, alpha, n, results_dir, boundaries,
                                   space=space, reg=reg)

    # -----------------------------------------------------------
    # The following is boiler-plate post-processing code
    # -----------------------------------------------------------
    file = HDF5File(mesh.mpi_comm(), name("phi.h5"), "w")
    file.write(phi, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("c.h5"), "w")
    file.write(c, "/function", 0)
    file.close()

    # Copy the mesh there as well
    shutil.copy(meshfile, name("sasflow%d.h5" % n))
    

    # Compute and output some key numbers
    phi_L2 = norm(phi, "L2")
    phi_Hdiv0 = norm(phi, "Hdiv0")
    omega = assemble(1*dx(domain=mesh))
    avg_phi = bar(phi, omega)
    c_L2 = norm(c, "L2")
    c1_L2 = norm(c1, "L2")

    values = (phi_L2, avg_phi, phi_Hdiv0, c_L2, c1_L2)
    header = ("|phi|_L2", "phi_avg (mm/h)", "|div phi|_L2" , "|c|_L2", "|c1|_L2")
    name = lambda s: os.path.join(results_dir, s)
    csvwrite(name("values.csv"), values, header)

    # Store images and phi to .pvd format as well
    name = lambda s: os.path.join(results_dir, "pvd", s)
    info("Storing images and data to %s" % name("..."))
    file = File(name("c_e.pvd"))
    file << c0
    file << c1

    file = File(name("c.pvd"))
    file << c

    file = File(name("phi.pvd"))
    if space != "CG":
        V = VectorFunctionSpace(mesh, "CG", 1)
        phi_CG = project(phi, V)
        file << phi_CG
    else:
        file << phi

    info("Reduced OCD computed successfully.\n")
    return (c, phi)

def run_ocd_reduced(num=None, n=None, beta=None, key=None):
    # Take patient number, mesh resolution parameter n (16, 32, 64),
    # beta > 0 and night/day as input

    # Handle arguments
    t0 = time.time()
    set_working_tape(Tape())

    if not num:
        num = sys.argv[1]
        n = int(sys.argv[2])
        beta = float(sys.argv[3])
        key = sys.argv[4]
        normalized = sys.argv[5]
        assert normalized in ["true", "false"]
        normalized = normalized == "true"
        assert isinstance(normalized, bool)
        imfolder = sys.argv[6]
        dispersion = float(sys.argv[7])
        

    info("Handling patient %s with mesh %d, beta=%.1e at %s with dispersion %s" % (num, n, beta, key, dispersion))

    # Use ./num as input directory and set output ./key as output
    # directory
    path = os.getcwd()
    patient = os.path.join(path, num)
    dirname = "results_red_ocd_pat%s_n%d_beta%.1e_%s_normalized%s_dispersion_%s" % (num, n, beta, key, normalized, dispersion)
    output_dir = os.path.join(patient, "constructed", dirname)
    info("Setting %s as output directory" % output_dir)

    # Create OCD(R) map
    (c, phi) = map_patient(patient, n, beta, key, output_dir, normalized, imfolder, dispersion)

    # Print timings and memory
    list_timings(TimingClear.keep, [TimingType.wall])
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t1 = time.time()

    info("\nMax memory usage (MB) = %g\n" % (mem/1024))
    info("\nTotal time (min) = %g\n" % ((t1 - t0)/60))

if __name__ == "__main__":

    # python3 compute_ocd_reduced.py 240 32 1.0e-04 night
    run_ocd_reduced()
