import sys
import resource
import os
import os.path
import shutil
import glob

import datetime

import numpy
import nibabel
from nibabel.affines import apply_affine

from dolfin import *

def compute_omt(I0, I1, tau, alpha, space="CG"):
    """Given two image intensity fields I0 and I1, a time step tau > 0 and
    a regularization parameter alpha > 0, return a vector field phi
    such that
    
      I_t + div (phi I) \approx 0  in \Omega
    where I(t = 0) = I0 and I(t = tau) = I1. 
    Mueller, M., Karasev, P., Kolesov, I. & Tannenbaum, A. Optical
    flow estimation for flame detection in videos. IEEE Transactions on
    image processing, 22, 2786-2797
    """

    # Check balance of mass assumption
    c0 = assemble(I0*dx())
    c1 = assemble(I1*dx())
    info("Is mass preserved (\int I0 == \in I1)? I0 = %.5g, I1 = %0.5g" % (c0, c1))

    Q = I0.function_space()
    mesh = Q.mesh()

    if space == "RT":
        V = FunctionSpace(mesh, "RT", 1)
    else:
        V = VectorFunctionSpace(mesh, "CG", 1)
        
    phi = TrialFunction(V)
    psi = TestFunction(V)

    # Define dI and barI
    dI = (I1 - I0)/tau
    barI = 0.5*(I0 + I1)
    I = barI
    a = (inner(div(I*phi), div(I*psi)) + alpha*inner(I*phi, psi))*dx()
    L = - inner(dI, div(I*psi))*dx()

    phi = Function(V)
    A = assemble(a)
    b = assemble(L)

    # Compensate for zero rows in A:
    A.ident_zeros()
    
    solve(A, phi.vector(), b, "mumps")

    info("V.dim() = %d" % V.dim())
    info("V.dim(), local = %d" % phi.vector().get_local().size)

    return phi

def map_image_onto_mesh(filename, mesh, maskfile, is_bc=False, tkr=True):
    """
    Parameters:
    filename:       File to read image data from
    mesh:           FEniCS Mesh for mapping image onto
    Optional parameters:
    tkr: ?? 
    Returns FEniCS Function c with image data represented on the given mesh
    """
    
    # Read image data from given file 
    info("Reading image data from %s" % filename)
    image = nibabel.load(filename)
    data = image.get_fdata()
    data = numpy.nan_to_num(data)

    mask = nibabel.load(maskfile)
    mask_data = mask.get_fdata()

    # Extract RAS-to-voxel-space mapping:
    if not tkr:
        ras2vox = image.header.get_ras2vox()
    else:
        ras2vox = numpy.linalg.inv(image.header.get_vox2ras_tkr())
    # Create continuous piecewise linear function space
    V = FunctionSpace(mesh, "CG", 1)
    
    # Map V's dof coordinates (vertices (x, y, z)) into voxel space
    # (i, j, k)
    xyz = V.tabulate_dof_coordinates()
    ijk = apply_affine(ras2vox, xyz).T
    i,j,k = numpy.rint(ijk).astype("int")
    
    if max(i) >= 256:
        i[i>=256] = 255
    if max(j) >= 256:
        j[j>=256] = 255
    if max(k) >= 256:
        k[k>=256] = 255

    if is_bc:
        j_boundary_height = 200
        z_boundary_height = xyz[numpy.where(j==j_boundary_height)[0][0]][2]
        bc_data = data
        Mmin = 3
        Mmax = 10
        for n in range(len(i)):
            i2 = i[n]
            j2 = j[n]
            k2 = k[n]
            if j2 > j_boundary_height:
                if mask_data[i2,j2,k2] == 0 and data[i2,j2,k2] == 0:
                    for m in range(Mmin, Mmax):

                        # Extract the data values from the neighborhood 
                        values = data[i2-m:i2+m+1, j2-m:j2+m+1, k2-m:k2+m+1]

                        # Reshape values from (2m+1, 2m+1, 2m+1) to list:
                        v = values.reshape(1, -1)

                        # Identify unique non-zero (positive) values and
                        # the number of each 
                        pairs, counts = numpy.unique(v[v > 0],
                                                    return_counts=True)  

                        # Return the most common non-zero tag:
                        success = counts.size > 0
                        if success:
                            new_value = numpy.average(pairs) # pairs[counts.argmax()]
                            break 
                    bc_data[i2,j2,k2] = new_value
            else:
                bc_data[i2,j2,k2] = 0
        c = Function(V)
        c.rename("c", "Concentration")
        c.vector()[:] = bc_data[i, j, k]
        return c, z_boundary_height

    # Create discrete function v on V and map image onto it
    # (vectorized)
    c = Function(V)
    c.rename("c", "Concentration")
    c.vector()[:] = data[i, j, k]

    return c


# def find_image_pair(images, night=True):

#     if night:
#         info("Finding two images and times overnight (~24h to 48h)")
#     else:
#         info("Finding two images and times during the first day (~0h to 2-12h)")

#     # Read dates/times of all images
#     times = []
#     print(images)
#     for image in images:
#         date, time = os.path.basename(image).split("_")
#         time, _ = os.path.basename(time).split(".")
#         d = datetime.date(int(date[0:4]), int(date[4:6]), int(date[6:8]))
#         t = datetime.time(int(time[0:2]), int(time[2:4]), int(time[4:6]))
#         times += [datetime.datetime.combine(d, t)]

#     if night:
#         # For day after to day day after 24h to 48h    
#         dt0 = datetime.timedelta(hours=8)
#         dt1 = datetime.timedelta(hours=30)
#         dt2 = datetime.timedelta(hours=60)
        
#         # Extract pair of images between dt0 and dt1, and dt1 and dt1
#         t0 = min(times)
#         for (i, t) in enumerate(times):
#             dt = t - t0
#             if dt0 < dt < dt1:
#                 t1 = t
#                 i1 = i
#             if dt1 < dt < dt2:
#                 t2 = t
#                 i2 = i
#     else:
#         # Morning to afternoon: 1st image and last image on Day 1 of exams.
#         t0 = min(times)
#         t2 = min(times)
#         dt0 = datetime.timedelta(hours=0.5)
#         dt1 = datetime.timedelta(hours=12)
#         for (i, t) in enumerate(times):
#             dt = t - t0 
#             if dt < dt0:
#                 t1 = t
#                 i1 = i
#             if (dt0 < dt < dt1) and t > t2: 
#                 t2 = t
#                 i2 = i
#     info("Images at %s and %s with dt = tau (h) %s\n" % (str(t1), str(t2), str(t2-t1)))
        
    # Return image filenames
    # return (images[i1], images[i2], t2-t1)

def find_image_pair(images, key):

    # Read dates/times of all images
    times = []
    for image in images:
        date, time = os.path.basename(image).split("_")
        time, _ = os.path.basename(time).split(".")
        d = datetime.date(int(date[0:4]), int(date[4:6]), int(date[6:8]))
        t = datetime.time(int(time[0:2]), int(time[2:4]), int(time[4:6]))
        times += [datetime.datetime.combine(d, t)]

    # Set all times to the default (the smallest time point, should be 0.0)
    t0 = min(times)
    t1 = min(times)
    t2 = min(times)

    # "0h" refers to pair of first image at 0h and image after 4.9h - 6.7h ('6h')
    #! Changed to refer to the image after 1h
    if key == "0h":
        dt0 = datetime.timedelta(hours=0)
        # dt1 = datetime.timedelta(hours=4.9)
        # dt2 = datetime.timedelta(hours=6.7)
        dt1 = datetime.timedelta(hours=0.75)
        dt2 = datetime.timedelta(hours=1)

        # Set t1 to be the first time and index
        t1 = min(times)
        i1 = 0

        # Find the second time and index
        for (i, t) in enumerate(times):
            dt = t - t0 
            if (dt1 < dt < dt2) and t > t2: 
                t2 = t
                i2 = i

    # "6h" refers to pair of image after 4.9h - 6.7h ('6h') and image 8-30 hours
    elif key == "6h":
        # For day after to day day after 24h to 48h    
        dt0 = datetime.timedelta(hours=4.9)
        dt1 = datetime.timedelta(hours=7.0)  #! 6.7 -> 7.0
        dt2 = datetime.timedelta(hours=26.5) #! 26 -> 26.5

        print("dt0: ", dt0)
        print("dt1: ", dt1)
        print("dt2: ", dt2)
        
        # Extract pair of images between dt0 and dt1, and dt1 and dt1
        for (i, t) in enumerate(times):
            dt = t - t0
            print("dt: ", dt)
            if dt0 < dt < dt1:
                t1 = t
                i1 = i
            elif dt1 < dt < dt2:
                t2 = t
                i2 = i
            else:
                pass

    elif key == "24h":
        # For day after to day day after 24h to 48h    
        dt0 = datetime.timedelta(hours=7)
        dt1 = datetime.timedelta(hours=30)
        dt2 = datetime.timedelta(hours=60)
        
        # Extract pair of images between dt0 and dt1, and dt1 and dt1
        for (i, t) in enumerate(times):
            dt = t - t0
            if dt0 < dt < dt1:
                t1 = t
                i1 = i
            elif dt1 < dt < dt2:
                t2 = t
                i2 = i
            else:
                pass

    else:
        info("Unrecognized image pair starting at %s, exiting" % key)
        exit()

    if (i2 == i1):
        info("Unable to identify pair of images, indices are: %d, %d" % (i1, i2))
        exit()

    info("Images at %s and %s with dt = tau (h) %s\n" % (str(t1), str(t2), str(t2-t1)))

    # Return image filenames
    return (images[i1], images[i2], t2-t1)

def map_patient(patient_dir, n, beta, key, results_dir):
    """
    patient_dir: Absolute patient directory e.g. "../../241"
    n:           Mesh resolution e.g. 16 (depending on available meshes)
    beta:        Percentage factor to use for regularization
    night:       True if comparison should be during night, false if during day 1.
    results_dir: Name for results directory (subdirectories will be created)
    """

    # Read FEniCS mesh 
    mesh = Mesh()
    meshfile = os.path.join(patient_dir, "mesh", "parenchyma%d_with_DTI.h5" % n)
    if not os.path.isfile(meshfile):
        meshfile = os.path.join(patient_dir, "mesh", "sasflow%d.h5" % n)
        assert os.path.isfile(meshfile)
        # raise Exception("Missing meshfile %s" % meshfile)

    info("Reading mesh from %s" % meshfile)
    hdf = HDF5File(mesh.mpi_comm(), meshfile, "r")
    hdf.read(mesh, "/mesh", False)
    #hdf.read(subdomains, "/subdomains", False) # Gray, white, brainstem
    # All FreeSurfer regions available in the stored mesh:
    #hdf.read(lookup_table, "/lookup_table", False)
    hdf.close()
    info("... with %d vertices and %d cells" % (mesh.num_vertices(), mesh.num_cells()))

    # Location of the mask file
    maskfile = os.path.join(patient_dir, "mesh", "sasparc.mgz")
    
    # Read set of images for this patient
    image_dir = os.path.join(patient_dir, "FIGURES_CONC_LUT")
    assert(os.path.isdir(image_dir)), "Missing image directory %s" % image_dir
    images = glob.glob(os.path.join(image_dir, "*.mgz"))

    # Find the images we want
    i0, i1, dt = find_image_pair(images, key)
    tau = dt.seconds/(60*60) # Compute rate in hours

    # Map images onto the mesh
    c0 = map_image_onto_mesh(i0, mesh, maskfile, is_bc=False)
    c1 = map_image_onto_mesh(i1, mesh, maskfile, is_bc=False)

    info("Normalizing")
    c0_norm = assemble(c0*dx)
    c1_norm = assemble(c1*dx)
    c0.vector()[:] = numpy.divide(c0.vector()[:], c0_norm)
    c1.vector()[:] = numpy.divide(c1.vector()[:], c1_norm)
    
    # Guess at some regularization parameter
    m = norm(c0, "L2")
    alpha = Constant(beta*m)

    space = "CG" # Alternatives: "RT"/"CG"
    phi = compute_omt(c0, c1, tau, alpha, space=space)
    info("||\phi||_L2 = %g" % norm(phi, "L2"))

    # -----------------------------------------------------------
    # The following is boiler-plate post-processing code
    # -----------------------------------------------------------
    
    # Store images and solution to .h5 format for postprocessing
    name = lambda s: os.path.join(results_dir, "hdf5", s)
    info("Storing images and data to %s" % name("..."))
    file = HDF5File(mesh.mpi_comm(), name("c0.h5"), "w")
    file.write(c0, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("c1.h5"), "w")
    file.write(c1, "/function", 0)
    file.close()

    file = HDF5File(mesh.mpi_comm(), name("phi.h5"), "w")
    file.write(phi, "/function", 0)
    file.close()

    # Copy the mesh there as well
    if "with_DTI" in meshfile:
        shutil.copy(meshfile, name("parencyma%d_with_DTI.h5" % n))
    else:
        shutil.copy(meshfile, name("parencyma%d.h5" % n))
    # FIXME: Replace above with this:
    # shutil.copy(meshfile, name("parenchyma%d_with_DTI.h5" % n))
    
    # Store images and phi to .pvd format as well
    name = lambda s: os.path.join(results_dir, "pvd", s)
    info("Storing images and data to %s" % name("..."))
    file = File(name("c.pvd"))
    file << c0
    file << c1

    file = File(name("phi.pvd"))
    if space != "CG":
        CG1 = VectorFunctionSpace(mesh, "CG", 1)
        phi = project(phi, CG1)
    file << phi

    info("OMT computed successfully.\n")
    return (c0, c1, phi)

def run_omt():
    # Take patient number, mesh resolution parameter n (16, 32, 64),
    # beta > 0 and night/day as input
    # FIXME: argparse
    num = sys.argv[1]
    n = int(sys.argv[2])
    beta = float(sys.argv[3])
    # is_night = sys.argv[4] == "night"
    # night = "night" if is_night else "day"
    key = sys.argv[4]
    info("Handling patient %s with mesh %d, beta=%.2f at %s" % (num, n, beta, key))

    # Use ./num as input directory and set output ./key as output
    # directory
    path = os.getcwd()
    patient = os.path.join(path, num)
    dirname = "results_omt_pat%s_n%d_beta%.1e_%s" % (num, n, beta, key)
    output_dir = os.path.join(patient, dirname)
    info("Setting %s as output directory" % output_dir)

    # Create OMT map
    (c0, c1, phi) = map_patient(patient, n, beta, key, output_dir)

    # Print timings and memory
    list_timings(TimingClear.keep, [TimingType.wall])
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    info("\nMax memory usage (MB) = %g\n" % (mem/1024))

    ## Some notes on memory considerations
    # mpirun -n 1 python3 compute_omt.py 240 16
    # Max memory usage (MB) =  3444.0
    
    # mpirun -n 2 python3 compute_omt.py 240 16
    # Max memory usage (MB) =  2155.27734375

    # mpirun -n 1 python3 compute_omt.py 240 32
    # Max memory usage (MB) =  11677.03125

    # mpirun -n 2 python3 compute_omt.py 240 32
    # Max memory usage (MB) =  7276.44140625

    # mpirun -n 4 python3 compute_omt.py 240 32
    # Max memory usage (MB) =  3733.19921875

if __name__ == "__main__":

    run_omt()