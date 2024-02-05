from fenics import *
import numpy as np
import csv
import os

patient = "PatID-004_cut"
imgs = ["not_scaled_z_vel", "not_scaled_stokes_vel"]
nums = ["4"]


def compute_magnitude_field(phi):

    # Split phi into its components (really split it)
    (phi0, phi1, phi2) = phi.split(deepcopy=True)

    # Get the map from vertex index to dof
    v2d = [vertex_to_dof_map(phii.function_space()) for phii in (phi0, phi1, phi2)]

    # Create array of component values at vertex indices so p0[0] is
    # value at vertex 0 e.g. 
    p0 = phi0.vector().get_local(v2d[0])
    p1 = phi1.vector().get_local(v2d[1])
    p2 = phi2.vector().get_local(v2d[2])

    # Take element-wise magnitude
    m = np.sqrt((np.square(p0) + np.square(p1) + np.square(p2)))

    # Map array of values at vertices back to values at dofs and
    # insert into magnitude field.
    M = phi0.function_space()
    magnitude = Function(M)
    d2v = dof_to_vertex_map(M)
    magnitude.vector()[:] = m[d2v]
    return magnitude

def bar(phi, volume=None):
    phi_mag = compute_magnitude_field(phi)
    if not volume:
        mesh = phi.function_space().mesh()
        volume = assemble(1.0*dx(domain=mesh))
    avg_phi_mag = assemble(phi_mag*dx())/volume
    return avg_phi_mag


for num in nums:
    mesh = Mesh()
    hdf = HDF5File(mesh.mpi_comm(), "%s/mesh/sasflow%s.h5" % (patient, num), "r")
    hdf.read(mesh, "/mesh", False)

    O = VectorFunctionSpace(mesh, 'CG', 1)

    phi = Function(O)


    xdmf = XDMFFile(MPI.comm_world, "%s/saga_stokes_vel_fields/stokes_vel_field_%s/velocity.xdmf" % (patient, num))
    xdmf.read_checkpoint(phi, "velocity", 0)

    # # Scale stokes to 20 mm/h max velocity magnitude
    # u, v, w = phi.split(deepcopy=True)
    # magnitude = np.sqrt(u.vector()[:]**2 + v.vector()[:]**2 + w.vector()[:]**2)
    # max_vel = 20
    # factor = max_vel/max(magnitude)
    # phi.vector()[:] = phi.vector()[:]*factor


    for img in imgs:
        if img == "data":
            dt = "0h"
            print(img, num, dt)
            folder_name = "%s/saga_files/%s/pat_%s_n_%s_beta_1_key_%s_disp_1_method_ocdr_stokes" \
                              % (patient, img, patient, num, dt)
            
            with open(os.path.join(folder_name, "optimization_values.csv"), "r") as infile, open(os.path.join(folder_name, "optimization_values_2.csv"), "w") as outfile:
                reader = csv.reader(infile)
                writer = csv.writer(outfile)
                # next(reader, None)
                i = 1
                for row in reader:
                    print("row", i)
                    if row[0] == "j":
                        writer.writerow(row)
                    else:
                        z_val = Constant(row[3])

                        new_phi = project(z_val*phi, O, solver_type='mumps')
                        # Should this be magnitude?
                        mphi = new_phi.vector().max()
                        aphi = bar(new_phi)
                        
                        new_row = row
                        new_row[3] = mphi
                        new_row[4] = aphi
                        writer.writerow(new_row)
                    i += 1
        else:
            dts = ["0_42", "0_83", "0_125", "0_165"]
            for dt in dts:
                print(img, num, dt)
                folder_name = "%s/saga_files/%s/pat_%s_n_%s_beta_1_dt_%s_disp_1_method_ocdr_stokes" \
                              % (patient, img, patient, num, dt)

                with open(os.path.join(folder_name, "optimization_values.csv"), "r") as infile, open(os.path.join(folder_name, "optimization_values_2.csv"), "w") as outfile:
                    reader = csv.reader(infile)
                    writer = csv.writer(outfile)
                    # next(reader, None)
                    i = 1
                    for row in reader:
                        print("row", i)
                        if row[0] == "j":
                            writer.writerow(row)
                        else:
                            z_val = Constant(row[3])

                            new_phi = project(z_val*phi, O, solver_type='mumps')
                            mphi = new_phi.vector().max()
                            aphi = bar(new_phi)
                            
                            new_row = row
                            new_row[3] = mphi
                            new_row[4] = aphi
                            writer.writerow(new_row)
                        i += 1
