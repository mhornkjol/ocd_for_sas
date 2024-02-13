import numpy as np
import csv
import sys
import os

patients = ["PatID-004_cut"]
imgs = ["data", "z_vel", "stokes_vel", "z_vel_noise_10"]
methods = ["ocdr", "ocdr_z_vel", "omt", "ocdr_stokes"]

print("The true max velocity for the manufactured solutions is: 20 mm/h")
print("The true average velocity for the manufactured stokes solutions is: 2.7716195960520817 mm/h")


for i, patient in enumerate(patients):
    print("#########################################")
    print("###    %s  Patient: %s       ###" % (i+1, patient))
    print("#########################################")
    for j, img in enumerate(imgs):
        print("#########################################")
        print("###    %s.%s  Problem: %s             ###" % (i+1, j+1, img))
        # if img == "z_vel_noise_10":
        #     methods = ["ocdr", "ocdr_z_vel", "omt"]
        print("#########################################")
        for k, method in enumerate(methods):
            print("#########################################")
            print("###    %s.%s.%s  Method: %s            ###" % (i+1, j+1, k+1, method))
            print("#########################################")
            # if method == "ocdr_z_vel":
            #     print("!!! NOTE: This method is using a different formulation. where we have not done partial integration in the forward problem. !!!")
            if method == "ocdr_stokes":
                print("!!! NOTE: This method is not done on the scaled manufacutred soluton. The images has a higher intesity than the others.")
            # DATA
            if img == "data":
                resolutions = ["1", "2", "4"]
                if method == "ocdr_stokes":
                    alphas = ["1"]
                else:
                    alphas = ["1e-2", "1e-4", "1e-6"]
                dt = "0h"
                disp = "1"

                mphis = []
                aphis = []
                j_reduc = []
                j_ratio = []
                for res in resolutions:
                    for alpha in alphas:
                        folder_name = "%s/saga_files/%s/pat_%s_n_%s_beta_%s_key_%s_disp_%s_method_%s" \
                                    % (patient, img, patient, res, alpha, dt, disp, method)
                        with open(os.path.join(folder_name, "optimization_values.csv"), "r") as f:
                            reader = csv.reader(f)
                            flag = True
                            for row in reader:
                                if row[0] != "j":
                                    if flag:
                                        j_init = float(row[0])
                                        flag = False
                                    mphi= float(row[3])
                                    aphi = float(row[4])
                                    j_final = float(row[0])
                                    jr_final = float(row[1])
                                    j0_final = float(row[2])

                            mphis.append(mphi)
                            aphis.append(aphi)
                            if method != "omt":
                                j_reduc.append((j_init-j_final)*100/j_init)
                            if method != "ocdr_stokes" and method != "omt":
                                j_ratio.append(j0_final/jr_final)

                # OCDR stokes vel
                if method == "ocdr_stokes":
                    print("Table Max Vel (mm/h), %s, %s" % (patient, img))
                    print("res \ alpha %s" % (alphas[0]))
                    print("--------------------------")
                    print("1           %0.2f" % (mphis[0]))
                    print("2           %0.2f" % (mphis[1]))
                    print("4           %0.2f" % (mphis[2]))

                    print("\n")

                    print("Table Avg Vel (mm/h), %s, %s" % (patient, img))
                    print("res \ alpha %s" % (alphas[0]))
                    print("--------------------------")
                    print("1           %0.2f" % (aphis[0]))
                    print("2           %0.2f" % (aphis[1]))
                    print("4           %0.2f" % (aphis[2]))

                    print("\n")

                    print("Table J reduction (%%), %s, %s" % (patient, img))
                    print("res \ alpha %s" % (alphas[0]))
                    print("--------------------------")
                    print("1           %0.2f" % (j_reduc[0]))
                    print("2           %0.2f" % (j_reduc[1]))
                    print("4           %0.2f" % (j_reduc[2]))

                    print("\n")
                
                # omt method
                elif method == "omt":
                    print("Table Max Vel (mm/h), %s, %s" % (patient, img))
                    print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f" % (mphis[0], mphis[1], mphis[2]))
                    print("2           %0.2f %0.2f %0.2f" % (mphis[3], mphis[4], mphis[5]))
                    print("4           %0.2f %0.2f %0.2f" % (mphis[6], mphis[7], mphis[8]))

                    print("\n")

                    print("Table Avg Vel (mm/h), %s, %s" % (patient, img))
                    print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f" % (aphis[0], aphis[1], aphis[2]))
                    print("2           %0.2f %0.2f %0.2f" % (aphis[3], aphis[4], aphis[5]))
                    print("4           %0.2f %0.2f %0.2f" % (aphis[6], aphis[7], aphis[8]))

                    print("\n")

                # Other methods
                else:
                    print("Table Max Vel (mm/h), %s, %s" % (patient, img))
                    print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f" % (mphis[0], mphis[1], mphis[2]))
                    print("2           %0.2f %0.2f %0.2f" % (mphis[3], mphis[4], mphis[5]))
                    print("4           %0.2f %0.2f %0.2f" % (mphis[6], mphis[7], mphis[8]))

                    print("\n")

                    print("Table Avg Vel (mm/h), %s, %s" % (patient, img))
                    print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f" % (aphis[0], aphis[1], aphis[2]))
                    print("2           %0.2f %0.2f %0.2f" % (aphis[3], aphis[4], aphis[5]))
                    print("4           %0.2f %0.2f %0.2f" % (aphis[6], aphis[7], aphis[8]))

                    print("\n")

                    print("Table J reduction (%%), %s, %s" % (patient, img))
                    print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f" % (j_reduc[0], j_reduc[1], j_reduc[2]))
                    print("2           %0.2f %0.2f %0.2f" % (j_reduc[3], j_reduc[4], j_reduc[5]))
                    print("4           %0.2f %0.2f %0.2f" % (j_reduc[6], j_reduc[7], j_reduc[8]))

                    print("\n")

                    print("Table J ratio (j0/j_r), %s, %s" % (patient, img))
                    print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f" % (j_ratio[0], j_ratio[1], j_ratio[2]))
                    print("2           %0.2f %0.2f %0.2f" % (j_ratio[3], j_ratio[4], j_ratio[5]))
                    print("4           %0.2f %0.2f %0.2f" % (j_ratio[6], j_ratio[7], j_ratio[8]))

                    print("\n")

            # Manufactured solutions
            else:
                resolutions = ["1", "2", "4"]
                if method == "ocdr_stokes":
                    alphas = ["1"]
                    if img == "z_vel":
                        img = "not_scaled_z_vel"
                    elif img == "stokes_vel":
                        img = "not_scaled_stokes_vel"
                    elif img == "z_vel_noise_10":
                        img = "not_scaled_z_vel_noise_10"
                else:
                    alphas = ["1e-4", "1e-5", "1e-6"]
                dts = ["0_42", "0_83", "0_125", "0_165"]
                disp = "1"

                # Set dt other methods
                if method != "ocdr_stokes":
                    dt = dts[3]
                    print("#------- For set dt %s ----------#" % dt)
                    mphis = []
                    aphis = []
                    j_reduc = []
                    j_ratio = []
                    for res in resolutions:
                        for alpha in alphas:
                            folder_name = "%s/saga_files/%s/pat_%s_n_%s_beta_%s_dt_%s_disp_%s_method_%s" \
                                        % (patient, img, patient, res, alpha, dt, disp, method)
                            with open(os.path.join(folder_name, "optimization_values.csv"), "r") as f:
                                reader = csv.reader(f)
                                flag = True
                                for row in reader:
                                    if row[0] != "j":
                                        if flag:
                                            j_init = float(row[0])
                                            flag = False
                                        mphi= float(row[3])
                                        aphi = float(row[4])
                                        j_final = float(row[0])
                                        jr_final = float(row[1])
                                        j0_final = float(row[2])

                                mphis.append(mphi)
                                aphis.append(aphi)
                                if method != "omt":
                                    j_reduc.append((j_init-j_final)*100/j_init)
                                    j_ratio.append(j0_final/jr_final)
                    
                    if method == "omt":
                        print("Table Max Vel (mm/h), %s, %s, %s" % (patient, img, dt))
                        print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                        print("--------------------------")
                        print("1           %0.2f %0.2f %0.2f" % (mphis[0], mphis[1], mphis[2]))
                        print("2           %0.2f %0.2f %0.2f" % (mphis[3], mphis[4], mphis[5]))
                        print("4           %0.2f %0.2f %0.2f" % (mphis[6], mphis[7], mphis[8]))

                        print("\n")

                        print("Table Avg Vel (mm/h), %s, %s, %s" % (patient, img, dt))
                        print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                        print("--------------------------")
                        print("1           %0.2f %0.2f %0.2f" % (aphis[0], aphis[1], aphis[2]))
                        print("2           %0.2f %0.2f %0.2f" % (aphis[3], aphis[4], aphis[5]))
                        print("4           %0.2f %0.2f %0.2f" % (aphis[6], aphis[7], aphis[8]))

                        print("\n")

                    else:
                        print("Table Max Vel (mm/h), %s, %s, %s" % (patient, img, dt))
                        print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                        print("--------------------------")
                        print("1           %0.2f %0.2f %0.2f" % (mphis[0], mphis[1], mphis[2]))
                        print("2           %0.2f %0.2f %0.2f" % (mphis[3], mphis[4], mphis[5]))
                        print("4           %0.2f %0.2f %0.2f" % (mphis[6], mphis[7], mphis[8]))

                        print("\n")

                        print("Table Avg Vel (mm/h), %s, %s, %s" % (patient, img, dt))
                        print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                        print("--------------------------")
                        print("1           %0.2f %0.2f %0.2f" % (aphis[0], aphis[1], aphis[2]))
                        print("2           %0.2f %0.2f %0.2f" % (aphis[3], aphis[4], aphis[5]))
                        print("4           %0.2f %0.2f %0.2f" % (aphis[6], aphis[7], aphis[8]))

                        print("\n")

                        print("Table J reduction (%%), %s, %s, %s" % (patient, img, dt))
                        print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                        print("--------------------------")
                        print("1           %0.2f %0.2f %0.2f" % (j_reduc[0], j_reduc[1], j_reduc[2]))
                        print("2           %0.2f %0.2f %0.2f" % (j_reduc[3], j_reduc[4], j_reduc[5]))
                        print("4           %0.2f %0.2f %0.2f" % (j_reduc[6], j_reduc[7], j_reduc[8]))

                        print("\n")

                        print("Table J ration (j0/j_r), %s, %s, %s" % (patient, img, dt))
                        print("res \ alpha %s %s %s" % (alphas[0], alphas[1], alphas[2]))
                        print("--------------------------")
                        print("1           %0.2f %0.2f %0.2f" % (j_ratio[0], j_ratio[1], j_ratio[2]))
                        print("2           %0.2f %0.2f %0.2f" % (j_ratio[3], j_ratio[4], j_ratio[5]))
                        print("4           %0.2f %0.2f %0.2f" % (j_ratio[6], j_ratio[7], j_ratio[8]))

                        print("\n\n")


                alpha = alphas[0]
                print("#------- For set alpha %s ----------#" % alpha)
                mphis = []
                aphis = []
                j_reduc = []
                j_ratio = []
                for res in resolutions:
                    for dt in dts:
                        folder_name = "%s/saga_files/%s/pat_%s_n_%s_beta_%s_dt_%s_disp_%s_method_%s" \
                                    % (patient, img, patient, res, alpha, dt, disp, method)
                        with open(os.path.join(folder_name, "optimization_values.csv"), "r") as f:
                            reader = csv.reader(f)
                            flag = True
                            for row in reader:
                                if row[0] != "j":
                                    if flag:
                                        j_init = float(row[0])
                                        flag = False
                                    mphi= float(row[3])
                                    aphi = float(row[4])
                                    j_final = float(row[0])
                                    jr_final = float(row[1])
                                    j0_final = float(row[2])

                            mphis.append(mphi)
                            aphis.append(aphi)
                            if method != "omt":
                                j_reduc.append((j_init-j_final)*100/j_init)
                            if method != "ocdr_stokes" and method != "omt":
                                j_ratio.append(j0_final/jr_final)

                # Stokes vel
                if method == "ocdr_stokes":
                    print("Table Max Vel (mm/h), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (mphis[0], mphis[1], mphis[2], mphis[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (mphis[4], mphis[5], mphis[6], mphis[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (mphis[8], mphis[9], mphis[10], mphis[11]))

                    print("\n")

                    print("Table Avg Vel (mm/h), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (aphis[0], aphis[1], aphis[2], aphis[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (aphis[4], aphis[5], aphis[6], aphis[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (aphis[8], aphis[9], aphis[10], aphis[11]))

                    print("\n")

                    print("Table J reduction (%%), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (j_reduc[0], j_reduc[1], j_reduc[2], j_reduc[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (j_reduc[4], j_reduc[5], j_reduc[6], j_reduc[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (j_reduc[8], j_reduc[9], j_reduc[10], j_reduc[11]))

                    print("\n")

                # omt method
                elif method == "omt":
                    print("Table Max Vel (mm/h), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (mphis[0], mphis[1], mphis[2], mphis[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (mphis[4], mphis[5], mphis[6], mphis[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (mphis[8], mphis[9], mphis[10], mphis[11]))

                    print("\n")

                    print("Table Avg Vel (mm/h), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (aphis[0], aphis[1], aphis[2], mphis[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (aphis[4], aphis[5], aphis[6], mphis[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (aphis[8], aphis[9], aphis[10], mphis[11]))

                    print("\n")

                # Other methods
                else:
                    print("Table Max Vel (mm/h), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (mphis[0], mphis[1], mphis[2], mphis[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (mphis[4], mphis[5], mphis[6], mphis[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (mphis[8], mphis[9], mphis[10], mphis[11]))

                    print("\n")

                    print("Table Avg Vel (mm/h), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (aphis[0], aphis[1], aphis[2], mphis[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (aphis[4], aphis[5], aphis[6], mphis[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (aphis[8], aphis[9], aphis[10], mphis[11]))

                    print("\n")

                    print("Table J reduction (%%), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (j_reduc[0], j_reduc[1], j_reduc[2], j_reduc[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (j_reduc[4], j_reduc[5], j_reduc[6], j_reduc[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (j_reduc[8], j_reduc[9], j_reduc[10], j_reduc[11]))

                    print("\n")

                    print("Table J ration (j0/j_r), %s, %s, %s" % (patient, img, alpha))
                    print("res \ dt %s %s %s %s" % (dts[0], dts[1], dts[2], dts[3]))
                    print("--------------------------")
                    print("1           %0.2f %0.2f %0.2f %0.2f" % (j_ratio[0], j_ratio[1], j_ratio[2], j_ratio[3]))
                    print("2           %0.2f %0.2f %0.2f %0.2f" % (j_ratio[4], j_ratio[5], j_ratio[6], j_ratio[7]))
                    print("4           %0.2f %0.2f %0.2f %0.2f" % (j_ratio[8], j_ratio[9], j_ratio[10], j_ratio[11]))

                    print("\n")