
head = r"""
\documentclass{article}
\usepackage{graphicx,subcaption}
\begin{document}

\textwidth=12cm
"""

tail = """ 
\end{document}
"""

pat_str = """
\\begin{figure}[!htb]
\\begin{center}
\\begin{subfigure}{0.50\\textwidth}
\includegraphics[width=0.95\\textwidth]{%s/phi.png}
\caption{Velocity field.}
\end{subfigure}
\\begin{subfigure}{0.50\\textwidth}
\includegraphics[width=0.95\\textwidth]{%s/velocities.png}
\caption{Max and average velocity as a function of iterations.}
\end{subfigure}
\\begin{subfigure}{0.50\\textwidth}
\includegraphics[width=0.95\\textwidth]{%s/j_error.png}
\caption{Evolution of the functional and regularization term over iterations.}
\end{subfigure}
\caption{Solving problem %s with patient %s, mesh size %s, regularization %s, dt %s with dispersion %s and method %s.}
\end{center}
\end{figure}
\\newpage
""" 

pat_str_omt = """
\\begin{figure}[!htb]
\\begin{center}
\includegraphics[width=0.95\\textwidth]{%s/phi.png}
\caption{Solving problem %s with patient %s, mesh size %s, regularization %s, dt %s with dispersion %s and method %s.}
\end{center}
\end{figure}
\\newpage
""" 

patient = "PatID-004_cut"
imgs = ["data", "z_vel", "stokes_vel", "z_vel_noise_10"]
nums = ["1"]
methods = ["ocdr", "ocdr_z_vel", "omt", "ocdr_stokes"]
dispersion = 1

file_str = head

for img in imgs:
    for method in methods:
        if method == "ocdr_stokes":
            if img == "z_vel":
                img = "not_scaled_z_vel"
            elif img == "stokes_vel":
                img = "not_scaled_stokes_vel"
            elif img == "z_vel_noise_10":
                img = "not_scaled_z_vel_noise_10"
            regs = ["1"]
        elif img == "data":
            regs = ["1e-2", "1e-6"]
        else:
            regs = ["1e-4", "1e-6"]
        for num in nums:
            for reg in regs:
                if img == "data":
                    dt = "0h"

                    pvd_folder_name = "%s/saga_pvd/%s/pat_%s_n_%s_beta_%s_key_%s_disp_1_method_%s" \
                            % (patient, img, patient, num, reg, dt, method)
                    
                    graph_folder_name = "%s/saga_files/%s/pat_%s_n_%s_beta_%s_key_%s_disp_1_method_%s" \
                            % (patient, img, patient, num, reg, dt, method)
                    
                    if method == "omt": 
                      file_str += pat_str_omt % (pvd_folder_name, \
                                            img, patient.split('_')[0], num, reg, dt, dispersion, method.replace('_', '\_'))

                    else:
                      file_str += pat_str % (pvd_folder_name, graph_folder_name, graph_folder_name, \
                                           img, patient.split('_')[0], num, reg, dt, dispersion, method.replace('_', '\_'))
                      
                else:
                  dt = "0_165"
                  pvd_folder_name = "%s/saga_pvd/%s/pat_%s_n_%s_beta_%s_dt_%s_disp_1_method_%s" \
                          % (patient, img, patient, num, reg, dt, method)
                  
                  graph_folder_name = "%s/saga_files/%s/pat_%s_n_%s_beta_%s_dt_%s_disp_1_method_%s" \
                          % (patient, img, patient, num, reg, dt, method)
                  
                  if method == "omt": 
                    file_str += pat_str_omt % (pvd_folder_name, \
                                          img.replace('_', '\_'), patient.split('_')[0], num, reg, dt.replace('_', '\_'), dispersion, method.replace('_', '\_'))

                  else:
                    file_str += pat_str % (pvd_folder_name, graph_folder_name, graph_folder_name, \
                                        img.replace('_', '\_'), patient.split('_')[0], num, reg, dt.replace('_', '\_'), dispersion, method.replace('_', '\_'))


with open("tables.txt", "r") as f:
    tables = f.read()

file_str += r"""
\begin{verbatim}
"""

file_str += tables

file_str += r"""
\end{verbatim}
"""
                    
file_str += tail 

f = open ("reportOCD.tex", "w") 
f.write(file_str)
f.close()




