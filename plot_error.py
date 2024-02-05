import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import os

folder_name = sys.argv[1]

j = []
jr = []
j0 = []
mphi = []
aphi = []
with open(os.path.join(folder_name, "optimization_values.csv"), "r") as f:
    reader = csv.reader(f)
    next(reader, None)
    for row in reader:
        if row[0] != "j":
            j.append(float(row[0]))
            jr.append(float(row[1]))
            j0.append(float(row[2]))
            mphi.append(float(row[3]))
            aphi.append(float(row[4]))

# new_j = []
# new_jr = []
# new_j0 = []
# for i, value in enumerate(j):
#     if value < j[0]+1000:
#         new_j.append(value)
#     else:
#         new_j.append(j[i-1])

# for i, value in enumerate(jr):
#     if value < j[0]+1000:
#         new_jr.append(value)
#     else:
#         new_jr.append(jr[i-1])

# for i, value in enumerate(j0):
#     if value < j[0]+1000:
#         new_j0.append(value)
#     else:
#         new_j0.append(j0[i-1])


# plt.plot(new_j, color='b', label='J')
# plt.plot(new_jr, color='r', label='Reg')
# plt.plot(new_j0, color='g', label='Func')
plt.figure(1)
plt.plot(j, color='b', label='J')
plt.plot(jr, color='r', label='Reg')
plt.plot(j0, color='g', label='Func')
plt.xlabel("# iters")
plt.ylabel("Errors")
plt.legend()
plt.ylim(0, j[0]+j[0]*0.1)
# plt.show()
plt.savefig(os.path.join(folder_name, "j_error.png"))


plt.figure(2)
plt.plot(mphi, color='b', label='Max phi')
plt.plot(aphi, color='r', label='Avg Phi')
plt.xlabel("# iters")
plt.ylabel("Velocities (mm/h)")
plt.legend()
plt.ylim(0, max(mphi))
# plt.show()
plt.savefig(os.path.join(folder_name, "velocities.png"))
