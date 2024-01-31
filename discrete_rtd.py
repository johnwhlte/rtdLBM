import numpy as np
import matplotlib.pyplot as plt

num_particles = 1000
deltaT = 0.01

with open("run_log","r") as f:

    lines = f.readlines()

stats_line = []

for line in lines:
    line_split = line.split()

    if line_split[0] == '[captureStatistics]' and len(line_split) == 4:
        stats_line.append(int(line_split[3].split('=')[1]))
        print(line_split)

e_t = np.zeros(shape=np.asarray(stats_line).shape)

for i, line in enumerate(stats_line):
    if i == 0:
        summer = 0
    else:
        summer+= line
        e_t[i] = ((line - stats_line[i-1])/num_particles)

t = np.arange(0, deltaT*len(e_t), deltaT)

print(np.sum(e_t))
print(np.mean(e_t*t) / deltaT)
print(1000000*0.020*(np.pi*(0.005**2)) / 1.985)

plt.plot(t, e_t)
plt.savefig("test.png")   
