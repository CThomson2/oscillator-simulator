import numpy as np
import matplotlib.pyplot as plt
from main import fermi_sim

f = open('userdata.txt', 'r')
data = f.readlines()
for i in range(len(data)):
    data[i] = data[i].strip()

data.append(False)

fig, ax = plt.subplots()

tau = 10

for i in range(10):
    fermi_sim(
        {'L': data[0], 'tf': data[1], 'N': data[2], 'k': data[3], 'rho': data[4],
            'alpha': data[5], 'amp_0': data[6], 'shape_0': data[7], 'anim': data[8]}, tau
    )