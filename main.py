import numpy as np
from calcul import random_step, number_neighbors
from visualisation import visu3D
import matplotlib as mpl

# Grid and length params
N = 10**5
n = 10000

# Random initial start
x, y, z = np.random.randint(0, N-1, 3)
pos = [[x,y,z]]

# Looping on the walk
for step in range(n):
   if number_neighbors(N, [x,y,z], pos) == 6:
      # Stoping the walk when it reaches a closed-loop of neighbors
      break
   x, y, z = random_step(N, pos[step])
   while [x, y, z] in pos:
      # Generating new step if the step is already present in the history of steps
      x, y, z = random_step(N, pos[step])
   pos.append([x,y,z])

pos = np.array(pos).T

# Global plotting parameters
mpl.rcParams['font.size'] = 15
mpl.rcParams['figure.figsize'] = (6,5)
# mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']
mpl.rcParams['axes.linewidth'] = 3

visu3D(pos[0], pos[1], pos[2])