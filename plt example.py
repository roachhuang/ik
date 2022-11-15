# Import Libraries

import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display, clear_output

# Create figure and subplot

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

# Define and update plot

for i in range(20):
    x = np.linspace(0, i, 100);
    y = np.cos(x)
    ax.set_xlim(0, i)
    ax.cla()
    ax.plot(x, y)
    display(fig)
    clear_output(wait = True)
    plt.pause(0.1)