from box import *
import numpy as np
import matplotlib.pyplot as plt


b = box(dim=np.array([10,10,10]))

b.insertNode(particle(), b.tree)
print(b.tree)
b.insertNode(particle(np.array([2,0,0])), b.tree)
print(b.tree)

fig, ax = plt.subplots(1,1)
b.createPlot()

fig.savefig('temp2.png')