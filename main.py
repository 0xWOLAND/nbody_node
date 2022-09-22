import numpy as np
import argparse
import matplotlib.pyplot as plt
import tkinter


from init_config import init_config

parser = argparse.ArgumentParser(description='Simulation Inputs')
parser.add_argument('count', metavar='N', type=int, help="Number of particles")

args = parser.parse_args()
__COUNT__ = int(args.count)

print("COUNT: ", __COUNT__)
ic = init_config(__COUNT__)
ic.init_rand()

plt.scatter(ic.xy()[0], ic.xy()[1])
plt.savefig("temp.png")