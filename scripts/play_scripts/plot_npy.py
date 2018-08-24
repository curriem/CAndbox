import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('input_fl', help='input file')
args = parser.parse_args()


input_fl  = np.load(args.input_fl)
time, n, m = input_fl.shape

for t in range(time):
    plt.figure()
    plt.imshow(input_fl[t, :, :])
    plt.savefig('../plots/output%s.png' % str(t).zfill(3))
    plt.close()
