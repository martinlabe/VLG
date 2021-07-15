import subprocess as sp
import numpy as np
import pandas as pd
import matplotlib
import os

matplotlib.use('TKAgg', force=True)
import matplotlib.pyplot as plt

n = 500
file = "../graphs/facebook_combined.txt"
args = ['0', file, '15']

output = "./out_" + file[10:]


def call(n, args):
    """ call n time the c program """
    os.remove(output)
    for i in range(n):
        if i % 10 == 0:
            print(i)
        sp.call(["../cmake-build-debug/VLG", *args])


call(n, args)

df = pd.read_csv(output, names=['diameter', 'bfs', 'time'])
print(df.head())
# the histogram of the data
x = df['bfs']
ni, bins, patches = plt.hist(x, 20, density=True, facecolor='b', alpha=0.75)

plt.xlabel('Number of BFS')
plt.ylabel('Density')
plt.legend(title=r'$mean={}, var={}$'.format(np.mean(x), np.var(x)))
plt.grid(True)
name = "0-" + file[10:-4] + "-15.png"
plt.savefig(name)
plt.show()
