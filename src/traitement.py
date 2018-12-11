import matplotlib.pyplot as plt
import numpy as np

plt.close()

f=open("../toPlot/Psis201.dat", "r")
n = f.read().split("\n")
A=[i.split(" ") for i in n]
B=[[] for _ in n]
for i in range(1,len(A)):
    B[i] = filter(lambda a : len(a) > 0, A[i])
B = filter(lambda a : len(a)>1, B)
##plt.axis([-1, 1, -50, 50])
C = [[float(B[i][j]) for i in range(len(B))] for j in range(len(B[0]))]
X = np.linspace(-1, 1, len(C[0]))
for i in C:
    print("Y : " + str(i))
    plt.plot(X, i)
plt.savefig('../plotted/bonjour.svg')