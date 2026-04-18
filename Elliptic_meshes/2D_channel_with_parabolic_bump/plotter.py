import numpy as np
import matplotlib.pyplot as plt

# load CSV
data = np.loadtxt("grid.csv", delimiter=",")

rows, cols2 = data.shape
cols = cols2 // 2  # because we stored (x,y) pairs

# reconstruct X and Y matrices
X = np.zeros((rows, cols))
Y = np.zeros((rows, cols))

for i in range(rows):
    for j in range(cols):
        X[i, j] = data[i, 2*j]
        Y[i, j] = data[i, 2*j + 1]

# plot grid lines
plt.figure(figsize=(10,5))

# horizontal lines
for i in range(rows):
    plt.plot(X[i, :], Y[i, :])

# vertical lines
for j in range(cols):
    plt.plot(X[:, j], Y[:, j])

plt.title("Grid in physical plane")
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")
plt.grid(True)

plt.show()
