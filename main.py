import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import griddata
from scipy.ndimage import zoom


from utils import *

Nx, Ny = 150, 250

Lx, Ly = 5, 8

# Fréquence de l'onde
F = 2.5e9

# Célérité de l'onde
c = 299792458   # 299792458m

# Pulsation
pulsation = 2 * np.pi * F

# Nombre d'ondes
k = pulsation / c


def fs(x, y):
    """
    Fonction chapeau f_s(x) qui approxime la distribution de Dirac autour de (x_s, y_s).

    Paramètres :
    - x, y : coordonnées où évaluer la fonction (peuvent être des tableaux NumPy).
    - x_s, y_s : coordonnées du centre du disque Ω_s.
    - eps_s : rayon du disque Ω_s.

    Retourne :
    - La valeur de f_s(x, y) sous forme d'un tableau NumPy.
    """
    
    eps_s = 1e-1
    x_s, y_s = x_s, y_s = 2 * eps_s * Lx , eps_s * Ly
    # Distance euclidienne entre (x, y) et (x_s, y_s)
    dist = np.sqrt((x - x_s)**2 + (y - y_s)**2)
    
    # Appliquer la définition de f_s
    inside = dist < eps_s  # Condition pour être dans Ω_s
    f_values = np.zeros_like(x)  # Initialisation avec 0
    
    f_values[inside] = (3 / (np.pi * eps_s**2)) * (1 - dist[inside] / eps_s)
    
    return f_values


nodes, connectivity = triangular_discretisation(nx=Nx, ny=Ny, Lx=Lx, Ly=Ly)

grid = discretiser_maison(Nx, Ny)

mef = MEF(grid=grid, nodes=nodes, connectivity=connectivity)


A, F = mef.assemblage(fs, k, Lx, Ly)

# Supposons que A est votre matrice globale (dense actuellement) et F votre vecteur.
# Convertir A en format sparse:
A_sparse = sp.csr_matrix(A)

# Résoudre le système avec spsolve
U = spla.spsolve(A_sparse, F)

# plot
fig = plt.figure(figsize=(10, 15))
ax = fig.add_subplot(111, projection='3d')

# Reshape u to match the grid shape
u_reshaped = U.reshape((Ny, Nx))

# sauvegarde de la solution
np.save('solution.npy', u_reshaped)

# # Create a meshgrid for plotting
# X, Y = np.meshgrid(np.linspace(0, Lx, Nx), np.linspace(0, Ly, Ny))

# u_interp = griddata(nodes, U, (X, Y), method='cubic' )

# print(u_reshaped.shape, X.shape, Y.shape)

# # Plot the surface
# ax.plot_surface(X, Y, np.abs(u_reshaped), cmap='viridis')

# # plt.contourf(X, Y, np.abs(u_interp), level=50, cmap='viridis')

# # ax.contour3D(X, Y, u_reshaped, 50, cmap='viridis')


# # Set labels
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('u')

# plt.show()
