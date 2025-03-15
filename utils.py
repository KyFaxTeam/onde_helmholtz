import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def triangular_discretisation(nx, ny, Lx, Ly):
      """
      Returns the coordinates of the nodes and the connectivity matrix for a triangular discretisation of a rectangle.
      """
      # Create a linearly spaced array of nx points between 0 and Lx
      x = np.linspace(0, Lx, nx)
      
      # Create a linearly spaced array of ny points between 0 and Ly
      y = np.linspace(0, Ly, ny)
      
      # Create a meshgrid from x and y arrays
      X, Y = np.meshgrid(x, y, sparse=False, copy=True)
      
      # Flatten the meshgrid arrays and combine them into a single array of node coordinates
      nodes = np.array([X.ravel(), Y.ravel()]).T
      
      # Calculate the number of elements in the mesh
      n_elements = 2 * (nx - 1) * (ny - 1)
      
      # Initialize the connectivity matrix with zeros
      connectivity = np.zeros((n_elements, 3), dtype=int)
      
      # 
      element = 0
      
      # Loop over the y direction
      for i in range(ny - 1):
            # Loop over the x direction
            for j in range(nx - 1):
                  # Define the connectivity of the first triangle
                  connectivity[element] = [i * nx + j, i * nx + j + 1, (i + 1) * nx + j]
                  # Define the connectivity of the second triangle
                  connectivity[element + 1] = [(i + 1) * nx + j, i * nx + j + 1, (i + 1) * nx + j + 1]
                  # Increment the element counter
                  element += 2
      return nodes, connectivity

def plot_triangulation(nodes, connectivity):
      """
      Plots the triangulation of a rectangle.
      """
      # Create a figure and axis
      fig, ax = plt.subplots()
      
      # Plot the nodes
      ax.plot(nodes[:, 0], nodes[:, 1], '.', label='Nodes')
      
      # Loop over the elements
      for element in connectivity:
            # Get the coordinates of the nodes of the element
            x = nodes[element, 0]
            y = nodes[element, 1]
            # Plot the element
            ax.plot(np.append(x, x[0]), np.append(y, y[0]), 'k-')
      
      # Set the aspect of the plot to be equal
      ax.set_aspect('equal')
      
      # Set the labels of the axes
      ax.set_xlabel('x')
      ax.set_ylabel('y')
      
      # Display the legend
      ax.legend()
      
      # Display the plot
      plt.show()

def discretiser_maison(Nx=1000, Ny=500):
    """Crée une grille discrète du domaine représentant la maison.
    
    Paramètre :
    - N : nombre de points de discrétisation sur chaque axe
    
    Retourne :
    - grid : matrice 2D représentant le domaine (1 = mur, 0 = espace libre)
    """
    grid = np.zeros((Nx, Ny), dtype=int)  # 0 = espace libre

    # Échelle en fonction de N
    def scale_x(val):  
        return int(val * Nx / 100)
    
    def scale_y(val):  
        return int(val * Ny / 100)
    
    # Epaisseur des murs
    d = 5
    # Aller
    ay = scale_y(10)
    ax = scale_x(10)

    # Murs extérieurs
    grid[0:d, :] = 1
    grid[-d:, :] = 1
    grid[:, 0:d] = 1
    grid[:, -d:] = 1

    # Murs intérieurs
    # grid[0:scale_x(35), scale_y(30):scale_y(30)+d] = 1  
    grid[scale_x(50):-1, scale_y(30):scale_y(30)+d] = 1
    
    grid[scale_x(50):-1, scale_y(40):scale_y(40)+d] = 1
    grid[scale_x(50):scale_x(50)+d, scale_y(40):scale_y(50)] = 1
    
    grid[scale_x(50):-1, scale_y(70):scale_y(70)+d] = 1
    
    grid[scale_x(50):scale_x(50)+d, scale_y(60):scale_y(70)+d] = 1
    grid[scale_x(50):scale_x(50)+d, scale_y(70)+ay:-1] = 1
     
    grid[scale_x(35)-d:scale_x(35), scale_y(30):scale_y(60)+d] = 1
    # grid[0:scale_x(30), scale_y(70)+ay:-1] = 1
    
    
    # Salle de bain
    grid[scale_x(35)-d:scale_x(35), scale_y(80):-1] = 1
    grid[0:scale_x(35), scale_y(70):scale_y(70)+d] = 1
    
    # Cuisine
    grid[0:scale_x(20), scale_y(50):scale_y(50)+d] = 1

    return grid


from typing import Tuple

def isInWall(grid: np.ndarray, nodes, connectivity, element: int, Nx: float, Ny: float, Lx:float, Ly:float) -> Tuple[bool, Tuple[int, ...]]:
    # Calcul des coordonnées physiques des nœuds de l'élément
    # xp = nodes[connectivity[element]] * np.array([Nx - 1, Ny - 1])
    xp = nodes[connectivity[element]] * np.array([(Nx-1)/Lx, (Ny-1)/Ly])

    inBord = []

    # Pour chaque nœud de l'élément, vérifier si le nœud est sur le bord du domaine
    for i in range(3):
        xx, yy = int(np.round(xp[i, 0])), int(np.round(xp[i, 1]))
        if xx == 0 or xx == Nx - 1 or yy == 0 or yy == Ny - 1:
            inBord.append(i)
    
    count = 0
        
    for i in range(3):
        xx, yy = int(np.round(xp[i, 0])), int(np.round(xp[i, 1]))
        if grid[xx, yy] == 1:
            count += 1
    
    if count == 3:
        return True, tuple(inBord)
    
    # Si plus d'un nœud est sur le bord, on considère l'élément comme "en mur"
    return False, tuple(inBord)


class MEF:
    def __init__(self, grid, nodes, connectivity):
        """
        grid : tableau 2D (ex. image) servant à définir le domaine
        nodes : tableau (n_nodes, 2) des coordonnées (normalisées dans [0,1] par exemple)
        connectivity : tableau (n_elements, 3) définissant les indices des nœuds pour chaque élément triangulaire (P1)
        """
        self.grid = grid
        self.nodes = nodes
        self.connectivity = connectivity
        self.n_elements = connectivity.shape[0]
        # Fonctions de base de référence pour l'élément P1 (triangle de référence)
        # Remarque : Ces expressions sont valables sur le triangle de référence (ex. (0,0), (1,0), (0,1))
        self.phi = {
            0: lambda x, y: 1 - x - y,
            1: lambda x, y: x,
            2: lambda x, y: y
        }
    
    def mass(self, element):
        """
        Calcule la matrice de masse élémentaire pour l'élément 'element' (P1).
        Utilise la formule consistante pour un triangle.
        """
        # Coordonnées des 3 nœuds de l'élément
        x = self.nodes[self.connectivity[element]]  # shape (3, 2)
        # Calcul de l'aire A = 0.5 * |det([x, y, 1])|
        A = 0.5 * np.abs(np.linalg.det(np.array([x[:, 0], x[:, 1], [1, 1, 1]])))
        # Matrice de masse consistante
        Me = A / 12 * np.array([[2, 1, 1],
                                [1, 2, 1],
                                [1, 1, 2]])
        return Me
    
    def rigidity(self, element):
        """
        Calcule la matrice de rigidité élémentaire (stiffness) pour l'élément 'element' (P1).
        On utilise les gradients des fonctions de base sur le triangle de référence.
        """
        x = self.nodes[self.connectivity[element]]  # shape (3,2)
        # Gradients dans l'élément de référence pour un triangle P1 :
        # phi_0 = 1 - x - y  => grad = (-1, -1)
        # phi_1 = x          => grad = (1, 0)
        # phi_2 = y          => grad = (0, 1)
        d_phi = np.array([[-1, -1], [1, 0], [0, 1]])
        
        # Jacobien de la transformation affine (2x2)
        Jacobian = d_phi.T.dot(x)
        detJ = np.linalg.det(Jacobian)
        # On prend l'inverse transpose pour transformer les gradients
        B = np.linalg.inv(Jacobian).T
        # Calcul de la matrice b = (Jacobian^{-1})(Jacobian^{-1})^T
        b = B.T.dot(B)
        
        Re = np.zeros((3, 3))
        # Formule : a_{ij} = (grad_phi_i)^T (Jacobian^{-1} Jacobian^{-T}) (grad_phi_j) * |det(Jacobian)|
        for i in range(3):
            for j in range(3):
                Re[i, j] = d_phi[i].dot(b).dot(d_phi[j])
        return Re * np.abs(detJ)
    
    def boundary(self, element, inBoundary, Nx, Ny, Lx, Ly):
        """
        Calcule la matrice d'intégration sur l'arête de l'élément 'element'
        qui est sur le bord (pour la condition de type Robin).
        
        inBoundary : liste des indices locaux (0, 1 ou 2) des nœuds de l'élément se trouvant sur le bord.
                     Pour un élément P1, s'il y a une arête sur le bord, inBoundary devrait contenir 2 indices.
        Nx, Ny : dimensions (nombre de points) de la grille, servant à convertir les coordonnées normalisées
                 en coordonnées physiques (si nécessaire).
        """
        # Initialisation de la matrice locale de bord (taille 3x3)
        Be = np.zeros((3, 3), dtype=complex)
        # Si aucun nœud de bord n'est spécifié ou moins de 2, on retourne zéro.
        if inBoundary is None or len(inBoundary) < 2:
            return Be
        
        # Récupération des coordonnées des nœuds de l'élément.
        # Ici, on suppose que self.nodes est défini dans [0,1].
        # On convertit en coordonnées physiques en multipliant par (Nx-1, Ny-1).
        xp = self.nodes[self.connectivity[element]] * np.array([(Nx-1)/Lx, (Ny-1)/Ly])

        
        # Pour une arête sur le bord, on suppose que inBoundary contient exactement 2 indices locaux.
        i_local, j_local = inBoundary[0], inBoundary[1]
        P = xp[i_local]  # Coordonnées du premier nœud de l'arête
        Q = xp[j_local]  # Coordonnées du second nœud de l'arête
        Le = np.linalg.norm(Q - P)  # Longueur de l'arête
        
        # Paramétrisation de l'arête sur l'intervalle de référence [-1, 1]
        # Les fonctions de forme linéaires sur l'arête sont :
        N1 = lambda xi: (1 - xi) / 2
        N2 = lambda xi: (1 + xi) / 2
        
        # Points et poids pour une quadrature de Gauss à 2 points sur [-1,1]
        gauss_points = np.array([-1/np.sqrt(3), 1/np.sqrt(3)])
        gauss_weights = np.array([1.0, 1.0])
        
        # Calcul de la matrice d'intégration sur l'arête (2x2) par quadrature de Gauss
        B_edge = np.zeros((2, 2))
        for xi, w in zip(gauss_points, gauss_weights):
            N1_val = N1(xi)
            N2_val = N2(xi)
            B_edge[0, 0] += w * (N1_val**2)
            B_edge[0, 1] += w * (N1_val * N2_val)
            B_edge[1, 0] += w * (N1_val * N2_val)
            B_edge[1, 1] += w * (N2_val**2)
        # Multiplication par le facteur de changement de variable ds = (Le/2) dξ
        B_edge *= (Le / 2)
        
        # Insertion de la contribution B_edge dans la matrice locale Be aux positions correspondant aux nœuds sur l'arête.
        Be[i_local, i_local] += B_edge[0, 0]
        Be[i_local, j_local] += B_edge[0, 1]
        Be[j_local, i_local] += B_edge[1, 0]
        Be[j_local, j_local] += B_edge[1, 1]
        
        return Be
    
    def external_force(self, element, f):
        """
        Calcule le vecteur de force externe élémentaire pour l'élément 'element'.
        fs est une fonction f(x,y).
        """
        x = self.nodes[self.connectivity[element]]
        # Coordonnées du centre de gravité de l'élément
        x_, y_ = x.sum(axis=0) / 3
        A = 0.5 * np.abs(np.linalg.det(np.array([x[:, 0], x[:, 1], [1, 1, 1]])))
        fe = A / 3 * np.array([1, 1, 1]) * f(x_, y_)
        return fe
    
    def assemblage(self, fs, k, Lx, Ly):
        """
        Assemble la matrice globale A et le vecteur global F pour le problème de Helmholtz.
        fs : fonction source f(x,y)
        k : nombre d'onde (fréquence)
        """
        n_nodes = self.nodes.shape[0]
        n_elements = self.connectivity.shape[0]
        
        # Initialisation de la matrice globale et du vecteur de force (complexes)
        A = np.zeros((n_nodes, n_nodes), dtype=complex)
        F = np.zeros(n_nodes, dtype=complex)
        
        Nx, Ny = self.grid.shape
        
        for element in tqdm(range(n_elements)):
            Me = self.mass(element)
            Re = self.rigidity(element)
            
            # La fonction isInWall doit retourner un booléen indiquant si l'élément est sur un mur
            # et inBord, la liste des indices locaux (0,1,2) des nœuds sur le bord.
            isInWall_, inBord = isInWall(self.grid, self.nodes, self.connectivity, element, Nx, Ny, Lx, Ly)
            
            # if isInWall_:
            #     print(f"Element {element} dans le mur")
            
            n = 4 if isInWall_ else 1
            
            Be = self.boundary(element, inBord, Nx, Ny, Lx, Ly)
            
            # Assemblage de la contribution de l'élément dans la matrice globale
            for i in range(3):
                for j in range(3):
                    A[self.connectivity[element, i], self.connectivity[element, j]] += (
                        -Re[i, j] + k**2 * n**2 * Me[i, j]) + 1j * k * n * Be[i, j]
                    
            fe = -self.external_force(element, fs)
            for i in range(3):
                F[self.connectivity[element, i]] += fe[i]
        return A, F