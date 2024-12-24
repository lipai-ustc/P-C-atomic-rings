Extra info for the following work:  
[ACS Appl. Mater. Interfaces 2022, 14, 18506](https://pubs.acs.org/doi/10.1021/acsami.2c01494)  
Title: Electrochemistry of P–C Bonds in Phosphorus–Carbon Based Anode Materials  
Authors: Pai Li, Zhenyu Li, _et. al._  
<img src="https://github.com/user-attachments/assets/a00dfe60-f0a4-4bf9-bc71-6584ac1dd4af" width="500">

### Strategies for counting atomic rings in a atomic structure model in `POSCAR` format
Note that you should first confine the cell, otherwise, you might overestimate the number of rings.  
Afterword: Consider all rings on the interface as periodic is an overestimation. A more accurate method is to calculate the number of rings including and excluding the interface, then take half the number of rings on the interface.
The following script generates an adjacency matrix (`NMat`) from a crystal structure defined in a `POSCAR` file. The script scales the unit cell to ensure sufficient separation between periodic images, removes certain types of atoms ('Na' and 'Li'), and creates a supercell. It then calculates the adjacency matrix based on interatomic distances, where atoms closer than 2.5 Å are considered connected. After saving this matrix, the script defines functions to perform a depth-first search (DFS) to count cycles of specified lengths within the graph represented by the adjacency matrix. Finally, it prints out the total number of cycles for lengths 3 through 8. 

```python
import ase
from ase.io import read
import numpy as np

# Read the atomic structure from the CONTCAR file
a = read("CONTCAR")

# Calculate the scaling factors for each cell vector to ensure a minimum distance of 16 Ångströms between periodic images
f0 = int(np.ceil(16 / np.linalg.norm(a.cell[0])))
f1 = int(np.ceil(16 / np.linalg.norm(a.cell[1])))
f2 = int(np.ceil(16 / np.linalg.norm(a.cell[2])))

# Remove atoms with symbols 'Na' and 'Li'
indices_to_remove = [atom.index for atom in a if atom.symbol in ['Na', 'Li']]
del a[indices_to_remove]

print("Cell vectors:")
print(a.cell)
print("Number of atoms:", len(a))
print("Scaling factors:", f0, f1, f2)

# Create a supercell
atoms = a * (f0, f1, f2)
print("Number of atoms in supercell:", len(atoms))

# Initialize the adjacency matrix with zeros
size = len(a) * f0 * f1 * f2
NMat = np.zeros((size, size))

# Populate the adjacency matrix based on interatomic distances
for i in range(len(atoms)):
    for j in range(i):
        if atoms.get_distance(i, j, mic=True) < 2.5:
            NMat[i, j] = 1
            NMat[j, i] = 1

# Save the adjacency matrix to a file
np.savetxt("NMat", NMat, fmt='%1d')

# Function to count the number of cycles using the adjacency matrix
# Python Program to count cycles of length n in a given graph.

# Number of vertices
graph = np.loadtxt("NMat")
V = graph.shape[0]

def DFS(graph, marked, n, vert, start, count):
    # Mark the vertex vert as visited
    marked[vert] = True
    
    # If the path of length (n-1) is found
    if n == 0:
        # Mark vert as unvisited to make it usable again
        marked[vert] = False
        
        # Check if vertex vert can end with vertex start
        if graph[vert][start] == 1:
            count += 1
            return count
        else:
            return count
    
    # For searching every possible path of length (n-1)
    for i in range(V):
        if not marked[i] and graph[vert][i] == 1:
            # DFS for searching path by decreasing length by 1
            count = DFS(graph, marked, n-1, i, start, count)
    
    # Marking vert as unvisited to make it usable again
    marked[vert] = False
    return count

# Counts cycles of length N in an undirected and connected graph
def countCycles(graph, n):
    # All vertices are marked unvisited initially
    marked = [False] * V
    
    # Searching for cycle by using v-n+1 vertices
    count = 0
    for i in range(V - (n - 1)):
        count = DFS(graph, marked, n-1, i, i, count)
        
        # ith vertex is marked as visited and will not be visited again
        marked[i] = True
    
    return int(count / 2)

# Main function to count cycles of various lengths
for n in [3, 4, 5, 6, 7, 8]:
    print(f"Total cycles of length {n} are {countCycles(graph, n)}")
```

