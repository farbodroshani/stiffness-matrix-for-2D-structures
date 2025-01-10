# Stiffness Matrix for 2D Structures
This repository contains Python scripts designed for the analysis of 2D truss structures and bending elements (frames) using the Direct Stiffness Matrix Method. The scripts calculate the displacements of each node and determine the strains and stresses in each element.
# Overview
The provided Python scripts analyze 2D structures by:

1-Calculating the displacement of each node.

2-Determining the strains and stresses in each element.

The method used is the Direct Stiffness Matrix Method.

# How to use 
# Truss Elements Analysis

This script analyzes 2D truss structures.

Instructions

1. Initialization and Input Data
Import the necessary libraries and set the print options for better readability.

2. Input Nodes and Elements
Input the total number of nodes and elements.

Collect the x and y coordinates of each node.

Input the modulus of elasticity of the material.

3. Define Elements
Input the start node, end node, and cross-sectional area for each element.

Calculate the length, constant, cosine, and sine for each element based on the node coordinates.

4. Construct Element Stiffness Matrices
Calculate and store the element stiffness matrices.

5. Global Stiffness Matrix Assembly
Map the element stiffness matrices to the global stiffness matrix.

6. Boundary Conditions and Loading
Define displacements and forces. Specify the support conditions and loading for each node.

7. Matrix Reduction
Reduce the global stiffness matrix and force matrix by removing rows and columns corresponding to the boundary conditions.

8. Solving the System
Calculate the displacements and forces for the nodes by solving the linear system.

9. Calculate New Coordinates and Lengths
Calculate the new coordinates of each node after displacement.

Calculate the new lengths of the elements.

# Bending Elements and Frames Analysis

This script analyzes 2D frame structures.

Instructions
1. Initialization and Input Data
Import the necessary libraries and set the print options for better readability.

2. Input Nodes and Elements
Input the total number of nodes and elements.

Collect the x and y coordinates of each node.

Input the modulus of elasticity of the material.

3. Define Elements
   
Input the start node, end node, cross-sectional area, and moment of inertia for each element.

Calculate the length, constant, cosine, and sine for each element based on the node coordinates.

4. Construct Element Stiffness Matrices
   
Calculate and store the local stiffness matrices for each element.

6. Transform Local Stiffness Matrices
   
Calculate the transformation matrices and transform the local stiffness matrices to the global coordinate system.

7. Global Stiffness Matrix Assembly
   
Integrate the transformed local stiffness matrices into the global stiffness matrix.

8. Boundary Conditions and Loading
   
Initialize displacement and force vectors.

Specify the forces and moments at each node.

Define the boundary conditions for each node (whether it is free to move or rotate).

8. Partition the Global Stiffness Matrix
   
Partition the global stiffness matrix and force vector into known and unknown displacements.

10. Solving the System
    
Solve for the unknown displacements.

Assign the solved displacements to the displacement vector.

10. Calculate Reactions and Element Forces
    
Calculate the reactions at each node with known displacements.

Calculate the forces in each element.

11. Calculate Strains and Stresses

Calculate the strains and stresses in each element.

# Output

The scripts will print the displacements and reactions for each node, as well as the strains and stresses for each element.
