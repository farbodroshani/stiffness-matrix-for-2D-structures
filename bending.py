import math
import numpy as np
import time

np.set_printoptions(linewidth=np.inf)
np.set_printoptions(precision=3, suppress=True)

print('Frame Analysis Code')



totalnodes = int(input('Enter the total number of nodes: '))  # total nodes
totalelements = int(input('Enter the total number of Elements: '))  # total elements

xco = []  # x co ordinate of nodes
yco = []  # y co ordinate of nodes

for i in range(totalnodes):
    x = float(input('Enter the x co-ordinate of node ' + str(i + 1) + ' in cm: '))
    y = float(input('Enter the y co-ordinate of node ' + str(i + 1) + ' in cm: '))
    xco.append(x)
    yco.append(y)

E = float(input('Enter the Modulus of Elasticity in ton/cm2: '))



elementarea = []  # area of the elements
snofel = []  # start node of elements
enofel = []  # end node of elements
lenofel = []  # length of the elements
elcon = []  # constant of the elements
cosofel = []  # cos of elements
sinofel = []  # sin of elements
elmomentofinteria = []  # moment of inertia of the elements

for i in range(totalelements):
    a = int(input('Enter the Start node of element ' + str(i + 1) + ': '))
    b = int(input('Enter the End node of element ' + str(i + 1) + ': '))
    area = float(input('Enter the Area of cross-section of the element ' + str(i + 1) + ' in cm2: '))
    momentofinteria = float(input('Enter the moment of inertia of element ' + str(i + 1) + ' in cm4: '))
    x1 = float(xco[a - 1])
    y1 = float(yco[a - 1])
    x2 = float(xco[b - 1])
    y2 = float(yco[b - 1])
    l = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    con = E / l
    cos = (x2 - x1) / l
    sin = (y2 - y1) / l

    snofel.append(a)
    enofel.append(b)
    lenofel.append(l)
    elcon.append(con)
    cosofel.append(cos)
    sinofel.append(sin)
    elementarea.append(area)
    elmomentofinteria.append(momentofinteria)



ellocalkmat = []  # local k matrix of elements

for i in range(totalelements):
    A = float(elementarea[i])
    I = float(elmomentofinteria[i])
    L = float(lenofel[i])

    mat = elcon[i] * np.array([[A, 0, -0, -A, 0, 0],
                               [0, 12 * I / L ** 2, 6 * I / L, 0, -12 * I / L ** 2, 6 * I / L],
                               [0, 6 * I / L, 4 * I, 0, -6 * I / L, 2 * I],
                               [-A, 0, 0, A, 0, 0],
                               [0, -12 * I / L ** 2, -6 * I / L, 0, 12 * I / L ** 2, -6 * I / L],
                               [0, 6 * I / L, 2 * I, 0, -6 * I / L, 4 * I]])
    ellocalkmat.append(mat)



eltransformmat = []

for i in range(totalelements):
    c = float(cosofel[i])
    s = float(sinofel[i])
    mat = np.array([[c, s, 0, 0, 0, 0],
                    [-s, c, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, c, s, 0],
                    [0, 0, 0, -s, c, 0],
                    [0, 0, 0, 0, 0, 1]])

    eltransformmat.append(mat)



globalkmat = np.zeros((3 * totalnodes, 3 * totalnodes))  # global stiffness matrix

for i in range(totalelements):
    kT = np.matmul(ellocalkmat[i], eltransformmat[i])
    Ttran = np.transpose(eltransformmat[i])
    globalmat = np.matmul(Ttran, kT)

    a = snofel[i] - 1
    b = enofel[i] - 1

    globalkmat[3 * a:3 * (a + 1), 3 * a:3 * (a + 1)] += globalmat[:3, :3]
    globalkmat[3 * a:3 * (a + 1), 3 * b:3 * (b + 1)] += globalmat[:3, 3:]
    globalkmat[3 * b:3 * (b + 1), 3 * a:3 * (a + 1)] += globalmat[3:, :3]
    globalkmat[3 * b:3 * (b + 1), 3 * b:3 * (b + 1)] += globalmat[3:, 3:]

print('\nGlobal Stiffness Matrix of the Frame\n')
print(globalkmat)




print('\n\n________________Boundry condition and Loading______________\n')

# Create a vector of nodal displacements and initialize it with zeros
u = np.zeros(3 * totalnodes)

# Create a vector of nodal forces and initialize it with zeros
f = np.zeros(3 * totalnodes)

# Loop over the nodes and assign the corresponding forces to the vector


for i in range(totalnodes):
    fx = float(input(f"Enter the force in the x-direction at node {i+1} in tons: "))
    fy = float(input(f"Enter the force in the y-direction at node {i+1} in tons: "))
    mz = float(input(f"Enter the moment at node {i+1} in ton.cm: "))
    f[3*i] = fx
    f[3*i + 1] = fy
    f[3*i + 2] = mz

# Create a list of indices for the known displacements

known_displacements_indices = []

# Loop over the boundary conditions and add the indices to the list

for i in range(totalnodes):
    x_free = input(f"Is node {i+1} free to move in the x-direction? (y/n): ")
    y_free = input(f"Is node {i+1} free to move in the y-direction? (y/n): ")
    rot_free = input(f"Is node {i+1} free to rotate? (y/n): ")
    if x_free.lower() == 'n':
        known_displacements_indices.append(3*i)
    if y_free.lower() == 'n':
        known_displacements_indices.append(3*i + 1)
    if rot_free.lower() == 'n':
        known_displacements_indices.append(3*i + 2)

# Create a list of indices for the unknown displacements

unknown_displacements_indices = [i for i in range(3*totalnodes) if i not in known_displacements_indices]

# Partition the global stiffness matrix and the forces vector

K_known = globalkmat[np.ix_(known_displacements_indices, known_displacements_indices)]
K_unknown = globalkmat[np.ix_(unknown_displacements_indices, unknown_displacements_indices)]
K_mixed = globalkmat[np.ix_(unknown_displacements_indices, known_displacements_indices)]

F_known = f[known_displacements_indices]
F_unknown = f[unknown_displacements_indices]

# Solve for the unknown displacements

u_unknown = np.linalg.solve(K_unknown, F_unknown - np.matmul(K_mixed, u[known_displacements_indices]))

# Assign the solved displacements to the u vector

u[unknown_displacements_indices] = u_unknown

# Calculate the reactions at each node
reactions = np.matmul(globalkmat, u)

# Print the displacements for each node and each axis

print('\n\nDisplacement of nodes in cm : \n')
for i in range(totalnodes):
    print(f"Node {i+1}:")
    print(f"  Displacement in x-direction: {u[3*i]} cm")
    print(f"  Displacement in y-direction: {u[3*i + 1]} cm")
    print(f"  Rotation: {u[3*i + 2]} rad")

# Print the reactions for each node with known displacements

print("\n\nReactions: \n")

for i in range(totalnodes):
    print(f"Node {i+1}:")
    print(f"  Reaction force in x-direction: {reactions[3*i]} ton")
    print(f"  Reaction force in y-direction: {reactions[3*i + 1]} ton")
    print(f"  Reaction moment: {reactions[3*i + 2]} ton.cm")


element_forces = []
for i in range(totalelements):
    a = snofel[i] - 1
    b = enofel[i] - 1
    u_element = np.array([u[3*a], u[3*a + 1], u[3*a + 2], u[3*b], u[3*b + 1], u[3*b + 2]])
    f_element = np.matmul(ellocalkmat[i], np.matmul(eltransformmat[i], u_element))
    element_forces.append(f_element)

# Calculate the strains and stresses in each element
strains = []
stresses = []
for i in range(totalelements):
    f_element = element_forces[i]
    A = elementarea[i]
    L = lenofel[i]
    strain = (f_element[3] - f_element[0]) / L
    stress = f_element[0] / A
    strains.append(strain)
    stresses.append(stress)

# Print the strains and stresses for each element
print('\n\nStrain and stress in the elements')
print('\n***Positive is Tensile\nNegetive is Compressive***\n')

for i in range(totalelements):
    print(f"Element {i+1}:")
    print(f"  Strain: {strains[i]}")
    print(f"  Stress: {stresses[i]} ton/cm2")



time.sleep(1200)





    
