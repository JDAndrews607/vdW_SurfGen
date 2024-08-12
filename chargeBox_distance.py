#This program is designed to be used in conjuction with chargebox.py
#This program will parse two files using command line arguments and compare the distance of charge coordinates to atom coordinates
#File1 will contain the molecule
#File2 will contain the charge coordinates for comparison

import sys
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


###MAIN##

#To open molecule and charge coordinate files
File1 = open(sys.argv[1], 'r')  #this will be the molecular coordinates
File2 = open(sys.argv[2], 'r')  #this will be the charge coordinates
mol_coord = File1.readlines()
chg_coord = File2.readlines()

print("===========================")
#computing distance
dist_array = []
for charge in chg_coord:
    charge = str(charge)
    charge = charge.split()
    for coord in mol_coord:
        coord = str(coord)
        coord = coord.split()
        x1 = float(coord[1])
        y1 = float(coord[2])
        z1 = float(coord[3])
        x2 = float(charge[1])
        y2 = float(charge[2])
        z2 = float(charge[3])
        dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**(1/2)
        print("Distance between charge ", charge[0], " and atom ", coord[0], " is ", dist, " Angstroms")
print("===========================")
print("\n")

#plotting molecule and charges on 3d scatterplot
figure = plt.figure()
plot = plt.axes(projection='3d')
chg_x_data = []
chg_y_data = []
chg_z_data = []
mol_x_data = []
mol_y_data = []
mol_z_data = []
for charge in chg_coord:
    charge = charge.split()
    chg_x_data.append(float(charge[1]))
    chg_y_data.append(float(charge[2]))
    chg_z_data.append(float(charge[3]))
    plot.text(float(charge[1]), float(charge[2]), float(charge[3]), charge[0])
for coord in mol_coord:
    coord = coord.split()
    chg_x_data.append(float(coord[1]))
    chg_y_data.append(float(coord[2]))
    chg_z_data.append(float(coord[3]))
    plot.text(float(coord[1]), float(coord[2]), float(coord[3]), coord[0])
plot.scatter3D(chg_x_data, chg_y_data, chg_z_data, color = 'blue')
plot.scatter3D(mol_x_data, mol_y_data, mol_z_data, color = 'red')
plot.set_aspect('equal')
plt.show()
