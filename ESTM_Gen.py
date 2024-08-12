'''
This program will take the input of a molecule as a collection of atoms in 3d space and
create an electrostatic spectral tuning map 
All mean vdw radii of atoms taken from ''A. Bondi. "van der waals Volumes and Radii". J. Chem. Phys. (1964) 68, 3'' 

Input file must have the following format:
atom1   x1  y1  z1
atom2   x2  y2  z2
atom3   x3  y3  z3
...    ... ... ...

All common organic atoms are supported along with several less common: H, B, C, N, O, F, P, S, Cl, Br, I
'''

import sys                    #for cli file input              
import matplotlib.pyplot as plt 
import numpy as np
import random

#Classes used to hold charge coordinates
class coordinates:
    def __chgCoord__(x_coord, y_coord, z_coord):    #charge coordinates around each atom
        return f"{x_coord} {y_coord} {z_coord}"

#Creates surface at vdw radius around atom
def createBubble(vdw_rad):
    fullAngle = np.linspace(0, 2*np.pi, 25)
    halfAngle = np.linspace(0, np.pi, 25)
    x_bub = vdw_rad * np.outer(np.cos(fullAngle), np.sin(halfAngle))    #np.outer() calculates outer product of two vectors 
    y_bub = vdw_rad * np.outer(np.sin(fullAngle), np.sin(halfAngle))
    z_bub = vdw_rad * np.outer(np.ones(np.size(fullAngle)), np.cos(halfAngle))  #np.ones() gives an array of only the number one
    return x_bub, y_bub, z_bub

    
#generates random point around selected atom
def createCharge(vdw_rad):
    fullAngle = random.uniform(0, 2*np.pi)
    halfAngle = random.uniform(0, np.pi)
    x_bub = vdw_rad * np.outer(np.cos(fullAngle), np.sin(halfAngle))    
    y_bub = vdw_rad * np.outer(np.sin(fullAngle), np.sin(halfAngle))
    z_bub = vdw_rad * np.outer(np.ones(np.size(fullAngle)), np.cos(halfAngle))
    return x_bub, y_bub, z_bub

#provides vdw radius for each atomic species from A. Bondi's paper
def vdwRad_assign(atom, figureType):
    if atom == 'H':
        vdw_rad = 1.20
    elif atom == 'B':
        vdw_rad = 1.65      #value is deBroglie wavelength, and thus will be generally overestimated for smaller atoms. No anisometrically derived vdw radius available for this atom from reference
    elif atom == 'C':
        vdw_rad = 1.70
    elif atom == 'N':
        vdw_rad = 1.55
    elif atom == 'O':
        vdw_rad == 1.52
    elif atom == 'F':
        vdw_rad = 1.47
    elif atom == 'P':
        vdw_rad = 1.80
    elif atom == 'S':
        vdw_rad = 1.80
    elif atom == 'Cl':
        vdw_rad = 1.75
    elif atom == 'Br':
        vdw_rad = 1.85
    elif atom == 'I':
        vdw_rad = 1.98
    else:
        print("Atom type not supported")
        sys.exit(0)
    if figureType.lower() == 'bubble':
        return vdw_rad          #using 1 vdw radius for bubble plot
    elif figureType.lower() == 'charge':
        return 2*vdw_rad    #using 2 x vdw radii as that is the cutoff for where vdw forces dominate for charge calculations

#Compares the individual charge to each atom in the molecule to ensure that it does not intersect the vdw radii of other atoms
def checkDistance(molecule, charge, figureType):
    hold = []   #this list will be populated when a coordinate lays within the vdw radius of a neighboring atom
    for atom in molecule.readlines():
        vdw_rad = vdwRad_assign(atom[0], figureType)
        dist = np.sqrt((float(atom[1])-float(charge[0]))**2 + (float(atom[2])-float(charge[1]))**2 + (float(atom[3])-float(charge[2]))**2)
        print(dist)
        if dist < vdw_rad:
            hold.append(dist)
    return len(hold)    #when len(hold) > 0, the atom being computed is within the vdw radius of the neighboring atom and will be discarded

###MAIN###

figureType = input("Would you like a bubble plot of your molecule or would you like to calculate charge coordinates for spectral tuning? (BUBBLE/CHARGE): ")

#creates visualized bubbles around atoms at their vdw radii
if figureType.lower() == 'bubble':
    fig = plt.figure()
    atomBubblePlot = fig.add_subplot(projection='3d')
    #opens coordinate file, parses and assigns vdw radius based on atomic species
    molecule = open(sys.argv[1], 'r')
    for atom in molecule.readlines():
        atom = str(atom)
        atom = atom.split()
        vdw_rad = vdwRad_assign(atom[0], figureType)
        bubble = createBubble(vdw_rad)
        atomBubblePlot.plot_surface(float(atom[1])+bubble[0], float(atom[2])+bubble[1], float(atom[3])+bubble[2])
    atomBubblePlot.set_aspect('equal')  #will scale axes to make perfect spheres. Does not graphically represent spheres. Do not use for now

#creates random charges around each atom
elif figureType.lower() =='charge':
    fig = plt.figure()
    chargeCoordPlot = fig.add_subplot(projection = '3d')
    charges = []    #will hold each individual charge
    molecule = open(sys.argv[1], 'r')
    for atom in molecule.readlines():
        atom = str(atom)
        atom = atom.split()
        chargeCoordPlot.scatter3D(float(atom[1]), float(atom[2]), float(atom[3]), color = 'blue')   #to plot molecule in blue
        chargeCoordPlot.text(float(atom[1]), float(atom[2]), float(atom[3]), atom[0])
        charge_count = 1
        while charge_count <= 2:   #assigns three charges per atom
            vdw_rad = vdwRad_assign(atom[0], figureType)
            initialCharge = createCharge(vdw_rad)
            dist_chk = checkDistance(molecule, initialCharge, figureType)
            if int(dist_chk) == 0:    
                charge_indv = coordinates.__chgCoord__(float(atom[1])+float(initialCharge[0]), float(atom[2])+float(initialCharge[1]), float(atom[3])+float(initialCharge[2]))
                charges.append(charge_indv)
                charge_count = charge_count + 1
            else:
                print("invalid")
    #for plotting charge points in red and printing to output file
    outfile = open("chargeCoord.out", 'w')
    charge_label_count = 1
    for charge in charges:
        charge = charge.split()
        chargeCoordPlot.scatter3D(float(charge[0]), float(charge[1]), float(charge[2]), color = 'red') 
        chargeCoordPlot.text(float(charge[0]), float(charge[1]), float(charge[2]), charge_label_count)
        print(charge_label_count, charge[0], charge[1], charge[2])
        #for printing to outfile:
        outfile.write(str(charge_label_count))
        outfile.write(" ")
        outfile.write(str(charge[0]))
        outfile.write(" ")
        outfile.write(str(charge[1]))
        outfile.write(" ")
        outfile.write(str(charge[2])) 
        outfile.write(" ")
        outfile.write("0.1")       
        outfile.write("\n")
        charge_label_count = charge_label_count + 1       
    chargeCoordPlot.set_aspect('equal')

else:
    print("Invalid response")
    sys.exit(0)

plt.show()
