#This program will parse the coordinates of a molecule in the form of a gaussian input file
#and then create a rough shape of the molecule, then give coordinates to place point charges
#around the molecule as a predetermined distance

#input file given in cli must have the following format (same as Gaussian coordinate input):
#   atom1   x1  y1  z1
#   atom2   x2  y2  z2
#   atom3   x3  y3  z3
#   ...     ... ... ...

#Functions will include:
#createBox, inputs: size of environment, minimum distance to place charges from molecule, number of charges to build, center of molecule

#Future update to replace cartesian charge placement and distance calculation with radial calculations instead

import sys  #for command line input upon program execution
import random   #to generate random values
import matplotlib.pyplot as plt #to plot charge and atom positions
from mpl_toolkits import mplot3d    #for 3d plot of charge and atom positions

#Class to hold atom type and coordinates.
class atomIdent:
    def __chgCoord__(x_coord, y_coord, z_coord):
        return f"{x_coord} {y_coord} {z_coord}"
    def __molCoord__(species, x_coord, y_coord, z_coord):
        return f"{species} {x_coord} {y_coord} {z_coord}"


#Function to assign randomized charges
def createBox(envSize, internalBoundary, numberCharge, centerMol):
    chargeCounter = []
    k = 0
    charges = []
    while k < numberCharge:
        num = random.randint(1,8)
        if num == 1:
            x_coord = random.uniform(centerMol[0], centerMol[0]+envSize)
            y_coord = random.uniform(centerMol[1], centerMol[1]+envSize)
            z_coord = random.uniform(centerMol[2], centerMol[2]+envSize)
        if num == 2:
            x_coord = random.uniform(centerMol[0], centerMol[0]+envSize)
            y_coord = random.uniform(centerMol[1], centerMol[1]+envSize)
            z_coord = random.uniform(centerMol[2], centerMol[2]-envSize)
        if num == 3:
            x_coord = random.uniform(centerMol[0], centerMol[0]+envSize)
            y_coord = random.uniform(centerMol[1], centerMol[1]-envSize)
            z_coord = random.uniform(centerMol[2], centerMol[2]-envSize)
        if num == 4:
            x_coord = random.uniform(centerMol[0], centerMol[0]-envSize)
            y_coord = random.uniform(centerMol[1], centerMol[1]-envSize)
            z_coord = random.uniform(centerMol[2], centerMol[2]-envSize)
        if num == 5:
            x_coord = random.uniform(centerMol[0], centerMol[0]-envSize)
            y_coord = random.uniform(centerMol[1], centerMol[1]-envSize)
            z_coord = random.uniform(centerMol[2], centerMol[2]+envSize)
        if num == 6:
            x_coord = random.uniform(centerMol[0], centerMol[0]-envSize)
            y_coord = random.uniform(centerMol[1], centerMol[1]+envSize)
            z_coord = random.uniform(centerMol[2], centerMol[2]+envSize)
        if num == 7:
            x_coord = random.uniform(centerMol[0], centerMol[0]-envSize)
            y_coord = random.uniform(centerMol[1], centerMol[1]+envSize)
            z_coord = random.uniform(centerMol[2], centerMol[2]-envSize)
        if num == 8:
            x_coord = random.uniform(centerMol[0], centerMol[0]+envSize)
            y_coord = random.uniform(centerMol[1], centerMol[1]-envSize)
            z_coord = random.uniform(centerMol[2], centerMol[2]+envSize)
        file = open(sys.argv[1])
        for line in file.readlines():
            line = str(line)
            line = line.split()
            if ((x_coord-float(line[1]))**2 + (y_coord-float(line[2]))**2 + (z_coord-float(line[3]))**2)**(1/2) >= internalBoundary and k not in chargeCounter:
                charges.append(atomIdent.__chgCoord__(x_coord, y_coord, z_coord))
                chargeCounter.append(k)
                break
            k = k + 1
    return charges


###MAIN###

##Greeting and user options
#cite ''A. Bondi. "van der waals Volumes and Radii". J. Chem. Phys. (1964) 68, 441'' 
#Default internal boundary set to vdw radius of carbon, 1.75 A  (from A. Bondi), to prevent charge overlap with atomic orbital
#If atoms other than C, H, N, or O are used may need to enter alternative minimum distance
print("This program will provide randomized charge coordinates within a distance from your provided molecule\n")
internalBoundary = input("In Angstroms, how close to your molecule would you like like the cutoff radius to be? (default 1.75 angstroms): ")
if internalBoundary == (""):
    internalBoundary = 1.75
envSize = input("How many Angstroms from the center of your molecule would you like your environment to reach? (default 5 angstroms): ")
if envSize == (""):
    envSize = 5
envSize = float(envSize)
internalBoundary = float(internalBoundary)
numberCharge = input("How many point charges would you like calculated? (default 100): ")
if numberCharge == (''):
    numberCharge = 100
numberCharge = float(numberCharge)
chargeAmt = input("What is the charge of your point charges in e? (default -0.1e): ")
if chargeAmt == (""):
    chargeAmt = -0.1
chargeAmt = float(chargeAmt)

print("\n")
print("==========")
print("Your solvent bubble is ", envSize*2, "Angstroms across")
print("Point charges will be placed ", internalBoundary, " Angstroms from the closest atom")
print("==========")
print("\n")

#Allows automatic input of coordinate file and places coordinates into a list within the program to be called later
#sys.argv[] used for command line interface entry of text file to be parsed
#also takes the sum of each cordinate set (x, y, z) to be used later when finding center of molecule
print("\n")
print("==========")
print("Your file will now be parsed")
print("==========")
x_sum = 0
y_sum = 0
z_sum = 0
file = open(sys.argv[1], 'r')
saved_coord = []
for line in file.readlines():
    line = str(line)
    line = line.split()
    print(line)
    saved_coord.append(atomIdent.__molCoord__(line[0], line[1], line[2], line[3]))
    x_sum = x_sum + float(line[1])
    y_sum = y_sum + float(line[2])
    z_sum = z_sum + float(line[3])
print("==========")

#Outputs coordinates placed into list
print("\n")
print("==========")
j = 0
while j < len(saved_coord):
    print("Atom ", j+1, " is:")
    print(saved_coord[j])
    j = j+1
print("==========")
print("\n")

#to compute the center of the molecule to form our solvent environment
print("==========")
x_avg = x_sum/len(saved_coord)
y_avg = y_sum/len(saved_coord)
z_avg = z_sum/len(saved_coord)
centerMol = [x_avg, y_avg, z_avg]
print("The center of your solvent environment will be at: ", centerMol)
print("==========")

print("\n")
print("==========")
print("Outputting random charge coordinates:")
print("==========")
print("\n")

#calls function to build solvent environment , then output charges
print("==========")
charges = createBox(envSize, internalBoundary, numberCharge, centerMol)
n = 0
while n < len(charges):
    print("Charge coordinate ", n+1, " is")
    print(charges[n], chargeAmt)
    n = n+1
print("==========")
print("\n")


#Checks if randomly placed charge is within vdw radius of atom (too close)
#or beyond the confines of the solvent box defined earlier (too far)
#Distance computed as: [(x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2]^(1/2)
print("==========")
print("Checking to see if any charges are within ", internalBoundary, " Angstroms of your molecule")
print("==========")
print("\n")
chargeCoord_final = []
outOfBounds =[]
final_list = []
molAtom_count = 1
file = open(sys.argv[1], 'r')       #needed to be opened again. not sure why... Don't remove or below will not work
for line in file.readlines():
    chargeCount = 1
    line = str(line)
    line = line.split()
    for charge in charges:
            charge = charge.split()
            x1 = float(charge[0])
            x2 = float(line[1])
            y1 = float(charge[1])
            y2 = float(line[2])
            z1 = float(charge[2])
            z2 = float(line[3])
            dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**(1/2)
            if dist < internalBoundary and dist > envSize:      #checking to see if charge is out of bounds  
                outOfBounds.append(chargeCount)
            if chargeCount not in outOfBounds and chargeCount not in final_list:
                final_list.append(chargeCount)
                chargeCoord_final.append(atomIdent.__chgCoord__(charge[0], charge[1], charge[2]))                    
            print("Charge ", chargeCount, " is ", dist, " Angstoms from atom ", molAtom_count)
            chargeCount = chargeCount + 1
    molAtom_count = + molAtom_count + 1
print("==========")



#prints final coordinates to screen
print("\n")
print("==========")
print("Final coordinates for point charges are below")
print("==========")
print("\n")
m = 0
while m < len(chargeCoord_final):
    print(chargeCoord_final[m], chargeAmt)
    m = m+1
print("\n")
print("==========")
print("you have ", m, " charge coordinates left")
print("==========")
print("\n")

outfile = open("chargeCoord.out", "w")
num = 1
for charge in chargeCoord_final:
    outfile.write(str(num))
    outfile.write(" ")
    outfile.write(charge)
    outfile.write(" 1.0")
    outfile.write('\n')
    num = num + 1
outfile.close()


#Plotting molecule and charge positions
figure = plt.figure()
plot = plt.axes(projection='3d')
chg_x_data = []
chg_y_data = []
chg_z_data = []
charge_num = 1
for charge in chargeCoord_final:
    charge = charge.split()
    chg_x_data.append(float(charge[0]))
    chg_y_data.append(float(charge[1]))
    chg_z_data.append(float(charge[2]))
    print("Charge ", charge_num, " is at: ", charge[0], charge[1], charge[2])
    plot.text(float(charge[0]), float(charge[1]), float(charge[2]), charge_num)
    charge_num += 1
print("\n")
print("==========")
print("Coordinates of your molecule are:")
print("==========")
print("\n")
mol_x_data = []
mol_y_data = []
mol_z_data = []
file = open(sys.argv[1], 'r')       #needed to be opened again. not sure why... Don't remove or below will not work
for line in file.readlines():
    line = line.split()
    mol_x_data.append(float(line[1]))
    mol_y_data.append(float(line[2]))
    mol_z_data.append(float(line[3]))
    print(line[0], line[1], line[2], line[3])
    plot.text(float(line[1]), float(line[2]), float(line[3]), line[0])
plot.scatter3D(chg_x_data, chg_y_data, chg_z_data, cmap = 'greens')
plot.scatter3D(mol_x_data, mol_y_data, mol_z_data, cmap = 'reds')
plt.show()
