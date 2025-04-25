import random
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d


###Classes###

class Atom:
  def __init__(self, species, x, y, z):
    self.x_coord = float(x)
    self.y_coord = float(y)
    self.z_coord = float(z)
    self.element = species
    self.getProperties()

  def atomMove(self, xM, yM, zM):
    self.x_coord += xM
    self.y_coord += yM
    self.z_coord += zM

  def ident(self):
    print(self.element, f"{self.x_coord:6f}", f"{self.y_coord:6f}", f"{self.z_coord:6f}", self.force, self.acceleration, self.velocity)

  def getProperties(self):
    match self.element:
      case "H":
        self.charge = 1
        self.mass = 1.008
        self.electronCount = 1
      case"C":
        self.charge = 6
        self.mass = 12.011
        self.electronCount = 4
      case "N":
        self.charge = 5
        self.mass = 14.007
        self.electronCount = 3
      case "O":
        self.charge = -2
        self.mass = 15.999
        self.electronCount = 2

  #basic properties
  element = " "
  x_coord = 0
  y_coord = 0
  z_coord = 0

  #element specific properties
  charge = 0
  mass = 0
  electronCount = 0 #valence + core
  coreCount = 0
  valenceCount = 0

  #MD required properties
  acceleration = 0
  velocity = 0
  force = 0

  pot_energy = 0
  kin_energy = 0
  total_energy = 0

###Functions###

#Function to randomly generate atoms within the cell
def generateAtoms(cellSize, atomSpecies_1, atomSpecies_2="", atomSpecies_3="", atomSpecies_4=""):
  cellVolume = cellSize[0] * cellSize[1] * cellSize[2]
  #cellDensity =

#Read a text file with each line containing species, x, y, z coordinates
def readInputFile():
  input = open(sys.argv[1], 'r')
  atomHold = []
  for atom in input.readlines():
    atom = str(atom)
    atom = atom.split()
    atomHold.append([atom[0], atom[1], atom[2], atom[3]])
  for atom in atomHold:
    atom = Atom(atom[0], atom[1], atom[2], atom[3])
    atom.ident()


#Molecular Dynamics
  #1) Provide initial conditions position, velocity, time, and compute potential energy
  #2) Compute forces: force = -(change in potential / change in position)
  #3) New configuration of particle from Newton's equations of motion (acceleration = force / mass)
  #4) iterate over all particles

#Use Verlet algorithm to integrate equations of motions

def runMD(atoms, dt, iterations = 5):
  count = 1
  while count <= iterations:
    print("Beginning MD run ", count)
    for i in range(len(atoms)):
      for j in range(len(atoms)):
        if i != j:
          calcForces(atoms[i], atoms[j], dt)
          calcAccelerations(atoms[i], atoms[j], dt)
          calcVelocities(atoms[i], atoms[j], dt)
          updatePositions(atoms[i], atoms[j], dt)
    count += 1
    for atom in atoms:
      atom.ident()
    plot(atoms)
    print()


def calcForces(atom1, atom2, dt):
  k = 8.99E9 #Newtons * meter^2 * coloumb^2
  electronChargeVal = 1.602176634E-19 #coloumbs
  distance = ((atom1.x_coord - atom2.x_coord)**2 + (atom1.y_coord - atom2.y_coord)**2 + (atom1.z_coord - atom2.z_coord)**2)**(1/2)
  force = k * (atom1.charge * atom2.charge * electronChargeVal**2) / distance
  atom1.force = force
  atom2.force = force


def calcAccelerations(atom1, atom2, dt):
  atom1.acceleration = atom1.force / atom1.mass
  atom2.acceleration = atom2.force / atom2.mass


def calcVelocities(atom1, atom2, dt):
  atom1.velocity = atom1.velocity + atom1.acceleration * dt
  atom2.velocity = atom2.velocity + atom2.acceleration * dt


def updatePositions(atom1, atom2, dt):
  #for atom 1
  atom1.x_coord = atom1.x_coord + atom1.velocity * dt + (atom1.acceleration / 2) * dt**2
  atom1.y_coord = atom1.y_coord + atom1.velocity * dt + (atom1.acceleration / 2) * dt**2
  atom1.z_coord = atom1.z_coord + atom1.velocity * dt + (atom1.acceleration / 2) * dt**2
  #for atom 2
  atom2.x_coord = atom2.x_coord + atom2.velocity * dt + (atom2.acceleration / 2) * dt**2
  atom2.y_coord = atom2.y_coord + atom2.velocity * dt + (atom2.acceleration / 2) * dt**2
  atom2.z_coord = atom2.z_coord + atom2.velocity * dt + (atom2.acceleration / 2) * dt**2


def plot(atoms):
  fgiure = plt.figure()
  plot = plt.axes(projection = '3d')
  x_coords = []
  y_coords = []
  z_coords = []
  for atom in atoms:
    x_coords.append(atom.x_coord)
    y_coords.append(atom.y_coord)
    z_coords.append(atom.z_coord)
    plot.text(float(atom.x_coord), float(atom.y_coord), float(atom.z_coord), atom.element)
  plot.scatter3D(x_coords, y_coords, z_coords)
  plt.show()


###For Testing on Work Computers###
def testRead():
  input = np.array([
    "H .01 .04 .05",
    "H .01 .06 .02",
    "C .01 .02 .02",
    "O .01 .01 0"
  ])
  atomHold = []
  for atom in input:
    atom = str(atom)
    atom = atom.split()
    atomHold.append(Atom(atom[0], atom[1], atom[2], atom[3]))
  return atomHold



###main()###

cellSize = [100, 100, 10]

atom1 = Atom("O", 0, 0, 0)
atom1.atomMove(random.uniform(0,float(cellSize[0])), random.uniform(0,float(cellSize[1])), random.uniform(0,float(cellSize[2])))
atom1.ident()

print("testing read function")
atomCoord = testRead()
runMD(atomCoord, 5)

