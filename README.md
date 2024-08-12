This repository contains three files for use in generating and analyzing point charges placed along the van der Waals surface of a molecule.

ESTM_Gen.py will parse a molecule and place point charges just outside the vdW radius of each atom, returning the coordinates for use in point charge calculations for an electrostatic tuning map. This program can also display the vdW surface around the molecule.

chargebox.py is an older version of ESTM_Gen.py, and only produces the point charges along the vdW surface

chargeBox_distance.py can be used to check the distance between the point charges produced in either chargebox.py or ESTM_Gen.py to verify that the point charges are beyond the vdW radii of the constituent atoms


Each file contains instructions on their use.
