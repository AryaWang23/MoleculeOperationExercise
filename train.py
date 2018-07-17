"""
     This is a molecule operation training.
     Date: 2018/07/09
     Wrote by Qianqian Lu

"""


import math
import pymolecule as pm


def train1():

    """Train1 : delete last two residue of protein and save as new pdb file"""

    molecule = pm.Molecule()

    molecule.load_pdb2('5wa5.pdb', True)

    molecule.print_out_info_pdb()

    select_res1 = molecule.selection(residue="EDO")

    molecule.delete_atoms(select_res1)

    select_res2 = molecule.selection(residue="4K4")

    molecule.delete_atoms(select_res2)

    molecule.print_out_info_pdb()

    molecule.save_pdb2('5wa5new.pdb', True)


def train2Andtrain3():

    """Train2 : get center cooridinate of geometry and center cooridinate of mass """
    #center cooridinate of geometry
    file=open('5wa5.pdb','r')
    lines = file.readlines()
    file.close()

    molecule3 = pm.Molecule()
    molecule3.load_pdb2('5wa5.pdb', True)


    for t in range(0, len(lines)):
        line = lines[t]
        if len(line) >= 7:
            if line[0:7] == "CRYST1 ":
                latticeA = float(line[8:16].strip())
                latticeB = float(line[17:25].strip())
                latticeC = float(line[26:34].strip())
                alpha = math.radians(float(line[35:41].strip()))
                belta = math.radians(float(line[42:48].strip()))
                gamma = math.radians(float(line[49:55].strip()))
    #get lattice vector from lattice constants and angles
    Ax = latticeA
    Ay = 0
    Az = 0
    Bx = latticeB*math.cos(gamma)
    By = latticeB*math.sin(gamma)
    Bz = 0
    Cx = latticeC*math.cos(belta)
    Cy = latticeC*((math.cos(alpha)-math.cos(belta)*math.cos(gamma))/math.sin(gamma))
    Cz = latticeC*(math.sqrt(1+2*math.cos(alpha)*math.cos(belta)*math.cos(gamma)-math.cos(alpha)**2 \
         -math.cos(belta)**2-math.cos(gamma)**2)/math.sin(gamma))
    centerX = (Ax+Bx+Cx)/2
    centerY = (Ay+By+Cy)/2
    centerZ = (Az+Bz+Cz)/2
    centerPoint =pm.Point(centerX,centerY,centerZ)
    centerPoint.print_coors()

    #center cooridinate of mass
    weightVectorSumX= 0
    weightVectorSumY= 0
    weightVectorSumZ= 0
    weightSum = 0
    for n in range(0, len(lines)):
        line = lines[n]

        if len(line) >= 7:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":  # Load atom data (coordinates, etc.)
                temp_atom = pm.Atom()
                temp_atom.read_pdb_line(line)
                weightVectorX=temp_atom.atomWeight()*temp_atom.coordinates.x
                weightVectorSumX += weightVectorX
                weightVectorY = temp_atom.atomWeight() * temp_atom.coordinates.y
                weightVectorSumY += weightVectorY
                weightVectorZ = temp_atom.atomWeight() * temp_atom.coordinates.z
                weightVectorSumZ += weightVectorZ
                weightSum += temp_atom.atomWeight()

    weightX = weightVectorSumX/weightSum
    weightY = weightVectorSumY/weightSum
    weightZ = weightVectorSumZ/weightSum

    weightPoint = pm.Point(weightX,weightY,weightZ)
    weightPoint.print_coors()

    """Train3 : calculate the distance between the first and last atom, and angles with any other atoms"""


    for m in range(0, len(lines)):
        line = lines[m]

        if len(line) >= 7:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":  # Load atom data (coordinates, etc.)
                atom = pm.Atom()
                atom.read_pdb_line(line)
                if atom.molecule_index == '1':
                    atomA = atom.coordinates
                if atom.molecule_index == str(len(molecule3.all_atoms)+1):
                    atomB = atom.coordinates
    distanceAB = atomA.distance_to_another_point(atomB)
    print distanceAB

    angleList = []
    for m in range(0, len(lines)):
        line = lines[m]
        if len(line) >= 7:
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":  # Load atom data (coordinates, etc.)
                atom = pm.Atom()
                atom.read_pdb_line(line)
                if atom.molecule_index != '1' and atom.molecule_index != str(len(molecule3.all_atoms)+1):
                    atomC = atom.coordinates
                    Angle = atomC.angle_between_three_points(point2=atomA,point3=atomB)
                    Angle = math.degrees(Angle)
                    angleList.append(Angle)
    print angleList

def train4():
    """Train4 : translate and rotate the molecule"""
    #translate
    molecule1 = pm.Molecule()
    molecule1.load_pdb2('5wa5.pdb', True)
    randomPoint = pm.Point(x=20,y=20,z=20)
    molecule1.translate_molecule(randomPoint)
    molecule1.print_out_info_pdb()
    molecule1.save_pdb2('5wa5TranslationNew.pdb', True)

    #rotate
    molecule1.undo()
    molecule1.rotate_molecule_around_pivot(1,1.0,1.0,1.0)
    molecule1.print_out_info_pdb()
    molecule1.save_pdb2('5wa5RotationNew.pdb', True)

def train5():

    """Train5 : Added new member fuction in Atom class of pymolecule.py: read_mol2_line(), create_mol2_line(),
                   Molecule class of pymolecule.py: loadmol2(), savemol2(), print_out_info_mol2() to finish load mol2 format file"""
    molecule2 = pm.Molecule()

    molecule2.loadmol2("00000001.mol2", False)

    molecule2.print_out_info_mol2()

    molecule2.save_mol2("00000001save.mol2", False)



if __name__ == '__main__':

    #delete last two residue of protein and save as new pdb file
    train1()
    #get center cooridinate of geometry and center cooridinate of mass, calculate the distance and angle
    train2Andtrain3()
    #translate and rotate the molecule
    train4()
    #load mol2 file
    train5()




















