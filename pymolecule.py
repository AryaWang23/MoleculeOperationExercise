"""pymolecule is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pymolecule is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2011 Jacob D. Durrant. If you have any questions, comments, or
suggestions, please don't hesitate to contact me at jdurrant [at] ucsd [dot] edu.

The latest version of pymolecule can be downloaded from 
http://sourceforge.net/projects/autoclickchem/

If you use pymolecule in your work, please cite [REFERENCE HERE]"""

import math
import os

version="1.0.1"


class Point:
    """A class describing a single point or vector."""
    
    x=99999.0
    y=99999.0
    z=99999.0
    
    def __init__ (self, x, y ,z):
        """Initializes the Point class."""
        
        self.x = x
        self.y = y
        self.z = z
        
    def __str__(self): return self.print_coors()

    def print_coors(self): # kept for backwards comparibility.
        """Prints the coordinates of the point."""
        
        print str(self.x)+"\t"+str(self.y)+"\t"+str(self.z)
    
    def copy_of(self):
        """Returns an identical copy of this point (pymolecule.Point)."""
        
        return Point(self.x, self.y, self.z)
    
    def magnitude(self):
        """Calculates the distance (float) between this point (pymolecule.Point) and the origin. If the current class is being used as a vector, returns the length (float) of the vector."""

        return self.distance_to_another_point(Point(0,0,0))
        
    def distance_to_another_point(self, other_point):
        """Calculates the distance (float) to another point (pymolecule.Point)."""
        
        deltax = self.x - other_point.x
        deltay = self.y - other_point.y
        deltaz = self.z - other_point.z
        
        return math.sqrt(math.pow(deltax,2) + math.pow(deltay,2) + math.pow(deltaz,2))
        
    def dot_product(self, other_point):
        """Calculates the dot product (float) of two vectors (pymolecule.Point)."""
        
        return self.x * other_point.x + self.y * other_point.y + self.z * other_point.z
    
    def angle_with_other_vector(self, other_point):
        """Calculates the angle (float) between two vectors (pymolecule.Point), in radians."""
        
        new_self = self.normalized()
        new_other_point =other_point.normalized()
        dot_prod = new_self.dot_product(new_other_point)
        if dot_prod > 1.0: dot_prod = 1.0 # to prevent errors that can rarely occur
        if dot_prod < -1.0: dot_prod = -1.0
        return math.acos(dot_prod)
    
    def cross_product(self, other_point):
        """Calculates the cross product (pymolecule.Point) of two vectors (pymolecule.Point)."""
        
        Response = Point(0, 0, 0)
    
        Response.x = self.y * other_point.z - self.z * other_point.y
        Response.y = self.z * other_point.x - self.x * other_point.z
        Response.z = self.x * other_point.y - self.y * other_point.x
    
        return Response
    
    def normalized(self):
        """Returns a normalized copy (pymolecule.Point) of the specified vector (pymolecule.Point)."""
        
        dist = self.distance_to_another_point(Point(0,0,0))
        return Point(self.x/dist, self.y/dist, self.z/dist)
    
    def scalar_multiply(self, scalar):
        """Returns a scaled copy (pymolecule.Point) of this point, scaled by the specified amount (float)."""
        
        return Point(self.x * scalar, self.y * scalar, self.z * scalar)
    
    def subtract(self, other_point): # self - other_point
        """Returns the vector (pymolecule.Point) given by subtracting the specified vector (pymolecule.Point) from this one."""
        
        return Point(self.x - other_point.x, self.y - other_point.y, self.z - other_point.z)
    
    def add(self, other_point):
        """Returns the vector (pymolecule.Point) given by adding the specified vector (pymolecule.Point) to this one."""
        
        return Point(self.x + other_point.x, self.y + other_point.y, self.z + other_point.z)

    def dihedral(self, point2, point3, point4):
        """Calculates  formed by four points (pymolecule.Point), where the first point is self."""
        
        point1 = self
    
        b1 = point2.subtract(point1)
        b2 = point3.subtract(point2)
        b3 = point4.subtract(point3)
    
        b2Xb3 = b2.cross_product(b3)
        b1Xb2 = b1.cross_product(b2)
        
        b1XMagb2 = b1.scalar_multiply(b2.magnitude())
        radians = math.atan2(b1XMagb2.dot_product(b2Xb3), b1Xb2.dot_product(b2Xb3))
        return radians

    def angle_between_three_points(self, point2, point3): # As in three connected atoms
        """Computes the angle formed by three points (pymolecule.Point), where the first point is self."""
        
        point1 = self
        
        vector1 = point1.subtract(point2)
        vector2 = point3.subtract(point2)
        return vector1.angle_with_other_vector(vector2)

    def planarity(self, point2, point3, point4):
        """Determines whether four points (pymolecule.Point) lie in a common plane, where the first point is self.
        
        Arguments:
        point2 -- A point.
        point3 -- A point.
        point4 -- A point.
        
        Returns:
        A float, the minimum distance between one point and the plane formed by the other three.
        
        """

        point1 = self

        x1 = point1.x
        y1 = point1.y
        z1 = point1.z
        x2 = point2.x
        y2 = point2.y
        z2 = point2.z
        x3 = point3.x
        y3 = point3.y
        z3 = point3.z
        x4 = point4.x
        y4 = point4.y
        z4 = point4.z
    
        A = (y1*(z2-z3))+(y2*(z3-z1))+(y3*(z1-z2))
        B = (z1*(x2-x3))+(z2*(x3-x1))+(z3*(x1-x2))
        C = (x1*(y2-y3))+(x2*(y3-y1))+(x3*(y1-y2))
        D = ((-x1)*((y2*z3)-(y3*z2)))+((-x2)*((y3*z1)-(y1*z3)))+((-x3)*((y1*z2)-(y2*z1)))
        distance=(math.fabs((A*x4)+(B*y4)+(C*z4)+D))/(math.sqrt(math.pow(A,2) + math.pow(B,2) + math.pow(C,2)))
        
        A1 = (y1*(z2-z4))+(y2*(z4-z1))+(y4*(z1-z2))
        B1 = (z1*(x2-x4))+(z2*(x4-x1))+(z4*(x1-x2))
        C1 = (x1*(y2-y4))+(x2*(y4-y1))+(x4*(y1-y2))
        D1 = ((-x1)*((y2*z4)-(y4*z2)))+((-x2)*((y4*z1)-(y1*z4)))+((-x4)*((y1*z2)-(y2*z1)))
        distance1=(math.fabs((A1*x3)+(B1*y3)+(C1*z3)+D1))/(math.sqrt(math.pow(A1,2) + math.pow(B1,2) + math.pow(C1,2)))
        
        A2 = (y1*(z4-z3))+(y4*(z3-z1))+(y3*(z1-z4))
        B2 = (z1*(x4-x3))+(z4*(x3-x1))+(z3*(x1-x4))
        C2 = (x1*(y4-y3))+(x4*(y3-y1))+(x3*(y1-y4))
        D2 = ((-x1)*((y4*z3)-(y3*z4)))+((-x4)*((y3*z1)-(y1*z3)))+((-x3)*((y1*z4)-(y4*z1)))
        distance2=(math.fabs((A2*x2)+(B2*y2)+(C2*z2)+D2))/(math.sqrt(math.pow(A2,2) + math.pow(B2,2) + math.pow(C2,2)))
        
        A3 = (y4*(z2-z3))+(y2*(z3-z4))+(y3*(z4-z2))
        B3 = (z4*(x2-x3))+(z2*(x3-x4))+(z3*(x4-x2))
        C3 = (x4*(y2-y3))+(x2*(y3-y4))+(x3*(y4-y2))
        D3 = ((-x4)*((y2*z3)-(y3*z2)))+((-x2)*((y3*z4)-(y4*z3)))+((-x3)*((y4*z2)-(y2*z4)))
        distance3=(math.fabs((A3*x1)+(B3*y1)+(C3*z1)+D3))/(math.sqrt(math.pow(A3,2) + math.pow(B3,2) + math.pow(C3,2)))
    
        final_dist = -1
    
        if (distance < distance1 and distance < distance2 and distance < distance3):
            final_dist = distance
        elif (distance1 < distance and distance1 < distance2 and distance1 < distance3):
            final_dist = distance1
        elif (distance2 < distance and distance2 < distance1 and distance2 < distance3):
            final_dist = distance2
        elif (distance3 < distance and distance3 < distance1 and distance3 < distance2):
            final_dist = distance3
    
        # Now normalize by the length of the longest bond
    
        return final_dist

    def create_pdb_line(self, index):
        """Create a string formatted according to the PDB standard from this point class.
        
        Arguments:
        index -- An int, the atomic index to be used in the string.
        
        Returns:
        A string, formatted according to the PDB standard.
        
        """
        
        output = "ATOM "
        output = output + str(index).rjust(6) + "XX".rjust(5) + "XXX".rjust(4) + "X".rjust(2) + str(1).rjust(4)
        output = output + ("%.3f" % self.x).rjust(12)
        output = output + ("%.3f" % self.y).rjust(8)
        output = output + ("%.3f" % self.z).rjust(8)
        output = output + "X".rjust(24) # + "   " + str(uniqueID) #This last part must be removed
        return output

class Atom:
    """A class containing atomic information."""
    
    def __init__ (self):
        """Initializes the Atom class."""
        
        self.atomname = ""
        self.residue = ""
        self.coordinates = Point(99999, 99999, 99999)
        self.undo_coordinates = Point(99999, 99999, 99999)
        self.element = ""
        self.molecule_index = ""
        self.line=""
        self.chain = ""
        self.resname = ""
        self.resid = 0
        self.indecies_of_atoms_connecting=[]
        self.charge = 0
        
    def copy_of(self):
        """Returns an identical copy (pymolecule.Atom) of this atom."""
        
        newatom = Atom()
        newatom.atomname = self.atomname
        newatom.residue = self.residue
        newatom.coordinates = self.coordinates.copy_of()
        newatom.undo_coordinates = self.undo_coordinates.copy_of()
        newatom.element = self.element
        newatom.molecule_index = self.molecule_index
        newatom.line = self.line
        newatom.chain = self.chain
        newatom.resname = self.resname
        newatom.resid = self.resid
        for index in self.indecies_of_atoms_connecting:
            newatom.indecies_of_atoms_connecting.append(index)
            
        return newatom
    
    def string_id(self):
        """Returns a string describing the current atom (chain, resname, resid, atomname, index)."""
        
        toreturn = ""
        if self.chain.strip() != '': toreturn = toreturn + self.chain.strip() + ":"
        toreturn = toreturn + self.residue.strip() + '(' + str(self.resid) + "):" + self.atomname.strip() + '(' + str(self.PDBIndex) + ')'
        return toreturn


    def change_residue(self, NewResName):
        """Changes the name of the residue (string) associated with this atom."""
        
        self.residue = NewResName
        self.resname = NewResName
    
    def number_of_neighbors(self):
        """Returns the number of neighbor atoms bonded to this atom (int)."""
        
        return len(self.indecies_of_atoms_connecting)
    
    def set_undo_point(self):
        """Sets ("saves") the undo point of all atoms. Any subsequent manipulations of atomic coordinates can be "undone" by reseting to this configuration."""
        
        self.undo_coordinates = self.coordinates.copy_of()

    def undo(self):
        """Resets the coordinates of this atom to those saved using the set_undo_point function."""

        self.coordinates = self.undo_coordinates.copy_of()
        
    def read_pdb_line(self, Line):
        """Reads atomic information from a string formatted according to the PDB standard.
        
        Arguments:
        Line -- A string formatted according to the PDB standard.
        
        """
        
        self.line = Line
        self.atomname = Line[11:16].strip()
        
        if len(self.atomname)==1:
            self.atomname = self.atomname + "  "
        elif len(self.atomname)==2:
            self.atomname = self.atomname + " "
        elif len(self.atomname)==3:
            self.atomname = self.atomname + " " # This line is necessary for babel to work, though many PDBs in the PDB would have this line commented out
            
        # now get the chain
        self.chain = Line[21:22]
        
        # now get the resid
        try:
            self.resid = int(Line[22:26])
        except: 
            self.resid = 0
        # print Line[30:38].strip()
        # print len(Line[30:38].strip())
        # print float(Line[30:38].strip())
        self.coordinates = Point(float(Line[30:38].strip()), float(Line[38:46].strip()), float(Line[46:54].strip()))
        self.undo_coordinates = self.coordinates.copy_of() # so you can set the PDB coordinates to some original value
        
        if len(Line) >= 79:
            self.element = Line[76:79].strip().upper() # element specified explicitly at end of life
        elif self.element == "": # try to guess at element from name
            two_letters = self.atomname[0:2].strip().upper()
            if two_letters=='BR':
                self.element='BR'
            elif two_letters=='CL':
                self.element='CL'
            elif two_letters=='BI':
                self.element='BI'
            elif two_letters=='AS':
                self.element='AS'
            elif two_letters=='AG':
                self.element='AG'
            elif two_letters=='LI':
                self.element='LI'
            elif two_letters=='HG':
                self.element='HG'
            elif two_letters=='MG':
                self.element='MG'
            elif two_letters=='RH':
                self.element='RH'
            elif two_letters=='ZN':
                self.element='ZN'
            else: #So, just assume it's the first letter.
                self.element = self.atomname[0:1].strip().upper()
                
        # Any number needs to be removed from the element name
        self.element = self.element.replace('0','')
        self.element = self.element.replace('1','')
        self.element = self.element.replace('2','')
        self.element = self.element.replace('3','')
        self.element = self.element.replace('4','')
        self.element = self.element.replace('5','')
        self.element = self.element.replace('6','')
        self.element = self.element.replace('7','')
        self.element = self.element.replace('8','')
        self.element = self.element.replace('9','')

        self.molecule_index = Line[6:12].strip()
        self.residue = Line[16:20]
        if self.residue.strip() == "": self.residue = " MOL"
        self.resname = self.residue

    def __str__(self): return self.create_pdb_line(0)

    def read_mol2_line(self, Line):

        """Reads atomic information from a string formatted according to the mol2 standard.

        Arguments:
        Line -- A string formatted according to the PDB standard.

        """
        self.line = Line
        self.molecule_index = Line[0:6].strip()
        self.atomname = Line[7:13].strip()

        if len(self.atomname) == 1:
            self.atomname = self.atomname + "   "
        elif len(self.atomname) == 2:
            self.atomname = self.atomname + "  "
        elif len(self.atomname) == 3:
            self.atomname = self.atomname + " "


        self.coordinates = Point(float(Line[13:22]), float(Line[22:31]), float(Line[31:40]))
        self.undo_coordinates = self.coordinates.copy_of()
        self.element = Line[40:44].strip()

        self.resid = int(Line[50:55])
        self.residue = Line[55:60].strip()
        self.charge = float(Line[62:69])








    def create_pdb_line(self, index):
        """Create a string formatted according to the PDB standard from the atomic information contained in this atom class.
        
        Arguments:
        index -- An int, the atomic index to be used in the string.
        
        Returns:
        A string, formatted according to the PDB standard.
        
        """
        
        output = "ATOM "
        output = output + str(index).rjust(6) + self.atomname.rjust(5) + self.residue.rjust(4) + self.chain.rjust(2) + str(self.resid).rjust(4)
        output = output + ("%.3f" % self.coordinates.x).rjust(12)
        output = output + ("%.3f" % self.coordinates.y).rjust(8)
        output = output + ("%.3f" % self.coordinates.z).rjust(8)
        output = output + self.element.rjust(24) # + "   " + str(uniqueID) #This last part must be removed
        return output

    def create_mol2_line(self, index):
        """Create a string formatted according to the mol2 standard from the atomic information contained in this atom class.

        Arguments:
        index -- An int, the atomic index to be used in the string.

        Returns:
        A string, formatted according to the PDB standard.

        """

        # output = "@<TRIPOS>ATOM\n"
        output = str(index).rjust(6) + " " + self.atomname.ljust(6)
        output = output + ("%.4f" % self.coordinates.x).rjust(9)
        output = output + ("%.4f" % self.coordinates.y).rjust(9)
        output = output + ("%.4f" % self.coordinates.z).rjust(9)
        output = output +"  "+ self.element.ljust(7) + str(self.resid).rjust(4) + " " + self.residue.rjust(6)
        output = output + "  " + ("%.4f" % self.charge).rjust(7)
        return output


    def add_neighbor_atom_index(self, index):
        """Add to the list of neighbor (bonded) atoms the index (int) of a recently identified new-neighbor atom."""
        
        if not (index in self.indecies_of_atoms_connecting):
            self.indecies_of_atoms_connecting.append(index)

    def atomWeight(self,):
            atom = self.element
            mass = 0.0
            if atom == 'H': mass = 1.00794
            elif atom == 'C': mass = 12.0107
            elif atom == 'CL': mass = 35.453
            elif atom == 'N': mass = 14.0067
            elif atom == 'O': mass = 15.9994
            elif atom == 'P': mass = 30.973762
            elif atom == 'S': mass = 32.065
            elif atom == 'BR': mass = 79.904
            elif atom == 'I': mass = 126.90447
            elif atom == 'F': mass = 18.9984032
            elif atom == 'B': mass = 24.3051
            elif atom == 'HG': mass = 200.59
            elif atom == 'BI': mass = 208.98040
            elif atom == 'AS': mass = 74.92160
            elif atom == 'AG': mass = 107.8682
            elif atom == 'K': mass = 39.0983
            elif atom == 'LI': mass = 6.941
            elif atom == 'MG': mass = 24.3050
            elif atom == 'RH': mass = 102.90550
            elif atom == 'ZN': mass = 65.38
            return mass


class Molecule:
    """Loads, saves, and manupulates molecuar models."""
    
    def __init__ (self):
        """Initializes the variables of the Molecule class."""
        
        self.all_atoms={}
        self.remarks = []
        self.babel_path = ["/usr/local/bin/babel", "/usr/bin/babel"]
        self.filename = ""
        self.bonds_created_by_distance = False
        
    def copy_of(self): # this makes a copy of the Molecule object.
        """Returns an exact copy (pymolecule.Molecule) of this Molecule object."""
        
        new_molecule = Molecule()
        for key in self.all_atoms.keys():
            new_molecule.all_atoms[key] = self.all_atoms[key].copy_of()
        for remark in self.remarks:
            new_molecule.remarks.append(remark)
        new_molecule.filename = self.filename
        new_molecule.bonds_created_by_distance = self.bonds_created_by_distance 
        
        # note that undo points are not copied
        return new_molecule
    
    def change_residue(self, NewResName):
        """Changes the residue name of all atoms.
        
        Arguments:
        NewResName -- a string, the new residue name
        
        """
        
        for atom_index in self.all_atoms:
            self.all_atoms[atom_index].change_residue(NewResName)

    def load_pdb(self, filename):
        """Loads a PDB file into the current Molecule object. This function is now depreciated. Used load_pdb2 instead.
        
        Arguments:
        filename -- a string, the name of the file to load

        """

        self.__init__()
        
        self.filename = filename
            
        # Now load the file into a list
        file = open(filename,"r")
        lines = file.readlines()
        file.close()
        
        for t in range(0,len(lines)):
            line=lines[t]
            if len(line) >= 7:
                if line[0:7]=="REMARK ":
                    self.remarks.append(line[7:-1])
                elif line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    temp_atom = Atom()
                    temp_atom.read_pdb_line(line)
                    self.all_atoms[int(temp_atom.molecule_index)] = temp_atom
        
        self.create_bonds_by_distance()

    def load_pdb2(self, filename, bonds_by_distance):
        """Loads a PDB file into the current Molecule object
        
        Arguments:
        filename -- a string, the name of the file to load
        bonds_by_distance -- A boolean, whether or not to determine atomic bonds by distance

        """

        self.__init__()
        
        self.filename = filename
            
        # Now load the file into a list
        file = open(filename,"r")
        lines = file.readlines()
        file.close()
        
        for t in range(0,len(lines)):
            line=lines[t]
            if len(line) >= 7:
                if line[0:7]=="REMARK ":
                    self.remarks.append(line[7:-1])
                elif line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    temp_atom = Atom()
                    temp_atom.read_pdb_line(line)
                    self.all_atoms[int(temp_atom.molecule_index)] = temp_atom
        
        if bonds_by_distance == True: self.create_bonds_by_distance()

    def loadmol2(self, filename, bonds_by_distance):
        """Loads a mol2 file into the current Molecule object

        Arguments:
        filename -- a string, the name of the file to load
        bonds_by_distance -- A boolean, whether or not to determine atomic bonds by distance

        """

        self.__init__()

        self.filename = filename

        # Now load the fiel into a list
        file = open(filename,"r")
        lines = file.readlines()
        file.close()

        for n in range(0,len(lines)):
            line = lines[n]
            if len(line) >= 7:
                if "@<TRIPOS>ATOM" in line:
                    i=0
                    for m in range(n+1,len(lines)):
                       atomline =lines[m]
                       i += 1

                       if len(atomline) <= 50:
                           break

                       if i ==int(atomline[0:6]):

                           temp_atom = Atom()
                           temp_atom.read_mol2_line((atomline))
                           self.all_atoms[int(temp_atom.molecule_index)] = temp_atom


                    break
        if bonds_by_distance == True: self.create_bonds_by_distance









    def divide_molecule_into_constituent_parts(self):
        """If this Molecule object contains multiple, isolated compounds (for example, a receptor and ligand that are not covalently bonded), separate them into distinct Molecule objects
        
        Returns:
        A list of Molecule objects, each one corresponding to an isolated molecule.
        
        """
        
        indices = self.all_atoms.keys()
        isolated = [] #the list that will contain multiple Molecule objects
        
        while len(indices) > 0:
        
            # molecule contains the indices of each isolated molecule
            molecule = [indices[0]]
            t = 0
            while t < len(molecule):
                current_index = molecule[t]
                atom = self.all_atoms[current_index]
                for neighbor in atom.indecies_of_atoms_connecting:
                    if not neighbor in molecule: molecule.append(neighbor)
                t = t + 1
            
            # now make the isolated pdb and saveit
            new_pdb = Molecule()
            for index in molecule: new_pdb.all_atoms[index] = self.all_atoms[index].copy_of()
            new_pdb.create_bonds_by_distance()
            isolated.append(new_pdb.copy_of())
            
            # now remove the indices from indices list so you can go on to the next isolated molecule
            for index in molecule: indices.remove(index)
        
        return isolated
    
    def steric_clashes(self):
        """Detects steric clashes.
        
        Returns:
        A boolean. True if steric clashes are present, False if they are not.
        
        """
        
        steric_clashes_present = False
        
        keys1 = self.all_atoms.keys()
        keys2 = self.all_atoms.keys()
        for index1 in range(len(keys1)-1):
            if steric_clashes_present == False:
                atomindex1 = keys1[index1]
                atom1 = self.all_atoms[atomindex1]
                
                if self.bonds_created_by_distance == True: # so only check number of bonds if you've determined those bonds previously
                    #print atom1.element, len(atom1.indecies_of_atoms_connecting)
                    steric_clashes_present = self.__atom_has_too_many_connections(atomindex1)
                    if steric_clashes_present == True: break
                    
                for index2 in range(index1+1, len(keys2)):
                    atomindex2 = keys2[index2]
                    atom2 = self.all_atoms[atomindex2]
                    
                    ideal_length = self.bond_length(atom1.element, atom2.element)
                    
                    dist = atom1.coordinates.distance_to_another_point(atom2.coordinates)
                    if dist < ideal_length * 0.5: # this is too short to be a bond with a hydrogen or anything.
                        steric_clashes_present = True
                        break # breaks out of first loop
        
                    # a second check
                    if (atom1.element != "H" and atom2.element != "H") or (atom1.element == "H" and atom2.element == "H"):
                        cut_off_dist = 1.0
                    else: # so it's a bond between hydrogen and some other atom
                        cut_off_dist = 0.5
                    if dist < cut_off_dist: # this is too short to be a bond with a hydrogen or anything.
                        steric_clashes_present = True
                        break # breaks out of first loop

        if steric_clashes_present == False: # check last atom too
            steric_clashes_present = self.__atom_has_too_many_connections(keys1[-1])

        '''# now an extra hydrogen check. If a hydrogen comes withn 1.6 of more than one atom, there's steric problems! # THIS IS NO GOOD!!! DON'T USE!!
        if steric_clashes_present == False:
            for index1 in self.all_atoms.keys():
                atom1 = self.all_atoms[index1]
                if atom1.element == "H":
                    num_within_radius = 0
                    for index2 in self.all_atoms.keys():
                        if index1 != index2:
                            atom2 = self.all_atoms[index2]
                            dist = atom1.coordinates.distance_to_another_point(atom2.coordinates)
                            if dist < 1.3:
                                num_within_radius = num_within_radius + 1
                                if num_within_radius > 1: break # no reason to keep counting... there's two or more other atoms within the radius
                    if num_within_radius > 1:
                        steric_clashes_present = True
                        break'''

        return steric_clashes_present
    
    def __atom_has_too_many_connections(self, atom_index):
        
        atom = self.all_atoms[atom_index]
        
        #response = False
        if atom.element == "C" and len(atom.indecies_of_atoms_connecting) > 4:
            #print "C more than 4"
            #response = True
            return True
        if atom.element == "N" and len(atom.indecies_of_atoms_connecting) > 4:
            #print "N more than 4"
            #response = True
            return True
        if atom.element == "O" and len(atom.indecies_of_atoms_connecting) > 2:
            #print "O more than 2"
            #response = True
            return True
        if atom.element == "H" and len(atom.indecies_of_atoms_connecting) > 1:
            #print "H more than 1"
            #response = True
            return True
        if atom.element == "F" and len(atom.indecies_of_atoms_connecting) > 1:
            #print "F more than 1"
            #response = True
            return True
        if atom.element == "Br" and len(atom.indecies_of_atoms_connecting) > 1:
            #print "Br more than 1"
            #response = True
            return True
        if atom.element == "Cl" and len(atom.indecies_of_atoms_connecting) > 1:
            #print "Cl more than 1"
            #response = True
            return True
        if atom.element == "BR" and len(atom.indecies_of_atoms_connecting) > 1:
            #print "BR more than 1"
            #response = True
            return True
        if atom.element == "CL" and len(atom.indecies_of_atoms_connecting) > 1:
            #print "CL more than 1"
            #response = True
            return True
        if atom.element == "C" and len(atom.indecies_of_atoms_connecting) <= 4: # additional geometry check
            response = self.__are_angles_with_neighbors_not_equal(atom_index)
            if response == True: return True
        if atom.element == "N" and len(atom.indecies_of_atoms_connecting) <= 4: # additional geometry check
            response = self.__are_angles_with_neighbors_not_equal(atom_index)
            if response == True: return True
        if atom.element == "P": # additional geometry check
            response = self.__are_angles_with_neighbors_not_equal(atom_index)
            if response == True: return True
        if atom.element == "S": # additional geometry check
            response = self.__are_angles_with_neighbors_not_equal(atom_index)
            if response == True: return True

        return False

    def __are_angles_with_neighbors_not_equal(self, index): # diagnoses additional structural problems
        
        response = False
        angles = []
        
        central_atom = self.all_atoms[index]
        for t1 in range(0,len(central_atom.indecies_of_atoms_connecting)-1):
            index1 = central_atom.indecies_of_atoms_connecting[t1]
            pt1 = self.all_atoms[index1].coordinates
            for t2 in range(t1+1,len(central_atom.indecies_of_atoms_connecting)):
                index2 = central_atom.indecies_of_atoms_connecting[t2]
                if index1 != index2:
                    pt2 = self.all_atoms[index2].coordinates
                    angles.append(pt1.angle_between_three_points(central_atom.coordinates,pt2) * 180 / math.pi)
                    

        #print angles
        
        if len(angles) < 2: return False
        
        for t1 in range(0,len(angles)-1):
            angle1 = angles[t1]
            for t2 in range(t1+1,len(angles)):
                angle2 = angles[t2]
                if math.fabs(angle1 - angle2) > 100.0: # why is this angle so great? Because angles can be legitimately strained in three-member rings, for example. A big tolerance is needed her. We'll miss some, but we'll catch a lot too. I looked at a group of 5000 ligands to come to this conclusion.
                    #print math.fabs(angle1 - angle2)
                    response = True
                    break
                
        return response
    
    def save_pdb(self, filename):
        """Saves data to a PDB file. This function is now depreciated. Use the Molecule.save_pdb2() function instead
        
        Arguments:
        filename -- A string, the filename to be written.
        
        """

        toprint = ""

        if len(self.all_atoms) > 0: # so the pdb is not empty (if it is empty, don't save)
            
            # check to see if there are any steric clashes
            if self.steric_clashes() == False:
                self.change_residue("LIG") # so only one residue name for all atoms in the ligand
            else: # so there are steric clashes
                # note that each fragment of the residue has a different name to aid with debugging
                filename = filename + '.possible_steric_clashes.pdb'
                toprint = toprint + "Lignad with possible steric clashes has been created: " + filename + "\n"
            
            file = open(filename,"w")
    
            # write comments. include derivation information. maybe include reaction information.
            for line in self.remarks:
                    file.write("REMARK "+line+"\n")
            
            # write coordinates
            for atomindex in self.all_atoms:
                file.write(self.all_atoms[atomindex].create_pdb_line(atomindex) + "\n")
    
            # add connect data
            if self.bonds_created_by_distance == True:
                for atomindex1 in self.all_atoms:
                    atom = self.all_atoms[atomindex1]
                    for atomindex2 in atom.indecies_of_atoms_connecting:
                        file.write("CONECT" + str(atomindex1).rjust(5, " ") + str(atomindex2).rjust(5, " ")+"\n")

            file.close()
            
            toprint = toprint + "New PDB created: " + filename
        else: print "Error: Cannot save a Molecule with no atoms (file name \"" + filename + "\")"
            
        return toprint
    
    def save_pdb2(self, filename, check_for_steric_clashes):
        """Saves data to a PDB file.
        
        Arguments:
        filename -- A string, the filename to be written.
        check_for_steric_clases -- A boolean, True if a check for steric clashes is to be performed, False otherwise
        
        """

        toprint = ""

        if len(self.all_atoms) > 0: # so the pdb is not empty (if it is empty, don't save)
            
            if check_for_steric_clashes == True:
                # check to see if there are any steric clashes
                if self.steric_clashes() == True: # so there are steric clashes
                    # note that each fragment of the residue has a different name to aid with debugging
                    filename = filename + '.possible_steric_clashes.pdb'
                    toprint = toprint + "Lignad with possible steric clashes has been created: " + filename + "\n"
            
            file = open(filename,"w")
    
            # write comments. include derivation information. maybe include reaction information.
            for line in self.remarks:
                    file.write("REMARK "+line+"\n")
            
            # write coordinates
            for atomindex in self.all_atoms:
                file.write(self.all_atoms[atomindex].create_pdb_line(atomindex) + "\n")
    
            file.close()
            
            toprint = toprint + "New PDB created: " + filename
        else: print "Error: Cannot save a Molecule with no atoms (file name \"" + filename + "\")"
            
        return toprint

    def save_mol2(self, filename, check_for_steric_clashes):

        """Saves data to a mol2 file.

        Arguments:
        filename -- A string, the filename to be written.
        check_for_steric_clases -- A boolean, True if a check for steric clashes is to be performed, False otherwise

        """
        toprint = ""

        if  len(self.all_atoms) >0 : # so the mol2 is not empty
            if check_for_steric_clashes == True:
                if self.steric_clashes() == True:

                    filename = filename + '.possible_steric_clashes.mol2'
                    toprint = toprint + "Lignad with possible steric clashes has been created: " +filename + "\n"

            file = open(filename, "w")

            for atomindex in self.all_atoms:
                file.write(self.all_atoms[atomindex].create_mol2_line(atomindex) + "\n")

            file.close()



    
    def __str__(self):
        representation = ""
        for atomindex in self.all_atoms: representation = representation + self.all_atoms[atomindex].create_pdb_line(atomindex) + "\n"
        return representation       

    def create_bonds_by_distance(self):
        """Determines which atoms are bound to each other based on their proximity."""
        
        # first remove old bond information
        for AtomIndex in self.all_atoms: self.all_atoms[AtomIndex].indecies_of_atoms_connecting = []
        
        # do a pairwise comparison
        all_atoms_keys = self.all_atoms.keys()
        for t in range(len(all_atoms_keys)-1):
            AtomIndex1 = all_atoms_keys[t]
            
            for s in range(t+1,len(all_atoms_keys)):
                AtomIndex2 = all_atoms_keys[s]

                if AtomIndex1 != AtomIndex2:
                    atom1 = self.all_atoms[AtomIndex1]
                    atom2 = self.all_atoms[AtomIndex2]
                    dist = atom1.coordinates.distance_to_another_point(atom2.coordinates)
                    
                    if (dist < self.bond_length(atom1.element, atom2.element) * 1.2):
                        atom1.add_neighbor_atom_index(AtomIndex2)
                        atom2.add_neighbor_atom_index(AtomIndex1)
                        #print str(AtomIndex1)+"\t"+str(AtomIndex2)
        self.bonds_created_by_distance = True
                        
    def set_undo_point(self): # you can restore all atom positions to some undo point. This sets that point.
        """Sets ("saves") the undo point of all atoms. Any subsequent manipulations of atomic coordinates can be "undone" by reseting to this configuration."""
        
        for atomindex in self.all_atoms:
            self.all_atoms[atomindex].set_undo_point()
            
    def undo(self):
        """Resets the coordinates of all atoms to those saved using the set_undo_point function."""
        
        for atomindex in self.all_atoms:
            self.all_atoms[atomindex].undo()
        
    def set_atom_location(self, AtomIndex, new_location): # translate the protein so that a specific atom has a specific location
        """Translates the entire molecular model so that an atom of specified index is located at a specified coordinate.
        
        Arguments:
        AtomIndex -- An int, the index of the target atom.
        new_location -- A pymolecule.Point object specifying the new location.
        
        Returns: A pymolecule.Point, the vector by which the Molecule was translated
        
        """
        
        currentloc = self.all_atoms[AtomIndex].coordinates
        delta = new_location.subtract(currentloc)
        self.translate_molecule(delta)
        
        return delta
        
    def translate_molecule(self, delta):
        """Translate all atoms of the molecular model by a specified vector.
        
        Arguments:
        delta -- A pymolecule.Point object specifying the amount to move each atom along the x, y, and z coordinates.
        
        """
        
        for atomindex in self.all_atoms:
            self.all_atoms[atomindex].coordinates = self.all_atoms[atomindex].coordinates.add(delta)

    def rotate_molecule_around_pivot(self, PivotIndex, thetax, thetay, thetaz):
        """Rotate the molecular model around a specified atom.
        
        Arguments:
        PivotIndex -- An int, the index of the atom about which the molecular model will be rotated
        thetax -- A float, the angle to rotate relative to the x axis, in radians
        thetay -- A float, the angle to rotate relative to the y axis, in radians
        thetaz -- A float, the angle to rotate relative to the z axis, in radians
        
        """
        
        # First, move the Molecule so the pivot is at the origin
        old_val = self.all_atoms[PivotIndex].coordinates.copy_of()
        self.set_atom_location(PivotIndex,Point(0,0,0))
        
        # now rotate
        for atomindex in self.all_atoms:
            vector = self.all_atoms[atomindex].coordinates
            
            sinx = math.sin(thetax)
            siny = math.sin(thetay)
            sinz = math.sin(thetaz)
            cosx = math.cos(thetax)
            cosy = math.cos(thetay)
            cosz = math.cos(thetaz)
            
            new_x = vector.x * cosy * cosz + vector.y * (sinx * siny * cosz + cosx * sinz) + vector.z * (sinx * sinz - cosx * siny * cosz)
            new_y = -vector.x * cosy * sinz + vector.y * (cosx * cosz - sinx * siny * sinz) + vector.z * (cosx * siny * sinz + sinx * cosz)
            new_z = vector.x * siny - vector.y * sinx * cosy + vector.z * cosx * cosy
    
            self.all_atoms[atomindex].coordinates = Point(new_x, new_y, new_z)
        
        # now move the pivot point back to it's old location
        self.set_atom_location(PivotIndex,old_val)

    def print_out_info_pdb(self):
        """Print out a quick representation of the molecular model, in PDB format."""
        
        print self.remarks
        for index in self.all_atoms:
            print self.all_atoms[index].create_pdb_line(index)

    def print_out_info_mol2(self):
        """Print out a quick representation of the molecular model, in mol2 format."""
        for index in self.all_atoms:
            print self.all_atoms[index].create_mol2_line(index)
    
    def rotate_molecule_around_a_line_use_atom_indicies(self, LinePoint1_Index, LinePoint2_Index, Rotate): # takes as input atom indicies
        """Rotate the molecular model about a line segment. The end points of the line segment are atoms of specified coordinates.
        
        Arguments:
        LinePoint1_Index -- An int, the index of the first atom at one end of the line segment.
        LinePoint2_Index -- An int, the index of the second atom at the other end of the line segment.
        Rotate -- A float, the angle of rotation, in radians.
        
        """
        
        pt1 = self.all_atoms[LinePoint1_Index].coordinates
        pt2 = self.all_atoms[LinePoint2_Index].coordinates
        self.rotate_molecule_around_a_line(pt1, pt2, Rotate)
    
    def rotate_molecule_around_a_line(self, LinePoint1, LinePoint2, Rotate): #Rotate is in radians
        """Rotate the molecular model about a line segment. The end points of the line segment are explicitly specified coordinates.
        
        Arguments:
        LinePoint1 -- A pymolecule.Point corresponding to one end of the line segment.
        LinePoint2 -- A pymolecule.Point corresponding to the other end of the line segment.
        Rotate -- A float, the angle of rotation, in radians.
        
        """

        a = LinePoint1.x
        b = LinePoint1.y
        c = LinePoint1.z
        d = LinePoint2.x
        e = LinePoint2.y
        f = LinePoint2.z
        u = d-a
        v = e-b
        w = f-c
                
        # Now rotate molecule
        for t in self.all_atoms: # so t is an atom index
            OrigAtom = self.all_atoms[t]
            x_not = OrigAtom.coordinates.x
            y_not = OrigAtom.coordinates.y
            z_not = OrigAtom.coordinates.z
            v_2_plus_w_2 = math.pow(v, 2) + math.pow(w, 2)
            u_2_plus_w_2 = math.pow(u, 2) + math.pow(w, 2)
            u_2_plus_v_2 = math.pow(u, 2) + math.pow(v, 2)
            u_2_plus_v_2_plus_w_2 = u_2_plus_v_2 + math.pow(w, 2)
            ux_plus_vy_plus_wz = u*x_not + v*y_not + w*z_not
            cos = math.cos(Rotate)
            sin = math.sin(Rotate)
            
            self.all_atoms[t].coordinates.x = (a*v_2_plus_w_2 + u*(-b*v-c*w+ux_plus_vy_plus_wz)+(-a*v_2_plus_w_2+u*(b*v+c*w-v*y_not-w*z_not)+v_2_plus_w_2*x_not)*cos+math.sqrt(u_2_plus_v_2_plus_w_2)*(-c*v+b*w-w*y_not+v*z_not)*sin)/u_2_plus_v_2_plus_w_2
            self.all_atoms[t].coordinates.y=(b*u_2_plus_w_2 + v*(-a*u-c*w+ux_plus_vy_plus_wz)+(-b*u_2_plus_w_2+v*(a*u+c*w-u*x_not-w*z_not)+u_2_plus_w_2*y_not)*cos+math.sqrt(u_2_plus_v_2_plus_w_2)*(c*u-a*w+w*x_not-u*z_not)*sin)/u_2_plus_v_2_plus_w_2
            self.all_atoms[t].coordinates.z=(c*u_2_plus_v_2 + w*(-a*u-b*v+ux_plus_vy_plus_wz)+(-c*u_2_plus_v_2+w*(a*u+b*v-u*x_not-v*y_not)+u_2_plus_v_2*z_not)*cos+math.sqrt(u_2_plus_v_2_plus_w_2)*(-b*u+a*v-v*x_not+u*y_not)*sin)/u_2_plus_v_2_plus_w_2

    def regenerate_bonds(self):
        """Recalculate all bonds based on proximity."""
        
        # first delete all the existing bond information
        for atomindex in self.all_atoms:
            del self.all_atoms[atomindex].indecies_of_atoms_connecting[:]
        
        # now recreate the bonds.
        self.create_bonds_by_distance()
    
    def in_same_ring(self, index1, index2): # indicies must be adjacent
        """Determine if two atoms are members of the same ring system.
        
        Arguments:
        index1 -- An int, the index of the first atom to be considered.
        index2 -- An int, the index of the second atom to be considered.
        
        Returns:
        A string, "TRUE" if the specified atoms are in the same ring system, "FALSE" otherwise.
        
        """
        
        # This is going to be a recursive algorithm that walks through the molecular model.
        AlreadyCrossed = [index2]
        
        for neighbor_index in self.all_atoms[index2].indecies_of_atoms_connecting:
            if neighbor_index != index1:
                self.__ring_recursive_walk(neighbor_index, AlreadyCrossed, index1)
    
        if index1 in AlreadyCrossed: response = 'TRUE'
        else: response = 'FALSE'
    
        return response

    def __ring_recursive_walk(self, index, AlreadyCrossed, starting_point):
        """A recursive function used by in_same_ring() to determine whether two atoms are in the same ring system."""
        
        # first, get a list of all the atoms taht atom:index is connected to.
        connecteds = self.all_atoms[index].indecies_of_atoms_connecting[:] # this creates a copy of the list
        
        # now delete any indicies from the list that you've already interated over
        for already in AlreadyCrossed:
            if already in connecteds:
                connecteds.remove(already)
        
        # The current index has now been interated over
        if not index in AlreadyCrossed: AlreadyCrossed.append(index)
        
        #Now go on to the neighbors
        for idx in connecteds:
            self.__ring_recursive_walk(idx, AlreadyCrossed, starting_point)

    
    def number_of_neighors_of_element(self, index, the_element):
        """Counts the number of atoms of a given element bonded to a specified atom of interest.
        
        Arguments:
        index -- An int, the index of the atom of interest.
        the_element -- A string specifying the element of the neighbors.
        
        Returns:
        An int, the number of neighboring atoms matching the criteria.
        
        """
        
        num = 0
        for index in self.all_atoms[index].indecies_of_atoms_connecting:
            if self.all_atoms[index].element == the_element: num = num + 1
        return num
    
    def index_of_neighbor_of_element(self,index,the_element): # returns the index of the first neighbor of specified element
        """For a given atom of interest, returns the index of the first neighbor of a specified element.
        
        Arguments:
        index -- An int, the index of the atom of interest.
        the_element -- A string specifying the desired element of the neighbor
        
        Returns:
        An int, the index of the first neighbor that is the specified element. If no such neighbor exists, returns -1.
        
        """
        
        for index in self.all_atoms[index].indecies_of_atoms_connecting:
            if self.all_atoms[index].element == the_element: return index
        return -1 # returns -1 if no match
    
    def number_of_atoms(self):
        """Returns the number of atoms (int) in the molecular model."""
        
        return len(self.all_atoms)
    
    def weight(self):
        """Returns the weight (float) of the molecular model."""
        
        total = 0.0
        for atomindex in self.all_atoms:
            atom = self.all_atoms[atomindex].element
            mass = 0.0
            if atom == 'H': mass = 1.00794
            elif atom == 'C': mass = 12.0107
            elif atom == 'CL': mass = 35.453
            elif atom == 'N': mass = 14.0067
            elif atom == 'O': mass = 15.9994
            elif atom == 'P': mass = 30.973762
            elif atom == 'S': mass = 32.065
            elif atom == 'BR': mass = 79.904
            elif atom == 'I': mass = 126.90447
            elif atom == 'F': mass = 18.9984032
            elif atom == 'B': mass = 24.3051
            elif atom == 'HG': mass = 200.59
            elif atom == 'BI': mass = 208.98040
            elif atom == 'AS': mass = 74.92160
            elif atom == 'AG': mass = 107.8682
            elif atom == 'K': mass = 39.0983
            elif atom == 'LI': mass = 6.941
            elif atom == 'MG': mass = 24.3050
            elif atom == 'RH': mass = 102.90550
            elif atom == 'ZN': mass = 65.38
            total += mass
        return total
    
    def smiles_string(self):
        """If the full path to the Open Babel executable has been specified, returns the SMILES string (string) of the molecular model."""
        
        for filename in self.babel_path:
            if os.path.isfile(filename):
                self.save_pdb('temp_pdb_from_pymolecule.pdb')
                SMILES = self.__run_file(filename + " -ipdb temp_pdb_from_pymolecule.pdb -ocan")
                SMILES = SMILES.split("\t")
                SMILES = SMILES[0]
                os.remove('temp_pdb_from_pymolecule.pdb')
                return SMILES

    def __run_file(self, command): # a private function
        """Executes a specified file (string) and returns the output (string)."""
        
        p = os.popen(command)
        s = p.readlines()
        p.close()
        total=""
        for item in s:
            total = total + item
        return total
    
    # Determining hybridization
    def hybridization(self, index):
        """Guesses at the hybridization of a specified atom.
        
        Arguments:
        index -- An int, the index of the atom whose hybridization is to be determined.
        
        Returns:
        An int, where 3 corresponds to sp3 hybridization, 2 corresponds to sp2 hybridization, and 1 corresponds to sp1 hybridization.
        
        """
        
        the_atom = Atom()
        the_atom = self.all_atoms[index]
        element = the_atom.element
        num_neighbors = the_atom.number_of_neighbors()
        ans = -1
        
        if element == 'S': 
            if num_neighbors == 2 or num_neighbors == 4: # don't know if this is right. I'm assuming an S in an SO3 group is SP3 hybridized
                ans = 3
            else: 
                ans = 2
        if element == 'C':
            if num_neighbors == 4:
                ans = 3
            elif num_neighbors == 3:
                ans = 2
            elif num_neighbors == 2:
                ans = 1
        elif element == 'O':
            if num_neighbors == 2:
                ans = 3
            elif num_neighbors == 1:
                ans = 2
        elif element == "N":
            if num_neighbors == 4:
                ans = 3
            elif num_neighbors == 3: # this is the tricky one.
                point1 = the_atom.coordinates
                point2 = self.all_atoms[the_atom.indecies_of_atoms_connecting[0]].coordinates
                point3 = self.all_atoms[the_atom.indecies_of_atoms_connecting[1]].coordinates
                point4 = self.all_atoms[the_atom.indecies_of_atoms_connecting[2]].coordinates
                plane_dist = point1.planarity(point2, point3, point4) # This is 0 for an ideal sp2 hybridized nitrogen, and 0.38 for an ideal sp3 hybridized nitrogen
                if plane_dist > 0.2:
                    ans = 3
                else:
                    ans = 2
            elif num_neighbors == 2: # If there are two, it could be sp2 or sp3 (what if it's in a benzene ring, for example?). Determine based on neighbors
                
                ans = 3 # by default, assume SP3 hybridized.
                
                index2 = the_atom.indecies_of_atoms_connecting[0]
                index3 = the_atom.indecies_of_atoms_connecting[1]
                
                if self.all_atoms[index2].element == 'C':
                    if self.hybridization(index2) == 2:
                        ans = 2
                
                if self.all_atoms[index3].element == 'C':
                    if self.hybridization(index3) == 2:
                        ans = 2
                
            else: #So connected to only 1
                ans = 1
                
        return ans
    
    def bond_length(self, element1, element2):
        """Returns the ideal length of a bond between two atoms of specified elements.
        
        Arguments:
        element1 -- A string corresponding to the element of the first atom.
        element2 -- A string corresponding to the element of the second atom.
        
        Returns:
        A float, corresponding to the ideal bond length.
        
        """
        
        '''Bond lengths taken from Handbook of Chemistry and Physics. The information provided there was very specific,
        so I tried to pick representative examples and used the bond lengths from those. Sitautions could arise where these
        lengths would be incorrect, probably slight errors (<0.06) in the hundreds.'''

	# occasionally element names might have + and - in them. Ignore anything after those.
	if "+" in element1: element1 = element1.split("+")[0]
	if "-" in element1: element1 = element1.split("-")[0]
	if "+" in element2: element2 = element2.split("+")[0]
	if "-" in element2: element2 = element2.split("-")[0]
       
        distance = 0.0
        if element1 == "C" and element2 == "C": distance = 1.53
        if element1 == "N" and element2 == "N": distance = 1.425
        if element1 == "O" and element2 == "O": distance = 1.469
        if element1 == "S" and element2 == "S": distance = 2.048
        if (element1 == "C" and element2 == "H") or (element1 == "H" and element2 == "C"): distance = 1.059
        if (element1 == "C" and element2 == "N") or (element1 == "N" and element2 == "C"): distance = 1.469
        if (element1 == "C" and element2 == "O") or (element1 == "O" and element2 == "C"): distance = 1.413
        if (element1 == "C" and element2 == "S") or (element1 == "S" and element2 == "C"): distance = 1.819
        if (element1 == "N" and element2 == "H") or (element1 == "H" and element2 == "N"): distance = 1.009
        if (element1 == "N" and element2 == "O") or (element1 == "O" and element2 == "N"): distance = 1.463
        if (element1 == "O" and element2 == "S") or (element1 == "S" and element2 == "O"): distance = 1.577
        if (element1 == "O" and element2 == "H") or (element1 == "H" and element2 == "O"): distance = 0.967
        if (element1 == "S" and element2 == "H") or (element1 == "H" and element2 == "S"): distance = 2.025/1.5 # This one not from source sited above. Not sure where it's from, but it wouldn't ever be used in the current context ("AutoGrow")
        if (element1 == "S" and element2 == "N") or (element1 == "H" and element2 == "N"): distance = 1.633
    
        if (element1 == "C" and element2 == "F") or (element1 == "F" and element2 == "C"): distance = 1.399
        if (element1 == "C" and element2 == "CL") or (element1 == "CL" and element2 == "C"): distance = 1.790
        if (element1 == "C" and element2 == "BR") or (element1 == "BR" and element2 == "C"): distance = 1.910
        if (element1 == "C" and element2 == "I") or (element1 == "I" and element2 == "C"): distance=2.162
    
        if (element1 == "S" and element2 == "BR") or (element1 == "BR" and element2 == "S"): distance = 2.321
        if (element1 == "S" and element2 == "CL") or (element1 == "CL" and element2 == "S"): distance = 2.283
        if (element1 == "S" and element2 == "F") or (element1 == "F" and element2 == "S"): distance = 1.640
        if (element1 == "S" and element2 == "I") or (element1 == "I" and element2 == "S"): distance= 2.687
    
        if (element1 == "P" and element2 == "BR") or (element1 == "BR" and element2 == "P"): distance = 2.366
        if (element1 == "P" and element2 == "CL") or (element1 == "CL" and element2 == "P"): distance = 2.008
        if (element1 == "P" and element2 == "F") or (element1 == "F" and element2 == "P"): distance = 1.495
        if (element1 == "P" and element2 == "I") or (element1 == "I" and element2 == "P"): distance = 2.490
        if (element1 == "P" and element2 == "C") or (element1 == "C" and element2 == "P"): distance = 1.841
        if (element1 == "P" and element2 == "N") or (element1 == "N" and element2 == "P"): distance = 1.730
        if (element1 == "P" and element2 == "O") or (element1 == "O" and element2 == "P"): distance = 1.662
        if (element1 == "P" and element2 == "S") or (element1 == "S" and element2 == "P"): distance = 1.954

        if (element1 == "N" and element2 == "BR") or (element1 == "BR" and element2 == "N"): distance = 1.843
        if (element1 == "N" and element2 == "CL") or (element1 == "CL" and element2 == "N"): distance = 1.743
        if (element1 == "N" and element2 == "F") or (element1 == "F" and element2 == "N"): distance = 1.406
        if (element1 == "N" and element2 == "I") or (element1 == "I" and element2 == "N"): distance= 2.2
    
        if (element1 == "SI" and element2 == "BR") or (element1 == "BR" and element2 == "SI"): distance = 2.284
        if (element1 == "SI" and element2 == "CL") or (element1 == "CL" and element2 == "SI"): distance = 2.072
        if (element1 == "SI" and element2 == "F") or (element1 == "F" and element2 == "SI"): distance = 1.636
        if (element1 == "SI" and element2 == "P") or (element1 == "P" and element2 == "SI"): distance= 2.264
        if (element1 == "SI" and element2 == "S") or (element1 == "S" and element2 == "SI"): distance= 2.145
        if (element1 == "SI" and element2 == "SI") or (element1 == "SI" and element2 == "SI"): distance= 2.359
        if (element1 == "SI" and element2 == "C") or (element1 == "C" and element2 == "SI"): distance= 1.888
        if (element1 == "SI" and element2 == "N") or (element1 == "N" and element2 == "SI"): distance= 1.743
        if (element1 == "SI" and element2 == "O") or (element1 == "O" and element2 == "SI"): distance= 1.631
        
        if element1 == "H" or element2 == "H": distance= 1.059
        
        if (element1 == "B" and element2 == "C") or (element1 == "C" and element2 == "B"): distance = 1.65
        if element1 == "B" and element2 == "B": distance = 1.7
        if (element1 == "B" and element2 == "F") or (element1 == "F" and element2 == "B"): distance = 1.37
        if (element1 == "B" and element2 == "N") or (element1 == "N" and element2 == "B"): distance = 1.42
        if element1 == "BR" and element2 == "BR": distance = 2.54
        if (element1 == "BR" and element2 == "O") or (element1 == "O" and element2 == "BR"): distance = 1.58
        if element1 == "CL" and element2 == "CL": distance = 2.26
        if (element1 == "CL" and element2 == "F") or (element1 == "F" and element2 == "CL"): distance = 0.001
        if (element1 == "CL" and element2 == "O") or (element1 == "O" and element2 == "CL"): distance = 1.41
        if (element1 == "C" and element2 == "SE") or (element1 == "SE" and element2 == "C"): distance = 1.9
        if element1 == "F" and element2 == "F": distance = 0.001
        if (element1 == "F" and element2 == "O") or (element1 == "O" and element2 == "F"): distance = 0.001
        if (element1 == "N" and element2 == "SE") or (element1 == "SE" and element2 == "N"): distance = 1.8
        if (element1 == "N" and element2 == "S") or (element1 == "S" and element2 == "N"): distance = 1.6
        if (element1 == "BR" and element2 == "CL") or (element1 == "CL" and element2 == "BR"): distance = 2.38
        if (element1 == "BR" and element2 == "F") or (element1 == "F" and element2 == "BR"): distance = 0.001
        if (element1 == "CL" and element2 == "I") or (element1 == "I" and element2 == "CL"): distance = 2.55
        if (element1 == "F" and element2 == "I") or (element1 == "I" and element2 == "F"): distance = 0.001
        if (element1 == "I" and element2 == "O") or (element1 == "O" and element2 == "I"): distance = 2.0
        if (element1 == "B" and element2 == "O") or (element1 == "O" and element2 == "B"): distance = 1.4
        if (element1 == "BR" and element2 == "I") or (element1 == "I" and element2 == "BR"): distance = 2.675
        if (element1 == "B" and element2 == "S") or (element1 == "S" and element2 == "B"): distance = 1.9
        if element1 == "P" and element2 == "P": distance = 2.1

        # X is a useful placeholder atom type you can use to protect certain groups
        if element1 == "X" and element2 == "X": distance = 1.53
        if (element1 == "X" and element2 == "C") or (element1 == "C" and element2 == "X"): distance = 1.53
        if (element1 == "X" and element2 == "H") or (element1 == "H" and element2 == "X"): distance = 1.059
        if (element1 == "X" and element2 == "N") or (element1 == "N" and element2 == "X"): distance = 1.469
        if (element1 == "X" and element2 == "O") or (element1 == "O" and element2 == "X"): distance = 1.413
        if (element1 == "X" and element2 == "S") or (element1 == "S" and element2 == "X"): distance = 1.819
        if (element1 == "X" and element2 == "F") or (element1 == "F" and element2 == "X"): distance = 1.399
        if (element1 == "X" and element2 == "CL") or (element1 == "CL" and element2 == "X"): distance = 1.790
        if (element1 == "X" and element2 == "BR") or (element1 == "BR" and element2 == "X"): distance = 1.910
        if (element1 == "X" and element2 == "I") or (element1 == "I" and element2 == "X"): distance=2.162
        if (element1 == "SI" and element2 == "X") or (element1 == "X" and element2 == "SI"): distance= 1.888
        
        if distance == 0.0:
            print "Warning: Bond length of " + element1 + "-" + element2 + " unknown, Guessing 2.2 A."
            distance = 2.2 # so the atom pairs couldn't be matched. Pick reasonable distance.
        
        return distance

    def align_another_molecule_to_this_one(self, molecule_to_align, tethers): # tethers conain pairs of self_index and molecule_to_align_index. The distances will be minimized. self stays static, molecule_to_align moves
        """Aligns one molecule file to another.
        
        Arguments:
        molecule_to_align -- A molecular model (pymolecule.Molecule object) to be aligned to this model.
        tethers -- A list of pairs, where the first element of each pair is the index (int) of an atom in self,
                   and the second element is the index (int) of an atom in molecule_to_align. molecule_to_align will be
                   rotated and translated so as to minimize the distance between these "tethers."
        
        Returns:
        A pymolecule.Molecule object corresponding to molecule_to_align, now aligned.
        
        """
    
        # The first tether will be exact. It's the pivot point.
        
        molecule_to_align.set_atom_location(tethers[0][1], self.all_atoms[tethers[0][0]].coordinates)
        molecule_to_align.set_undo_point()
        
        if len(tethers) == 2: # If there's only one tether point, you're already done
            
            # Now rotate the second Molecule around that pivot point so that the distance between all other tethers is minimized
            
            variables = {'thetax':0, 'thetay':0, 'thetaz':0}
            constants = {'molecule1':self, 'molecule2':molecule_to_align, 'tethers':tethers}
            
            self.__find_minimum(self.__align_minimization_function, variables, constants) # This should change molecule_to_align to the aligned function.
        
        elif len(tethers) == 3: # so there's three points
    
            variables = {'thetax':0, 'thetay':0, 'thetaz':0}
            constants = {'molecule1':self, 'molecule2':molecule_to_align, 'tethers':tethers[0:2]} # so last tether is missing. This just aligns first two points.
            
            self.__find_minimum(self.__align_minimization_function, variables, constants) # This should change molecule_to_align to the aligned function.
    
            # now do a second minimization, rotating around the line just minimized
            molecule_to_align.set_undo_point()
            variables = {'theta':0}
            constants = {'molecule1':self, 'molecule2':molecule_to_align, 'tethers':tethers} 
            
            self.__find_minimum(self.__align_minimization_function_2, variables, constants) # This should change molecule_to_align to the aligned function.
            
        else:
            print "Error because more than 3 tethers specified!!!"
            
        return molecule_to_align
    
    def __align_minimization_function(self, variables, constants): #molecule1, molecule2, tethers, thetax, thetay, thetaz):
        """A function that calculates the sum of all distances between tethered atoms [see align_another_molecule_to_this_one()].
        This is the function that is minimized to effectuate the alignment of two molecular models."""
        
        thetax = variables['thetax']
        thetay = variables['thetay']
        thetaz = variables['thetaz']
            
        molecule1 = constants['molecule1']
        molecule2 = constants['molecule2']
        
        tethers = constants['tethers']
        
        molecule2_PivotIndex = tethers[0][1]
        
        molecule2.undo()
        molecule2.rotate_molecule_around_pivot(molecule2_PivotIndex, thetax, thetay, thetaz)
    
        dist = 0
        
        for t in range(1,len(tethers)): # not 0, because that's the pivot point that's already right on top of one of the molecule1 atoms
            pair = tethers[t]
            coor1 = molecule1.all_atoms[pair[0]].coordinates
            coor2 = molecule2.all_atoms[pair[1]].coordinates
            
            dist = dist + coor1.distance_to_another_point(coor2)
            
        return dist
            
    def __align_minimization_function_2(self, variables, constants): 
        """A function that calculates the distance between the last of the tethered-atom pairs [see align_another_molecule_to_this_one()].
        This is one of the functions that is minimized to effectuate the alignment of two molecular models when three thethers are present."""
        
        theta = variables['theta']
    
        molecule1 = constants['molecule1']
        molecule2 = constants['molecule2']
        
        tethers = constants['tethers']
        
        molecule2_RotateLine_Index1 = tethers[0][1]
        molecule2_RotateLine_Index2 = tethers[1][1]
        
        molecule2.undo()
        molecule2.rotate_molecule_around_a_line_use_atom_indicies(molecule2_RotateLine_Index1, molecule2_RotateLine_Index2, theta)
        
        dist = 0
        
        pair = tethers[-1] # just look at last point, since first two already aligned.
        coor1 = molecule1.all_atoms[pair[0]].coordinates
        coor2 = molecule2.all_atoms[pair[1]].coordinates
        
        dist = dist + coor1.distance_to_another_point(coor2)
            
        return dist

    def __find_minimum(self, func_name, variables, constants): # This is the function you use
            """Optimizes variables in order to identify the minimum of a specified function.
            
            Arguments:
            func_name -- The name (string) of the function to be minimized.
            variables -- A list of variables (float) that can be optimized.
            constants -- A list of constants (float) that cannot be optimized.
            
            Returns:
            A list of the variables (float), with values optimized.
            
            """
            
            step = 1
            current_val = 100000000
            last_val = 100000001
            
            for t in range(0,1000):
    
                    if math.fabs(current_val - last_val) < 1e-2 and step < 1e-2: # terminate criteria
                            break
                            
                    last_val = current_val
    
                    changes = 'FALSE'
                    
                    for k in variables.keys():
                            center_x = variables[k]
                            center_val_x = func_name(variables, constants)
                            
                            right_x = center_x + step
                            right_val_x = self.__func_in_one_direction(func_name, variables, constants, k, right_x)
                            
                            left_x = center_x - step
                            left_val_x = self.__func_in_one_direction(func_name, variables, constants, k, left_x)
                            
                            if center_val_x < right_val_x and center_val_x < left_val_x:
                                    variables[k] = center_x
                                    current_val = center_val_x
                            elif right_val_x < center_val_x and right_val_x < left_val_x:
                                    variables[k] = right_x
                                    changes = 'TRUE'
                                    current_val = right_val_x
                            elif left_val_x < center_val_x and left_val_x < right_val_x:
                                    variables[k] = left_x
                                    current_val = left_val_x
                                    changes = 'TRUE'
                                    
                    if changes == 'FALSE':
                            step = step / 2.0
                            
            return variables
            
    def __func_in_one_direction(self, func_name, variables, constants, movable_variable, movable_variable_value): # variables is a dictionary containing all the domain variables
            """A function to support __find_minimum() above."""
            
            variables[movable_variable] = movable_variable_value
            return func_name(variables, constants)
            
    def merge_with_another_molecule(self, other_molecule):
        """Merges two molecular models into a single model.
        
        Arguments:
        other_molecule -- A molecular model (pymolecule.Molecule object).
        
        Returns:
        
        A single pymolecule.Molecule object containing the atoms of this model combined with the atoms of other_molecule.
        
        """
        
        # first, what's the maximum index on self?
        maxindex = 0 
        for atomindex in self.all_atoms:
            if atomindex > maxindex: maxindex = atomindex
        
        # okay, add the atoms of self to new molecule
        new_molecule = Molecule()
        for atomindex in self.all_atoms:
            new_molecule.all_atoms[atomindex]=self.all_atoms[atomindex].copy_of()
            
        # now add the atoms of the second Molecule, updating the indicies
        for atomindex in other_molecule.all_atoms:
            new_molecule.all_atoms[atomindex + maxindex] = other_molecule.all_atoms[atomindex].copy_of()
            
        # now go through the molecule and renumber the indecies
        newer_molecule = Molecule()
        count = 1
        for atomindex in new_molecule.all_atoms:
            newer_molecule.all_atoms[count] = new_molecule.all_atoms[atomindex].copy_of()
            count = count + 1
    
        # now rebuild the bonds
        newer_molecule.regenerate_bonds()
    
        return newer_molecule

    def get_branch(self, index1, index2): # Get the branch starting with atom:index1 and moving in the direction of atom:index2
        """Identify an isolated "branch" of this molecular model.
        
        Arguments:
        index1 -- An int, the index of the first atom in the branch (the "root").
        index2 -- An int, the index of the second atom in the branch, used to establish directionality
        
        Returns:
        A pymolecule.Molecule object containing the atoms of the branch.
        
        """
        
        # This is going to be a recursive algorithm that walks through the molecular model.
        AlreadyCrossed = [index1]
        self.__branch_recursive_walk(index2, AlreadyCrossed)
    
        # Now sort the list of AlreadyCrossed
        AlreadyCrossed.sort()
    
        # Create a new self containing just the branch
        new_molecule = Molecule()
        for index in AlreadyCrossed:
            new_molecule.all_atoms[index] = self.all_atoms[index].copy_of()
            
        # Recreate all the bonds
        new_molecule.regenerate_bonds()
        
        return new_molecule
        
    def __branch_recursive_walk(self, index, AlreadyCrossed): # used with teh GetselfBranch function. This is the recursive part.
        """A recursive function to support pymolecule.Molecule.get_branch()."""
        
        # first, get a list of all the atoms taht atom:index is connected to.
        connecteds = self.all_atoms[index].indecies_of_atoms_connecting[:] # this creates a copy of the list
        
        # now delete any indicies from the list that you've already interated over
        for already in AlreadyCrossed:
            if already in connecteds:
                connecteds.remove(already)
        
        # The current index has now been interated over
        if not index in AlreadyCrossed: AlreadyCrossed.append(index)
        
        #Now go on to the neighbors
        for idx in connecteds:
            self.__branch_recursive_walk(idx, AlreadyCrossed)

    def distance_to_another_molecule(self, other_molecule, self_index_to_ignore, other_molecule_index_to_ignore): # computes the minimum distance between two molecules
        """Returns the minimum distance between any of the atoms of this molecular model and any of the atoms of a second, specified model.
        
        Arguments:
        other_molecule -- A molecular model (pymolecule.Molecule).
        self_index_to_ignore -- An int, the index of an atom from this pymolecule.Molecule that should not be considered when calculating the minimum distance.
        other_molecule_index_to_ignore -- An int, the index of an atom from other_molecule that should not be considered when calculating the minimum distance.
        
        Returns:
        A float, the minimum distance between any two atoms of the two specified molecular models.
        
        """
        
        min_dist = 1e10
        for atomindex1 in self.all_atoms:
            if atomindex1 != self_index_to_ignore:
                for atomindex2 in other_molecule.all_atoms:
                    if atomindex2 != other_molecule_index_to_ignore:
                        atom1 = self.all_atoms[atomindex1]
                        atom2 = other_molecule.all_atoms[atomindex2]
                        dist = atom1.coordinates.distance_to_another_point(atom2.coordinates)
                        if dist < min_dist: min_dist = dist
        return min_dist
    
    def bounding_box(self, padding):
        """Specifies a box that bounds the molecules
        
        Arguments:
        padding -- The bounding box will extend this many angstroms beyond the molecule
        
        Returns:
        A six-tuple of the form (min_x, max_x, min_y, max_y, min_z, max_z) that bounds the molecule.
        
        """
        
        # get the coordinates of the first atom to initalize the min and max variables
        atom = Atom()
        initialized = False
        for index in self.all_atoms.keys():
            atom = self.all_atoms[index]
            
            if initialized == False:
                initialized = True
                min_x = atom.coordinates.x
                max_x = atom.coordinates.x
                min_y = atom.coordinates.y
                max_y = atom.coordinates.y
                min_z = atom.coordinates.z
                max_z = atom.coordinates.z
            else:
                if min_x > atom.coordinates.x: min_x = atom.coordinates.x
                if max_x < atom.coordinates.x: max_x = atom.coordinates.x
                if min_y > atom.coordinates.y: min_y = atom.coordinates.y
                if max_y < atom.coordinates.y: max_y = atom.coordinates.y
                if min_z > atom.coordinates.z: min_z = atom.coordinates.z
                if max_z < atom.coordinates.z: max_z = atom.coordinates.z
    
        min_x = min_x - padding
        min_y = min_y - padding
        min_z = min_z - padding

        max_x = max_x + padding
        max_y = max_y + padding
        max_z = max_z + padding

        return (min_x, max_x, min_y, max_y, min_z, max_z)
    
    def two_atoms_in_same_residue(self,index1, index2):
        """Determines whether two atoms are in the same residue
        
        Arguments:
        index1 -- The index of the first atom
        index2 -- The index of the second atom
        
        Returns:
        A boolean, True if they are in the same residue, False if they're not
        
        """
        
        # two residues are the same if their chains, resnames, and resid's are the same.
        atom1 = self.all_atoms[index1]
        atom2 = self.all_atoms[index2]
        
        if atom1.chain == atom2.chain and atom1.residue == atom2.residue and atom1.resid == atom2.resid: return True
        return False
        
    def remove_rotamers(self):
        """Removes rotamers. Rotamers should be specified by placing a letter before the residue name, should be comprised of 4 letters, and should appear next to each other in the sequence. For example, ATRP vs. BTRP. Keeps the first rotamer."""
        
        # note that this function removes rotamers, but if there are multiple distinct residues at a given location, it will not eliminate the distinct residues
        
        # let's start by grouping all residues
        residues = {}
        keys = self.all_atoms.keys()
        for t in range(len(keys)):
            index = keys[t]
            atom = self.all_atoms[index]
            
            key = atom.chain + "_" + atom.residue[-3:] + "_" + str(atom.resid) # so only counting the last three letters of residue name. The first letter might be an rotamer identifier
            if not key in residues: residues[key] = [[], []] # first list contains unique residue names, second list contains all the indices of the residue atoms
            if not atom.residue in residues[key][0]: residues[key][0].append(atom.residue)
            residues[key][1].append(index)
            
        # Now go through each residue and see if any of them are associated with multiple residues
        residues_with_rotamers = []
        for key in residues.keys():
            val = residues[key]
            if len(val[0]) > 1: residues_with_rotamers.append(val)
        
        # Now identify which atoms need to be deleted, and which need to be renamed
        rename_delete_candidates = []
        for residue_with_rotamers in residues_with_rotamers:
            
            # find the first rotamer name. this is the one to keep
            first_rotamer_name = ""
            rotamer_names = residue_with_rotamers[0]
            for resname in rotamer_names:
                if len(resname.strip()) == 4: # so the first four-letter resname, the first rotamer
                    first_rotamer_name = resname
                    break
            
            # now find the names of the ones to delete
            names_to_delete = []
            for resname in rotamer_names:
                if len(resname.strip()) == 4 and resname != first_rotamer_name: names_to_delete.append(resname)
                
            indices_of_residue = residue_with_rotamers[1]
            rename_delete_candidates.append((first_rotamer_name, names_to_delete,indices_of_residue))

        # now look through the candidates
        delete_atoms_with_these_indices = []
        
        for rename_delete_candidate in rename_delete_candidates:

            indices = rename_delete_candidate[2]

            rename_name = rename_delete_candidate[0]
            # rename residues as appropriate
            for index in indices:
                atom = self.all_atoms[index]
                if atom.residue == rename_name: self.all_atoms[index].change_residue(" " + rename_name[-3:])
            
            delete_names = rename_delete_candidate[1]
            # identify the ones to delete
            for index in indices:
                for delete_name in delete_names:
                    for index in indices:
                        atom = self.all_atoms[index]
                        if atom.residue == delete_name:
                            if not index in delete_atoms_with_these_indices: delete_atoms_with_these_indices.append(index)
                        
        # now delete the atoms that have been identified as rotamers
        for index in delete_atoms_with_these_indices:
            self.delete_atom(index)
            
    def selection(self, chain="", residue="", resid=-9999, name="", element="", index= -9999):
        """Select multiple atoms from the molecule. Set any of the arguments to "" or -9999 to ignore. Logical AND used to combine.
        
        Arguments:
        chain -- A string, the letter name of the chain to be deleted
        residue -- A string, the name of the residue to be deleted (equivalent to resname)
        resid -- An integer, the residue id of the residue to be deleted
        name -- A string, the name of the atom to be deleted
        element -- A string, the name of the element to be deleted
        index -- An integer, the atom index of the atom to be deleted
        
        Returns:
        A list of integers containing the indices of the selection.

        """

        chain = chain.strip()
        residue = residue.strip()
        name = name.strip()
        element = element.strip()
        
        keys_to_select = []
        
        for aindex in self.all_atoms.keys():
            atom = self.all_atoms[aindex]

            if chain == atom.chain.strip() or chain == "":
                if residue == atom.residue.strip() or residue == "":
                    if resid == atom.resid or resid == -9999:
                        if name == atom.atomname.strip() or name == "":
                            if element == atom.element.strip() or element == "":
                                if index == aindex or index == -9999:
                                    keys_to_select.append(aindex)
                                    
        return keys_to_select

    def invert_selection(self, selection):
        """Return a new selection that contains all atoms that are not in a given selection
        
        Arguments:
        selection -- A list of integers containing the indices of the atoms
        
        Returns:
        A list of integers containing the indices of those atoms that are not in selection.

        """
        
        inverted_selection = []
        for index in self.all_atoms.keys():
            if index not in selection: inverted_selection.append(index)
        
        return inverted_selection

    def delete_atom(self, index):
        """Remove an atom from the molecular model.
        
        Arguments:
        index -- An int, the index of the atom to remove.
        
        """
        
        # first, delete the index from the bondedatoms lists
        for atomindex in self.all_atoms:
            if index in self.all_atoms[atomindex].indecies_of_atoms_connecting: self.all_atoms[atomindex].indecies_of_atoms_connecting.remove(index)

        # now, delete the atom
        self.all_atoms.pop(index)
        
    def delete_atoms(self, selection):
        """Remove multiple atoms from the molecule. Set any of the arguments to "" or -9999 to ignore. Logical AND used to combine.
        
        Arguments:
        selection -- A list of integers containing the indices of the atoms to be deleted
        
        """
        
        for index in selection:
            self.delete_atom(index)
            
    def rename_residues(self, selection, new_resname):
        """Rename the resnames of all the atoms in a selection.
        
        Arguments:
        selection -- A list of integers containing the indices of the atoms to be renamed
        new_resname -- A string, the new residue name.
        
        """
        
        # residue names should be 4 letters long
        while len(new_resname) < 4: new_resname = " " + new_resname
        
        for index in selection:
            self.all_atoms[index].residue = new_resname
            self.all_atoms[index].resname = new_resname
            
if __name__ == '__main__':
    print "\npymolecule " + version + "\n"



