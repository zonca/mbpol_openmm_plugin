import mbpolplugin
from simtk.openmm.app import forcefield
from simtk import unit

## @private
class MBPolOneBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolOneBodyForceGenerator()
        forceField._forces.append(generator)

        # <MBPolOneBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolOneBodyForce>
        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for MBPolOneBodyForce_template in element.findall('MBPolOneBodyForce'):
            types = forceField._findAtomTypes(MBPolOneBodyForce_template.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolOneBodyForce_template.attrib['class1'],
                                    MBPolOneBodyForce_template.attrib['class2'],
                                    MBPolOneBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolOneBodyForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolOneBodyForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolOneBodyForce and match their name
                force.addOneBody(atom2.index, atom1.index, atom3.index)

forcefield.parsers["MBPolOneBodyForce"] = MBPolOneBodyForceGenerator.parseElement

## @private
class MBPolTwoBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolTwoBodyForceGenerator()
        forceField._forces.append(generator)

        # <MBPolTwoBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolTwoBodyForce>
        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for MBPolTwoBodyForce_template in element.findall('MBPolTwoBodyForce'):
            types = forceField._findAtomTypes(MBPolTwoBodyForce_template, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolTwoBodyForce_template.attrib['class1'],
                                    MBPolTwoBodyForce_template.attrib['class2'],
                                    MBPolTwoBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolTwoBodyForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolTwoBodyForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolTwoBodyForce and match their name
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

                force.addParticle(v)

forcefield.parsers["MBPolTwoBodyForce"] = MBPolTwoBodyForceGenerator.parseElement

## @private
class MBPolThreeBodyForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolThreeBodyForceGenerator()
        forceField._forces.append(generator)

        # <MBPolThreeBodyForce>
        #     <Residue class1="OW" class2="HW" class3="HW" />
        # </MBPolThreeBodyForce>
        # <AmoebaStretchBendForce stretchBendUnit="1.0">
        # <StretchBend class1="2" class2="1" class3="3" k1="5.25776946506" k2="5.25776946506" />
        # <StretchBend class1="2" class2="1" class3="4" k1="3.14005676385" k2="3.14005676385" />

        for MBPolThreeBodyForce_template in element.findall('MBPolThreeBodyForce'):
            types = forceField._findAtomTypes(MBPolThreeBodyForce_template.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolThreeBodyForce_template.attrib['class1'],
                                    MBPolThreeBodyForce_template.attrib['class2'],
                                    MBPolThreeBodyForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolThreeBodyForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolThreeBodyForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolThreeBodyForce and match their name
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

                force.addParticle(v)

forcefield.parsers["MBPolThreeBodyForce"] = MBPolThreeBodyForceGenerator.parseElement

## @private
class MBPolDispersionForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []

    @staticmethod
    def parseElement(element, forceField):
        generator = MBPolDispersionForceGenerator()
        forceField._forces.append(generator)

        # <MBPolDispersionForce>
        #     <Residue name="HOH" class1="O" class2="H" class3="H" />
        #     <Parameters type1="MBPol-O" type2="MBPol-O" c6="9.92951990e+08" d6="9.29548582e+01" />
        #     <Parameters type1="MBPol-O" type2="MBPol-H" c6="3.49345451e+08" d6="9.77520243e+01" />
        #     <Parameters type1="MBPol-H" type2="MBPol-H" c6="8.40715638e+07" d6="9.40647517e+01" />
        # </MBPolDispersionForce>

        for MBPolDispersionForce_template in element.findall('MBPolDispersionForce'):
            types = forceField._findAtomTypes(MBPolDispersionForce_template.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolDispersionForce_template.attrib['class1'],
                                    MBPolDispersionForce_template.attrib['class2'],
                                    MBPolDispersionForce_template.attrib['class3'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolDispersionForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolDispersionForce()
            force.setCutoff(float(nonbondedCutoff.value_in_unit(unit.nanometer)))
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            atom1 = data.atoms[angle[0]]
            atom2 = data.atoms[angle[1]]
            atom3 = data.atoms[angle[2]]
            if atom1.residue.name == 'HOH' and atom2.residue.name == 'HOH' and atom3.residue.name == 'HOH':
                # FIXME loop through all residues of MBPolDispersionForce and match their name
                v = mbpolplugin.vectori()
                v.push_back(atom2.index)
                v.push_back(atom1.index)
                v.push_back(atom3.index)

                force.addParticle(v)

forcefield.parsers["MBPolDispersionForce"] = MBPolDispersionForceGenerator.parseElement

## @private
class MBPolElectrostaticsForceGenerator:

    def __init__(self):

        self.types1 = []
        self.types2 = []
        self.types3 = []
        self.typeMap = {}

    @staticmethod
    def parseElement(element, forceField):

        # <MBPolElectrostaticsForce>
        #     <Residue class1="OW" class2="HW" class3="HW" thole-charge-charge="0.4" thole-charge-dipole="0.4" thole-dipole-dipole-intermolecules="0.055" thole-dipole-dipole-1-2="0.055" thole-dipole-dipole-1-3="0.626" thole-dipole-dipole-2-3="0.626" /> 
        #     <Atom type="MBPol-O" charge="-5.1966000e-01" damping-factor="0.00131" polarizability="0.00131" />
        #     <Atom type="MBPol-H" charge="2.5983000e-01" damping-factor="0.000294" polarizability="0.000294" />
        #     <Atom type="MBPol-M" charge="0" damping-factor="0.00131" polarizability="0" />
        # </MBPolElectrostaticsForce>

        generator = MBPolElectrostaticsForceGenerator()
        forceField._forces.append(generator)

        for MBPolElectrostaticsForce_template in element.findall('MBPolElectrostaticsForce'):
            types = forceField._findAtomTypes(MBPolElectrostaticsForce_template.attrib, 3)
            if None not in types:

                generator.types1.append(types[0])
                generator.types2.append(types[1])
                generator.types3.append(types[2])

            else:
                outputString = self.__class__ + ": error getting types: %s %s %s" % (
                                    MBPolElectrostaticsForce_template.attrib['class1'],
                                    MBPolElectrostaticsForce_template.attrib['class2'],
                                    MBPolElectrostaticsForce_template.attrib['class3'])
                raise ValueError(outputString)

        for residue in element.findall('Residue'):

        #     <Residue name="HOH" class1="OW" class2="HW" class3="HW" thole-charge-charge="0.4" thole-charge-dipole="0.4" thole-dipole-dipole-intermolecules="0.055" thole-dipole-dipole-1-2="0.055" thole-dipole-dipole-1-3="0.626" thole-dipole-dipole-2-3="0.626" /> 

            name = residue.attrib["name"]

            v = mbpolplugin.vectord()
            thole_components = ["thole-charge-charge", "thole-charge-dipole", "thole-dipole-dipole-intermolecules", "thole-dipole-dipole-OH", "thole-dipole-dipole-HH", ]
            for each in thole_components:
                v.push_back(float(residue.attrib[each]))

            generator.typeMap[name] = dict(thole = v)

        for atom in element.findall('Atom'):
        #     <Atom type="MBPol-H" charge="2.5983000e-01" damping-factor="0.000294" polarizability="0.000294" />
            types = forceField._findAtomTypes(atom.attrib, 1)
            if None not in types:

                for t in types[0]:
                    generator.typeMap[t] = dict(  polarizability = float(atom.attrib['polarizability']),
                                             charge         = float(atom.attrib['charge']),
                                             damping_factor = float(atom.attrib['damping-factor']) )

            else:
                outputString = "MBPolElectrostaticsForceGenerator: error getting type for atom: %s" % (atom.attrib['type'])
                raise ValueError(outputString)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):

        existing = [sys.getForce(i) for i in range(sys.getNumForces())]
        existing = [f for f in existing if type(f) == mbpolplugin.MBPolElectrostaticsForce]

        if len(existing) == 0:
            force = mbpolplugin.MBPolElectrostaticsForce()
            sys.addForce(force)
        else:
            force = existing[0]

        for i in range(len(data.angles)):
            angle = data.angles[i]
            # FIXME loop through all residues of MBPolElectrostaticsForce and match their name

            # FIXME cheating! get virtual site index by max(otheratom indices)  + 1
            global_atoms_indices = set([angle[local_atom_index] for local_atom_index in [0,1,2]])
            virtual_site_index = max(global_atoms_indices) + 1
            global_atoms_indices.add(virtual_site_index)

            for atomIndex in global_atoms_indices:
                atom = data.atoms[atomIndex]
                t = data.atomType[atom]
                if t in self.typeMap:
                    other_atoms = global_atoms_indices - set([atomIndex])
                    force.addElectrostatics(self.typeMap[t]['charge'], data.atoms[other_atoms.pop()].index, data.atoms[other_atoms.pop()].index, data.atoms[other_atoms.pop()].index, self.typeMap[atom.residue.name]['thole'], self.typeMap[t]['damping_factor'], self.typeMap[t]['polarizability'])

                else:
                    raise ValueError('No type for atom %s %s %d' % (atom.name, atom.residue.name, atom.residue.index))

forcefield.parsers["MBPolElectrostaticsForce"] = MBPolElectrostaticsForceGenerator.parseElement
