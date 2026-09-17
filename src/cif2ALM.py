"""
Definition of classes for each ALM input block
Basic idea: For each input block all possible tags are stored in a class, where objects can be created for each input block; Most attributes are set as NONE as default
-> To be able to set each attribute individually when creating an input. Some are set by functions availabe for ase ATOM objects

Input blocks for cutoff, cell and positions are created by individual functions (below) since they do not require selection of keywords


"""




from ase.io import cif
from ase.data import chemical_symbols, atomic_numbers
from ase.build import make_supercell
from ase.io.crystal import write_crystal
from ase.units import Bohr

import os  #  to get cwd
import argparse  #  to be able to set variables when calling the script
import numpy as np
import sys

from ase.io.jsonio import write_json, read_json



#  to create the dimesnion matrix needed for the supercell creation
def create_supercell_matrix(dimList):
    dimMatrix = []
    for num, dim in enumerate(dimList):
        vec = np.zeros(3)
        vec[num] = dim
        dimMatrix.append(vec)
    npMatrix = np.array(dimMatrix)
    return npMatrix




#  class for the &general input block

class ALMgeneral():
    
    def __init__(self, aseAtom, **kwargs):
        #List of all variables that can be set for the general block input
        #descriptions taken from ALAMODE webpage
        self.PREFIX = kwargs.get('PREFIX', None)  #  Job prefix to be used for names of output files
        self.MODE = kwargs.get('MODE', None)  #  options: optimize | suggest
        self.NAT = len(aseAtom.positions)  #  Number of atoms in the supercell -> created by atoms object
        self.NKD = len(set(aseAtom.numbers))  #  Number of atomic species
        #self.KD = [chemical_symbols[el] for el in set(aseAtom.numbers)]  #  Chemical symbols for each atom species; due to set() method they are ordered according to the alphabet...
        self.KD = [chemical_symbols[el] for el in list(dict.fromkeys(aseAtom.numbers))] 
        self.TOLERANCE = kwargs.get('TOLERANCE', None)  #  TOLERANCE-tag : Tolerance for finding symmetry operations
        self.PRINTSYM = kwargs.get('PRINTSYM', None)  #  PRINTSYM-tag = 0 | 1
        self.FCSYM_BASIS = kwargs.get('FCSYM_BASIS', None)  #  FCSYM_BASIS-tag = Cartesian | Lattice
        self.MAGMOM = kwargs.get('MAGMOM', None)  #  List of magnetic moments
        self.NONCOLLINEAR = kwargs.get('NONCOLLINEAR', None)  #  0 | 1 if 1: one can set MAGNOM with noncolliniar meganetic structure
        self.PERIODIC = kwargs.get('PERIODIC', None)  #  PERIODIC[1], PERIODIC[2], PERIODIC[3] -> array if one item is 0 periodic boundary condition is turned off along that direction
        self.NMAXSAVE = kwargs.get('NMAXSAVE', None)  #  The maximum order of anharmonic force constants printed out in PREFIX.xml
        self.HESSIAN = kwargs.get('HESSIAN', None)  #  0 | 1 if 1 Hessian matrix is saved
        self.FC3_SHENGBTE = kwargs.get('FC3_SHENGBTE', None)  #  0 | 1 if 1 3rd order force constants saved in ShengBTE code
        self.FC_ZERO_THR = kwargs.get('FC_ZERO_THR', None)  #  Threshold value used when trimming force constants when creating PREFIX.xml
        
    def get_all_att(self):

        List = vars(self)

        return List
    

    

#  class for the &interaction input block
class interaction():
    
    def __init__(self, **kwargs):
        self.NORDER = kwargs.get('NORDER', None)  #  The order of force constants to be calculated. Anharmonic terms up to (m+1) th order will be considered with NORDER = m
        self.NBODY = kwargs.get('NBODY', None)  #  Entry for excluding multiple-body interactions from anharmonic force constants
        
    def get_all_att(self):

        List = vars(self)

        return List
      

#  class for the &optimize input block
class optimize():
    
    def __init__(self, **kwargs):
        self.LMODEL = kwargs.get('LMODEL', None)  #  Choice of the linear model used for estimating force constants
        self.DFSET = kwargs.get('DFSET', None)  #  File name containing displacement-force datasets for training
        self.NDATA = kwargs.get('NDATA', None)  #  Number of displacement-force data sets
        self.NSTART = kwargs.get('NSTART', None)
        self.NEND = kwargs.get('NEND', None)  #  Specifies the range of data to be used for training
        self.SKIP = kwargs.get('SKIP', None)  #  Specifies the range of data to be skipped for training
        self.ICONST = kwargs.get('ICONST', None)  #  0 | 1 | 2 | 3 | 11  set constrains
        self.PERIODIC_IMAGE_CONV = kwargs.get('PERIODIC_IMAGE_CONV', None)  #  0 | 1 whether to consider periodic images for constarints
        self.ROTAXIS = kwargs.get('ROTAXIS', None)  #  Rotation axis used to estimate constraints for rotational invariance. This entry is necessary when ICONST = 2, 3.
        self.FC2XML = kwargs.get('FC2XML', None)  #  XML file to which the harmonic terms are fixed upon training
        self.FC3XML = kwargs.get('FC3XML', None)  #  XML file to which the cubic terms are fixed upon training
        self.SPARSE = kwargs.get('SPARSE', None)  #  0 | 1 to set sparse solver use 1
        self.SPARSESOLVER = kwargs.get('SPARSESOLVER', None)  #  Type of the sparse solver to use
        self.MAXITER = kwargs.get('MAXITER', None)  #  Number of maximum iterations in iterative algorithms
        self.CONV_TOL = kwargs.get('CONV_TOL', None)  #  Convergence criterion of iterative algorithms
        self.L1_RATIO = kwargs.get('L1_RATIO', None)  #  The ratio of the L1 regularization term
        self.L1_ALPHA = kwargs.get('L1_ALPHA', None)  #  The coefficient of the L1 regularization term
        self.CV = kwargs.get('CV', None)  #  Cross-validation mode for elastic net
        self.DFSET_CV = kwargs.get('DFSET_CV', None)  #  File name containing displacement-force datasets used for manual cross-validation
        self.NDATA_CV = kwargs.get('NDATA_CV', None)  #  Number of displacement-force validation datasets
        self.NSTART_CV = kwargs.get('NSTART_CV', None)
        self.NEND_CV = kwargs.get('NEND_CV', None)  #  Specifies the range of data to be used for validation
        self.CV_MINALPHA = kwargs.get('CV_MINALPHA', None)
        self.CV_MAXALPHA = kwargs.get('CV_MAXALPHA', None)
        self.CV_NALPHA = kwargs.get('CV_NALPHA', None)  #  Options to specify the L1_ALPHA values used in cross-validation
        self.STANDARDIZE = kwargs.get('STANDARDIZE', None)  #  0 | 1 set 1 to standardize each column for the sensing matrix
        self.ENET_DNORM = kwargs.get('ENET_DNORM', None)  #  Normalization factor of atomic displacements
        self.SOLUTION_PATH = kwargs.get('SOLUTION_PATH', None)  #  0 | 1 set 1 to save solution path
        self.DEBIAS_OLS = kwargs.get('DEBIAS_OLS', None)  #  0 | 1 set 1 to collect only non-zero coefficients and fit again before saving
        self.STOP_CRITERION = kwargs.get('STOP_CRITERION', None)  #  The scan over L1_ALPHA stops when the cross-validation score keeps increasing in STOP_CRITERION consecutive steps
        
    def get_all_att(self):
        
        List = vars(self)
        
        return List



#  function to define the cutoff block
#  cutoffs can be set as list of lists, with first and second element of interaction as seperate items eg. sepcCutoff_anharm = ['Ti', 'Sr', 5, 5] will only reduce cutoff for Ti-Sr interaction
def set_cutoff_block(aseAtom, specCutoffs_harm = [], specCutoffs_anharm = []):
    cutoffInput = []
    
    #  first get list with all possible interactions
    interactionList = []
    
    #  get the elements from aseAtoms
    #atomKinds = [chemical_symbols[el] for el in set(aseAtom.numbers)]
    
    atomKinds = [chemical_symbols[el] for el in list(dict.fromkeys(aseAtom.numbers))]
    
    for elNum, element in enumerate(atomKinds):
        for curEl in atomKinds[elNum:]:
            interactionList.append(str(element) + '-' + str(curEl))
    
    for interaction in interactionList:
        curInteraction = []
        spec_harm = []
        spec_anharm = []
        curInteraction.append(interaction)
        
        if specCutoffs_harm != []:
            for specHarm in specCutoffs_harm:
                if specHarm[0] + '-' + specHarm[1] in interaction or specHarm[1] + '-' + specHarm[0] in interaction:
                    #print('Setting ' + str(interaction))
                    spec_harm.append(specHarm[2])
        
        if spec_harm != []:
            curInteraction.append(spec_harm[0])
        else:
            curInteraction.append('None')
            
        
        if specCutoffs_anharm != []:
            for specAnharm in specCutoffs_anharm:
                if specAnharm[0] + '-' + specAnharm[1] in interaction or specAnharm[1] + '-' + specAnharm[0] in interaction:
                    #print('Setting ' + str(interaction))
                    spec_anharm.append(specAnharm[2])
                    spec_anharm.append(specAnharm[3])
                    
        if spec_anharm != []:
            curInteraction.append(spec_anharm[0])
            curInteraction.append(spec_anharm[1])
        else:
            curInteraction.append(8)
            curInteraction.append(8)
                    
        cutoffInput.append(curInteraction)
    
    return cutoffInput
   
    
#  get the cell from ase atoms object, converts Ångström to bohr
#  multiplication factor for the cell can be set optionally
def set_cell_block(aseAtom, multFactor = 1):
    cell = [[multFactor]]
    
    for line in aseAtom.cell:
        curLine = [round(a/Bohr, 10) for a in line]
        cell.append(curLine)
    
    return cell


#  Getting fractional coordinates from ase atoms object
#  option to set number of digits to round
#  coordinates will be sorted by element, sorted by order that is given for the KD keyword
#  numDig can be set optionally, it represents the number of digits the coordinates are rounded by. The default of 6 should be accurate eneough to match the accuracy of possible experimental values
def set_positions_block(aseAtom, numDig = 10):
    fracCoord = aseAtom.get_scaled_positions()
    #atomKind = [chemical_symbols[el] for el in set(aseAtom.numbers)]  #  Same command as for the KD comand in general block
    
    atomKind = [chemical_symbols[el] for el in list(dict.fromkeys(aseAtom.numbers))]
        
    #sort atoms by numbers and combine atom kind number and fractional coordinates
    sortedCoordList = []
    
    atomKindNum = []

    for pos in aseAtom.numbers:
        name = chemical_symbols[pos]
        for elNum, element in enumerate(atomKind):
            if name == element:
                atomKindNum.append(elNum+1)

    for num, element in enumerate(atomKind):
        for line, atomNum in enumerate(atomKindNum):
            if num+1 == atomNum:
                curLine = []
                curLine.append(atomNum)
                for coord in fracCoord[line]:
                    curLine.append(round(coord, numDig))
                sortedCoordList.append(curLine)
    
    return sortedCoordList


  #  the parameter modCPatoms is used to modify the atom numbers for those that normally have a core potential in the basis set -> external geometry file is opened and modified accordingly
def write_CRYSTAL_supercell(aseAtom, modCPatoms=True):
    #since ALAMODE CRYSTAL interface needs one output with only the supercell create an external geometry .ext and a .d12 input to quickly run with runcrys
    write_crystal('supercell.ext', aseAtom)

    #  check the external file for atoms with higher atomic numbers of Rb and above to change their number sto 
    with open('supercell.ext', 'r') as file:
        data = file.readlines()
        newData = []

        for lineNum, line in enumerate(data):
            line = line.split()
            
            if len(line) > 1 and line[0].isdigit() == True and int(line[0]) > 36.:
                line[0] = str(int(line[0])+200)

            newData.append(line)
        file.close()
        
        
        #  This will sort the coordinates in the same way as the write_aLM_input() later on does, so that the crystal parser can read the right coordinates
        newDataSorted = []
        foundCoords = False        
        #AtomList = [chemical_symbols[el] for el in set(aseAtom.numbers)]   
        AtomList = [chemical_symbols[el] for el in list(dict.fromkeys(aseAtom.numbers))]
        
        for lineNum, line in enumerate(newData):
            if foundCoords == False:
                if line[0] == str(len(aseAtom.positions)):
                    foundCoords = True
                    coordStartNum = lineNum + 1
                newDataSorted.append(line)
          
                    
        for atom in AtomList:
            ATnum = atomic_numbers[atom]

            if ATnum > 36:
                ATnum += 200

            for line in newData[coordStartNum::]:
                if int(line[0]) == ATnum:
                    newDataSorted.append(line)                    
            
        
    #  writing the supercell.ext again with modified atomic numbers
    with open('supercell.ext', 'w') as file:
        for line in newDataSorted:
            file.writelines(' '.join(line))
            file.write('\n')
        file.close()
    
    #quick and dirty get .d12 for runcrys / needs same name as external structure file
    outCrystalInput = 'supercell.d12'            
    
    with open(outCrystalInput, 'w') as file:
        file.writelines('Created with ALAMODE input generator, please run with runcrys and rename output to .o\n')
        file.writelines('EXTERNAL\n')
        file.writelines('TESTGEOM\n')
        file.writelines('END\n')
        file.close()    

    
    print('Creating .ext and .d12 for supercell, please run with runcrys and rename output to .o')
    print()
    
    return


  #  this functions writes an ALM input based on one ASE atoms object
  #  by default it adds the &cell and &positions based on the passed object and sets the cutoff parameters with default values
  
""" maybe add option to load the supercell from a json file so test if the aseAtom is an ase atoms object or a json file"""

def write_ALM_input(aseAtom, generalInput=None, interactionInput=None, optimizeInput=None, Cutoffs_harm=[], Cutoffs_anharm=[], multFactor=1, numDig=10, fileName='', **kwargs):
    
    #  function that writes the ALM input based on one ASE atoms object and instances of the general and interaction class, that hold the corresponding keywords for the calculation
    #  optional possibility to set a filename
    
    if fileName != '' and fileName.endswith('.in'):
        ALMinput = fileName
    elif fileName != '' and fileName.endswith('.in') == False:
        ALMinput = fileName + '.in'
    else:
        ALMinput = 'ALMinput.in'
        
    with open(ALMinput, 'w') as file:
        
        if generalInput != None:
            file.write('&general\n')
            
            #  get all overview of all attributes of the general object
            generalCom = generalInput.get_all_att()
    
            #  iterate over all attributes of the general object, only consider those that were set (thus are not NONE), get their value and print it
            for com in generalCom:
                if generalCom.get(com) != None:
                    file.write(' ' + com + ' = ')
                    
                    if type(generalCom.get(com)) == list:
                        for val in generalCom.get(com):
                            file.write(str(val) + ' ')
                    else:
                        file.write(str(generalCom.get(com)))
                    file.write('\n')
                    
            file.write('/\n\n\n')
        
        
        if interactionInput != None:
            file.write('&interaction\n')
            
            interactionCom = interactionInput.get_all_att()
            
            for com in interactionCom:
                if interactionCom.get(com) != None:
                    file.write(' ' + com + ' = ')
                    
                    if type(interactionCom.get(com)) == list:
                        for val in interactionCom.get(com):
                            file.write(str(val) + ' ')
                    else:
                        file.write(str(interactionCom.get(com)))
                    file.write('\n')
            
            file.write('/\n\n\n')
            
            
        if optimizeInput != None:
            file.write('&optimize\n')
            
            optimizeCom = optimizeInput.get_all_att()
            
            for com in optimizeCom:
                if optimizeCom.get(com) != None:
                    file.write(' ' + com + ' = ')
                    
                    if type(optimizeCom.get(com)) == list:
                        for val in optimizeCom.get(com):
                            file.write(str(val) + ' ')
                    else:
                        file.write(str(optimizeCom.get(com)))
                    file.write('\n')
            
            file.write('/\n\n\n')
        
        
        file.write('&cutoff\n')
        cutoffInput = set_cutoff_block(aseAtom, specCutoffs_harm=Cutoffs_harm, specCutoffs_anharm=Cutoffs_anharm)
        
        # from list generated for cutoff get the lines by parsing throug them
        
        for line in cutoffInput:
            file.write(' ')
            for el in line:
                file.write(str(el) + ' ')
            file.write('\n')
        
        file.write('/\n\n\n')
        
        
        
        file.write('&cell\n')
        cellInput = set_cell_block(aseAtom, multFactor=multFactor)
        
        for line in cellInput:
            file.write(' ')
            for el in line:
                #file.write(str(el) + ' ')
                file.write("%.10f" % round(el, 10) + ' ')
            file.write('\n')
        
        file.write('/\n\n\n')

        
        
        
        file.write('&position\n')
        positionInput = set_positions_block(aseAtom, numDig = numDig)
        
        for line in positionInput:
            file.write(' ')
            file.write(str(line[0]))
            file.write(' ')
            for el in line[1::]:
                #file.write(str(el) + ' ')
                file.write("%.10f" % round(el, 10) + ' ')
            file.write('\n')
        
        file.write('/\n\n\n')
        file.close()

        
    print('ALM input saved as %s\n'%ALMinput)
    
    return


#  loop to sepparate additional arguments passed to the script that could be recognized by any of the blocks
def process_additional_arguments(unknownArgs):
    kwargs = {}
    key = ''
    curAttValues = []
    
    for el in unknownArgs:
        if el.find('--') != -1:
            if key != '':
                kwargs[key] = curAttValues
            curAttValues = []
            key = el.replace('--', '')
        else:
            curAttValues.append(el)

    if key != '':
        kwargs[key] = curAttValues
     
    return kwargs

    #  Small function to check for possible fc_X_xml-files so that they can get passed as varaibles to the input automatically
    #  Maybe replace the -fc2 and -fc3 variables with just a "if set invoke the function and look for the right files; gives warning if no file was found"
    #  fcs2 and fcs3 manage whether warnings should be displayed if multiple files with harmonic or cubic terms are found. To display warnings set 'True'


def get_fcsXML(fcs2 = False, fcs3 = False, fcs4 = False):
        
    cwd = os.getcwd()
    
    FC2XML = ''
    FC3XML = ''
    FC4XML = ''
    
    #  Loop through all files in that directory, ignores any files in subdirectory
    for file in os.listdir(cwd):
        if file.endswith('.xml'):  #  only get the .xml files
            foundFCS2 = False
            foundFCS3 = False
            foundFCS4 = False

            with open(file) as file:  #  open each file and search for FC'X' in each line to determnine whether they have harmonic, cubic or quartic terms
                while all(f for f in [foundFCS2, foundFCS3, foundFCS4]) == False:
                    line = file.readline()
                    if line == file.readline(-1):
                        break
                    
                    if line.find('FC2') != -1:
                        foundFCS2 = True
                        
                    elif line.find('FC3') != -1:
                        foundFCS3 = True
                    
                    elif line.find('FC4') != -1:
                        foundFCS4 = True
                
                if foundFCS2 == True and foundFCS3 == False and foundFCS4 == False and FC2XML == '' and fcs2 == True:
                    FC2XML = file.name
                    print('File %s set as value for FC2XML.\n'%FC2XML)

                
                elif foundFCS2 == True and foundFCS3 == False and foundFCS4 == False and FC2XML != '' and fcs2 == True:  #  only the first file found will be saved, if there is more than one, a it will be printed in the console
                   print('Found more than one file for harmonic terms, please check input, %s was selected.\n'%FC2XML)

                elif foundFCS2 == True and foundFCS3 == True and foundFCS4 == False and FC3XML == '' and fcs3 == True:
                    FC3XML = file.name
                    print('File %s set as value for FC3XML.\n'%FC3XML)

                
                elif foundFCS2 == True and foundFCS3 == True and foundFCS4 == False and FC3XML != '' and fcs3 == True:
                    print('Found more than one file for cubic terms, please check input, %s was selected.\n'%FC3XML)
                    
                elif foundFCS2 == True and foundFCS3 == True and foundFCS4 == True and FC4XML == '' and fcs4 == True:
                    FC4XML = file.name
                    print('File %s set as value for FCSXML.\n'%FC4XML)

                    
                elif foundFCS2 == True and foundFCS3 == True and foundFCS4 == True and FC4XML != '' and fcs4 == True:
                    print('Found more than one file for quartic terms, please check input, %s was selected.\n'%FC4XML)
                    
                
            file.close()
        
    if FC2XML == '' and fcs2 == True:
        print('No .xml file with harmoic terms found.\n')
        FC2XML = '!! add harmonic .xml !!'
    if FC3XML == '' and fcs3 == True:
        print('No .xml file with cubic terms found.\n')
        FC3XML = '!! add cubic terms .xml !!'
    if FC4XML == '' and fcs4 == True:
        print('No .xlm file with quartic terms found\n')
        FC4XML = '!! add quartic terms .xml !!'
    
    
    return FC2XML, FC3XML, FC4XML  #  Return the file names, so that they can be used (or are left blank if nothing is found)


  #  get keywords defined in a template file (option to set the name of the template for diffenrent steps) (since each keyword is unique to a block within ALM or ANPHON, there would be no need to add the blocks keywords)
  
def get_keywords_from_template(tempFile):
    
    with open(tempFile, 'r') as inFile:
        lines = [l.split() for l in (line.strip() for line in inFile) if l]  # gets rid of blank lines and "\n" at the end of each line like magic
        inFile.close()
    
    keys = {}
    
    for line in lines:
        if line[0].find('&') == -1:
            keyword = line[0]
            attributes = line[2::]
            
            keys[keyword] = attributes
                
    
    return keys

   
    
   
if __name__ == '__main__':

    print(r'/\___/\ ')
    print(r':O: :O:   Script to create ALM.in input-files from cif-files.')
    print(r' ¨   ¨')


    """
    
    This script is designed for our use, so if other parameters are needed or other workflows are used, please customize the script or just import the classes to write own file

    variables to be set by the user:
        optional:
            cif -> cif-file with (optimized) crystal structure, this will also store the cell and supercell as json files, so setting the cif file is only neccesary for the first step
            fn  -> file-name of the output
            s   -> user specific step (in our case the standard workflow for each step)
            p   -> prefix to add to created files by ALAMODE
            fc2 -> file with stored harmonic terms (mode 2 and 3)
            fc3 -> file with stored cubic terms (mode 3)
            d   -> dimension of the supercell
    """
    

    parser = argparse.ArgumentParser(description='Wizard to create ALM inputs from cif files', formatter_class=argparse.RawTextHelpFormatter, epilog='Additional keywords of ALM variables can be given, they will be added to the right blocks automatically.')
    
    parser.add_argument("-cif", help="CIF-file from which the structure is read, saves the unit cell and created supercell using json; If not set cell and supercell are loaded if available", nargs='?', default='', type=str)
    
    parser.add_argument("-fn", "--FNAME", help="set the filename of the created ALM input", nargs='?', default='', type=str)
    
    parser.add_argument("-s", "--STEP", help='set the step of the calculation\n' + 
                        '0: create ALM0.in with default parameters\n' + 
                        '1: create ALM1.in with default parameters\n' + 
                        '2: create ALM2.in with default parameters (set FC2XML if not in directory)\n' + 
                        '3: create ALM3.in with default parameters (set FC2XML and/or FC3XML if not in directory)', type=int, nargs='?', default=0)
    
    parser.add_argument("-p", "--PREFIX", nargs='?', type=str, help="set the prefix for ALM output files", default='prefix')
    
    parser.add_argument("-fc2", "--FC2XML", help='xml file with stored harmonic terms', nargs='?', type=str, default='')
    
    parser.add_argument("-fc3", "--FC3XML", help='xml file with stored cubic terms', nargs='?', type=str, default='')
    
    parser.add_argument("-d", "--DIM", help="set the dimension of the created supercell, default is [2,2,2]", nargs='+', type=int)
    
    parser.add_argument("-temp", "--TEMPLATE", help="give a template that includes additional arguments, they will be added to the corresponding blocks (&block arguments are ignored, but can be included for better fromating of the template)", nargs='?', type=str, default='')

    
    
    #  parse the arguments and save all unknown attributes in a dictionary that gets passed to the instances of the objects
    args, unknown = parser.parse_known_args()
    kwargs = process_additional_arguments(unknown)  

    
    ###########################################################################
    
    #  get possible additional keywords from a template (change name to variable)
    if args.TEMPLATE != '':
        add_keys = get_keywords_from_template(args.TEMPLATE)
        kwargs = add_keys | kwargs
    
    
    #  check if a cif is set and create aseAtoms objects from that cif
    #  if left empty, search for saved ase atoms objects in directory; exit the script if no objects are found
    
    if args.cif != '':
        cifIn = args.cif
        compound = cif.read_cif(cifIn)        
        write_json(open('cell', 'w'), compound)
    
        #use dimension variable to create matrix
        if args.DIM != None:
            matrix = create_supercell_matrix(args.DIM)
            print('Supercell set to custom %s dimension'%args.DIM)
            print()
        else:
            matrix = create_supercell_matrix([2,2,2])
            print('Supercell set to [2,2,2] (default) dimension')
            print()
        
        #  create supercell
        compound_sc = make_supercell(compound, matrix, order='atom-major')
        
        write_json(open('supercell', 'w'), compound_sc)
        
    
    else:
        try:
            compound = read_json(open('cell', 'r'))
            compound_sc = read_json(open('supercell', 'r'))
            print('Cell and supercell loaded from previous step.\n')
        
        except FileNotFoundError:
            print('No saved cell and/or supercell found.\n')
            print('Please rerun script with -cif option set!')
            sys.exit()
    

    #  Look for FCSXML files in directory and pass them to kwargs, in case -fcs2 and/or -fcs3 are not set
    
    if args.STEP == 2:
        if args.FC2XML == '':
            
            print('Searching for file with harmonic terms in directory.\n')
            
            FCS2, FCS3, FCS4 = get_fcsXML(fcs2=True)
            args.FC2XML = FCS2
            
    if args.STEP == 3:
        if args.FC2XML == '':
            
            print('Searching for file with harmonic terms in directory.\n')

            FCS2, FCS3, FCS4 = get_fcsXML(fcs2=True)
            args.FC2XML = FCS2
            
        if args.FC3XML == '':
            
            print('Searching for file with cubic terms in directory.\n')
            
            FCS2, FCS3, FCS4 = get_fcsXML(fcs3=True)
            args.FC3XML = FCS3
    
    #  This block defines the instances of the classes with which the general, interaction and optimize block are generated. Change this if other defaults are needed.
    #  The mode parameter sets the calculation step
        #  0 -> estimate interatomic force constants by least squares fitting
        #  1 -> optimizing of harmonic force constants
        #  2 -> optimizing of cubic force constants 
        #  3 -> optimizing of quartic force constants
    
    #  In the future it should be possible to give a template for the commands to set in the each block
    
    if args.STEP == 0:
        
        if args.FNAME == '':
            args.FNAME = 'ALM0.in'
        
        #  This will produce an .ext and .d12 file to run with runcrys once, so that the CRYSTAL ALAMODE interface can get the appropriate information
        write_CRYSTAL_supercell(compound_sc)

        
        generalIn = ALMgeneral(compound_sc, PREFIX = args.PREFIX, MODE = 'suggest', **kwargs)
        
        interactionIn = interaction(NORDER='1', **kwargs)
        
        
        #  This part determines whether kwargs has gotten a keyword from the optimize block
        #  If the bool is set to true an optimize block is included for the write_ALM_input function        
        optimizeIn = optimize(**kwargs)

        optList = optimizeIn.get_all_att()

        oneOPTset = False        
        for key in optimizeIn.get_all_att():
            if optList.get(key) != None:
                oneOPTset = True
        
        if oneOPTset == True:  #  is set True if one keyword in the bash line is from the optimize block
            write_ALM_input(compound_sc, generalIn, interactionIn, optimizeIn, fileName = args.FNAME)
        
        else:
            write_ALM_input(compound_sc, generalIn, interactionIn, fileName = args.FNAME)        
        
    elif args.STEP == 1:
        
        if args.FNAME == '':
            args.FNAME = 'ALM1.in'

        
        generalIn = ALMgeneral(compound_sc, PREFIX = args.PREFIX, MODE = 'opt', **kwargs)
        
        interactionIn = interaction(NORDER='1', **kwargs)
        
        optimizeIn = optimize(DFSET = 'DFSET_harmonic', **kwargs)
        
        write_ALM_input(compound_sc, generalIn, interactionIn, optimizeIn, fileName = args.FNAME)
        
    elif args.STEP == 2:
        
        if args.FNAME == '':
            args.FNAME = 'ALM2.in'

        
        if args.FC2XML == '':
            print('Please pass filename for stored harmonic terms')
        
        generalIn = ALMgeneral(compound_sc, PREFIX = args.PREFIX, MODE = 'opt', **kwargs)
        
        interactionIn = interaction(NORDER='2', **kwargs)
        
        optimizeIn = optimize(DFSET = 'DFSET_random', FC2XML = args.FC2XML, **kwargs)
                
        write_ALM_input(compound_sc, generalIn, interactionIn, optimizeIn, fileName = args.FNAME)  
    
    elif args.STEP == 3:
        
        if args.FNAME == '':
            args.FNAME = 'ALM3.in'

        
        if args.FC2XML == '':
            print('Please pass filename for stored harmonic terms')
            print()
        
        if args.FC3XML == '':
            print('Please pass filename for stored cubic terms')
            print()
        
        generalIn = ALMgeneral(compound_sc, PREFIX = args.PREFIX, MODE = 'opt', **kwargs)
        
        interactionIn = interaction(NORDER='3', **kwargs)
        
        optimizeIn = optimize(DFSET = 'DFSET_random', FC2XML = args.FC2XML, FC3XML = args.FC3XML, **kwargs)
                
        write_ALM_input(compound_sc, generalIn, interactionIn, optimizeIn, fileName = args.FNAME)
        
    else:
        print('invalide mode selected')
    
    
    
    
    
    
    #  End of code gimmick
    print()
    print()  
    print(r'      /\       /\ ')
    print(r'     /  \_____/  \      .')
    print(r'    . :::     ::: .    .:.')
    print(r'   : ::O::: :::O:: :  .:::.')
    print(r'   :  :::  o  :::  : .......')
    print('   mmm          mmm   ..... ')
    print('---------------------------------')
    print('FINISHED')


    
    
    
    
    
    
    
    
    """
        idea: one can give a template that loads the needed options for ALM0, ALM1 or ALM2; type of file can be specified by one input command or one programs a standard input themselve
        
        for us: give Template with all needed variables for ALM0 to ALM2
        loop over classes and add cell & positions part to write input
        
        if nothing is set, it just creates all input blocks with commands empty/comment that one can add the commands here
    
    """
    














