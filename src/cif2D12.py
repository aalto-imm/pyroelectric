#!/usr/bin/env python 

"""
This function will write a .d12 Inputs for our standard workflows
The function can either read in a cif file (and standardize it) or an ase atoms object.
Options for TESTGEOM, FINDSYM and optimize are avalable (later on also frequency or other)

An option is given to give a template with further keywords which will be read in.
An option is given to add the basis set block for selected atoms (library like mayebe through a dictionary)


"""


from ase.io import cif
from ase.data import chemical_symbols, atomic_numbers
from ase.build import make_supercell
from ase.io.crystal import write_crystal
from ase.spacegroup import *
from ase.spacegroup.symmetrize import refine_symmetry

from copy import deepcopy

from ase import Atoms
import ase


import os  #  to get cwd
import argparse  #  to be able to set variables when calling the script
import numpy as np

import spglib as spg


def correct_annoying_CIFs(cifFile):
        
    dataCIF = []
    
    #  reading the original cif-file, encoding is set to UTF-8 to handle possible errors arising from e.g. greek letters in titles; this will skip all empty lines
    with open(cifFile, 'r', encoding='utf-8') as file:
        for line in file:
            if line != '\n':
                line = line.replace("\n", "")
            
                dataCIF.append(line)
        file.close()
        
    #  The corDataCIF will be used to correct the data later
    corDataCIF = deepcopy(dataCIF)
    
    lineReplace = []
    semikolonMarker = False  #  will be set to True, in case one multi-line item is found
    insertNum = 0
    
    replacedSmthg = False  #  will be set True if anything was changed in the cif-file through the function and used as bool to determine which file is returned
    
    #  parsing through the read data to look for ; ... ; encapsulated items -> corresponding lines will be stored temporarily before being changed in corDataCIF
    for lineNum, line in enumerate(dataCIF):
        
        #  The split() is used to distinguish between semi-colons used for marking an extended line and semi-colons that are for example part of the articles's or journal's title
        spLine = line.split(sep=' ')
        
        if ';' in spLine and semikolonMarker == False:
            replacedSmthg = True
            semikolonMarker = True
            insertNum = lineNum - 1
            lineReplace.append(dataCIF[lineNum-1])
            lineReplace.append(dataCIF[lineNum])
           
            continue

        #  Saves all lines within the ; ... ; block
        if semikolonMarker == True and ';' not in spLine:
            lineReplace.append(line)
            
        #  Finds the end of the code block
        if semikolonMarker == True and ';' in spLine:
            lineReplace.append(line)
            semikolonMarker = False
            
            #  Changes the semi-colon to " to make it readable by ase.io.read_cif()
            for itemNum, item in enumerate(lineReplace):
                splits = item.split(sep=' ')
                
                if item == ';':
                    lineReplace[itemNum] = '"'
                elif splits[0] == ';':
                    lineReplace[itemNum] = lineReplace[itemNum].replace(';', '"')
            
            NewLine = ' '.join(lineReplace)

        #  deletes all lines that will be replaced by a single line argument for the ; ... ; block
            for i in range(insertNum, insertNum+len(lineReplace), 1):               
                corDataCIF[i] = ""
            
            corDataCIF[insertNum] = NewLine
            lineReplace = []  #  clears the list for next ; ... ; block
                
    #  writes the correcte data in a new cif-file with a "_corrected" extension
    if replacedSmthg == True:
        newFileName = cifFile[:-4] + '_corrected' + '.cif'
        
        #print('Saving corrected cif as: %s.\n'%newFileName)                    

        with open(newFileName, 'w', encoding='utf-8') as file:
            for line in corDataCIF:
                if line != '':
                    file.write(line + '\n')
        
        return newFileName  #  returns file name for correted cif
    
    else:
        
        return cifFile  #  if nothing was changes it returns the original file name



  #  this is a copy of the ase method _get_reduced_indices which is needed to get a list of indices, which can be used to get the reduced list of numbers/symbols for the basis of a crystal structure
  
def _get_reduced_indices(atoms, tol: float = 1e-5):# -> List[int]:
    """Get a list of the reduced atomic indices using spglib.
    Note: Does no checks to see if spglib is installed.

    :param atoms: ase Atoms object to reduce
    :param tol: ``float``, numeric tolerance for positional comparisons
    """
    from ase.spacegroup.symmetrize import spglib_get_symmetry_dataset

    # Create input for spglib
    spglib_cell = (atoms.get_cell(), atoms.get_scaled_positions(),
                   atoms.numbers)
    symmetry_data = spglib_get_symmetry_dataset(spglib_cell,
                                                symprec=tol)

    return list(set(symmetry_data.equivalent_atoms))


  #  This function takes an ase atoms object created from a cif file (which is equal to the unit cell created by that cif) and returns the ase atoms object of the standard setting of said unit cell by utilizing spglib

def standardize_atomsAtom_crstalStructure(aseAtom, precision=0.001, verbosity=1, **kwargs):    
    
    #  get the variables to build the spglib cell
    fracCoords = aseAtom.get_scaled_positions()
    cell = aseAtom.get_cell()
    numbers = aseAtom.get_atomic_numbers()
    
    #  get the space group before standardizing; symprec is set to 0.01 so that hexagonal cells can be dealt with (otherwise due to rounding errors and wrong space groups might occurr)
    sg_before = spg.get_spacegroup([cell, fracCoords, numbers], symprec=precision)

    #  use spglib.standardize_cell -> returns spglib cell object with standardized coodrinates and unit cell matrix
    standardCell = spg.standardize_cell([cell, fracCoords, numbers], symprec=precision, **kwargs)
    
    
    #if len(fracCoords) != len(standardCell[1]):  #  if the sizes of the unit cell before and after standardization do not match up (e.g. if the size of the unit cell was decreased due to a different setting), a new ase Atoms object will be created based on the retrived standard cell
    
    #  skip copying the original cell, just dump the standardized cell, positions and numbers in a new object
    newStandAseAtom = Atoms(numbers=standardCell[2], scaled_positions=standardCell[1], cell=standardCell[0], pbc=[True, True, True])
    
    # else:  #  make a deepcopy of the original cell, to change the cell parameters and the scaled_positions aka fractional coordinates
    #     newStandAseAtom = deepcopy(aseAtom)
            
    #     newStandAseAtom.set_cell(standardCell[0])       
    #     newStandAseAtom.set_scaled_positions(standardCell[1])  
    #     newStandAseAtom.set_atomic_numbers(standardCell[2])

    
    #  get the space group after standardizing to check if something is wrong    
    #sg_after = spg.get_spacegroup(standardCell)
    sg_after = spg.get_spacegroup([newStandAseAtom.get_cell(), newStandAseAtom.get_scaled_positions(), newStandAseAtom.get_atomic_numbers()], symprec=precision)
    
    if verbosity > 1:
        print('Spacegroup before standardizing: %s'%sg_before)
        print('Spacegroup after standardizing: %s\n'%sg_after)
    
    #  print an Error that something if something went wrong
    if sg_before != sg_after:
        print('Something went wrong, spacegroups before and after standardizing do not match!\n')
                
    

    
    return newStandAseAtom


def correct_sorting_of_stand_to_prev_ase(aseAtom, stand_aseAtom):
    
    atomSpecies = [chemical_symbols[el] for el in list(dict.fromkeys(aseAtom.numbers))]
    
    xf_sort = []
    num_sort = []
    
    for atom in atomSpecies:
        for numNum, num in enumerate(stand_aseAtom.numbers):
            if atomic_numbers[atom] == num:
                xf_sort.append(stand_aseAtom.positions[numNum])
                num_sort.append(num)
                
                
    sorted_aseAtom = Atoms(positions=xf_sort, numbers=num_sort, cell=stand_aseAtom.cell, pbc=[True, True, True])
    
    
    return sorted_aseAtom


def get_spacegroup_from_primitive(aseAtom):
    
    fracCoords = aseAtom.get_scaled_positions()
    cell = aseAtom.get_cell()
    numbers = aseAtom.get_atomic_numbers()

    primCell = spg.find_primitive([cell, fracCoords, numbers], )

    spacegroup = spg.get_spacegroup(primCell)
    
    part1, part2 = spacegroup.split("(", 1)
    
    sgNum = int(part2.split(")", 1)[0])
    
    return sgNum

def write_crystal_mod(filename, atoms):
    """Method to write atom structure in crystal format
       (fort.34 format)
    """

    ispbc = atoms.get_pbc()
    box = atoms.get_cell()
    
    with open(filename, 'w') as fd:

        # here it is assumed that the non-periodic direction are z
        # in 2D case, z and y in the 1D case.
    
        if ispbc[2]:
            fd.write('%2s %2s %2s %23s \n' %
                     ('3', '1', '1', 'E -0.0E+0 DE 0.0E+0( 1)'))
        elif ispbc[1]:
            fd.write('%2s %2s %2s %23s \n' %
                     ('2', '1', '1', 'E -0.0E+0 DE 0.0E+0( 1)'))
            box[2, 2] = 500.
        elif ispbc[0]:
            fd.write('%2s %2s %2s %23s \n' %
                     ('1', '1', '1', 'E -0.0E+0 DE 0.0E+0( 1)'))
            box[2, 2] = 500.
            box[1, 1] = 500.
        else:
            fd.write('%2s %2s %2s %23s \n' %
                     ('0', '1', '1', 'E -0.0E+0 DE 0.0E+0( 1)'))
            box[2, 2] = 500.
            box[1, 1] = 500.
            box[0, 0] = 500.
    
        # write box
        # crystal dummy
        fd.write(' %.10E %.10E %.10E \n'
                 % (box[0][0], box[0][1], box[0][2]))
        fd.write(' %.10E %.10E %.10E \n'
                 % (box[1][0], box[1][1], box[1][2]))
        fd.write(' %.10E %.10E %.10E \n'
                 % (box[2][0], box[2][1], box[2][2]))
    
        # write symmetry operations (not implemented yet for
        # higher symmetries than C1)
        fd.write(' %2s \n' % (1))
        fd.write(f' {1:.10E} {0:.10E} {0:.10E} \n')
        fd.write(f' {0:.10E} {1:.10E} {0:.10E} \n')
        fd.write(f' {0:.10E} {0:.10E} {1:.10E} \n')
        fd.write(f' {0:.10E} {0:.10E} {0:.10E} \n')
    
        # write coordinates
        fd.write(' %8s \n' % (len(atoms)))
        coords = atoms.get_positions()
        tags = atoms.get_tags()
        atomnum = atoms.get_atomic_numbers()
        for iatom, coord in enumerate(coords):
            fd.write('%5i  %19.10f %19.10f %19.10f \n'
                     % (atomnum[iatom] + tags[iatom],
                        coords[iatom][0], coords[iatom][1], coords[iatom][2]))


  #  This function writes .d12 inputs from ase atoms, which should include the unit cell. To make sure that the unit cell is in standard setting, the ase atoms object will standardized beforehand
  #  Up to this point space groups with S/Z and thus possible origin shifts are not handled properly -> errors or wrong space groups are obtained if the .d12 is run with CRTYSATL -> this can often be avoided by changing the "1 0 0" line to "1 0 1" which in the .d12 so please double check the input.

def write_d12_GEOM_from_aseAtoms(aseAtom, d12_filename='', P1 = False,  numDig= 6, verbosity=1, template=None, add_basisSets=[], add_standard=False, sort=False, external=False, **kwargs):
    
    #  Depending on the setting of the primitive varaible, the unit cell will either be printed out as the primitive cell (True) or the normal unit cell with only the symmetry reduced coordinates (False). The latter will automatically standardize the unit cell according to spglib.
    
    cwd = __file__
    path = os.path.dirname(cwd)    

    
    if P1 == True:
        #  Set the spacegroup to P1 instead of the standard space group
        sg = 1
        sg_ase = ase.spacegroup.Spacegroup(sg)
        
        cell = aseAtom.cell.cellpar()

        
        #  since it's P1, all coordinatates are 1a positions
        basis = aseAtom.get_scaled_positions() 
        
        
        basis_numbers = [i for i in range(0, len(basis), 1)]
        allNumbers = aseAtom.numbers
        StandAseAtom = aseAtom
        
        if verbosity > 0:
            print('P1 symmetry chosen, unit cell will be created directly from the input object.\n')
        

    
    elif P1 == False:
        #  standardize the unit cell
        StandAseAtom = standardize_atomsAtom_crstalStructure(aseAtom, **kwargs)
        
        if verbosity > 0:
            print('Standardizing the unit cell with spglib.standardize().\n')
        
        #  Testing if the symmetry standardisation did something weird and give a warning (although not helpful, I guess)
        # from ase.utils.structure_comparator import SymmetryEquivalenceCheck

        # comp = SymmetryEquivalenceCheck()
        
        # if comp.compare(aseAtom, StandAseAtom) == False:
        #     print(r'Oh no, something went wrong... ¨\ O_Q /¨ \n')

        
        #  getting the spacegroup of the ase atoms object as an ase Spacegroup object (since they already have the right format for the space group name to add to the d12 input)
        sg = ase.spacegroup.symmetrize.check_symmetry(StandAseAtom).number
        sg_ase = ase.spacegroup.Spacegroup(sg)
        
        cell = StandAseAtom.cell.cellpar()
        
        #  get the basis of the ase atoms object and the right atom numbers
        basis = get_basis(StandAseAtom)
        basis_numbers = _get_reduced_indices(StandAseAtom)  #  hidden function from ase, to get idices of the numbers list that hold unique positions    
        allNumbers = StandAseAtom.numbers
        
            
    #  Check which kind of crystal system the spacegroup belongs to, to set which information is to be set in the cell parameter line
    
    if sg <= 2:  #  triclinic
        cellParam = cell
        
    elif sg >= 3 and sg <= 74:  #  monoclinic and orthorhombic
        cellParam = []
        for param in cell:
            if param != 90.:
                cellParam.append(param)
    
    elif sg >= 75 and sg <= 142:  #  tetragonal
        cellParam = []
        for param in cell:
            if param != 90. and param not in cellParam:
                cellParam.append(param)
    
    elif sg >= 143 and sg <= 194:  # trigonal and hexagonal
        cellParam = []
        for param in cell:
            if param not in [90., 120.] and param not in cellParam:
                cellParam.append(param)
    
    elif sg >= 195:  #  cubic:
        cellParam = [cell[0]]
        
    else:
        print('No space group specified or something went wrong, please check cell parameters.\n')
        cellParam = []
        
    
    #  get the standardized space group symbolfrom the ase spacegroup; modify in case of cubic groups with no > 220 to change "-3" to "3"
    #  replace the "e" in space groups 39 and 41 with b; in 64, 67 and 68 with a
    
    if sg > 220 or sg in [200, 201, 202, 203, 204, 205, 206]:
        sg_name = sg_ase.symbol.replace('-3', '3')
    elif sg == 39 or sg == 41:
        sg_name = sg_ase.symbol.replace('e', 'b')
    elif sg in [64, 67, 68]:
        sg_name = sg_ase.symbol.replace('e', 'a')
    
    else:
        sg_name = sg_ase.symbol
        
    
    #  for atoms starting at Rb (37) all number are changed to 2XX so that CRYSTAL can read them properly (core potential possibilities, I guess)
        
    
    fracCoordBlock = []
    
    for num, index in enumerate(basis_numbers):
        curLine = []
        
        if allNumbers[index] > 36:
            curLine.append(allNumbers[index] + 200)
            
        else:
            curLine.append(allNumbers[index])
        
        for coord in basis[num]:
            number = round(coord, numDig)
            curLine.append(number)        
        
        fracCoordBlock.append(curLine)
    
    #  Sort the coordinates by atomic number (lightest first) to make stuff more human readable
    if sort == True:
        fracCoordBlock.sort()    
        
    
    #  check if a filename is set, otherwise set it to the compound name; file ending of .d12 is added if not given
    
    if d12_filename == '':
        d12_filename += StandAseAtom.get_chemical_formula(empirical=True)
        #d12_filename += str(StandAseAtom.symbols)

    if d12_filename.find('.d12') == -1:
        d12_filename += '.d12'
        
    #  write the .d12 file with the geometry block and ENDGEOM/END to close the block (maybe in the Future are other options to set more keywords/options)
    
        
    
    
    with open(d12_filename, 'w') as d12_file:
        d12_file.write('converted and standardized from ase Atoms() object for %s in space group %s\n'%(StandAseAtom.get_chemical_formula(mode='metal'), sg))
        
        if external == True:
            d12_file.write('EXTERNAL\n')
            write_crystal_mod(d12_filename.replace('d12', 'ext'), StandAseAtom)
            
        else:
        
            d12_file.write('CRYSTAL\n')
            d12_file.write('1 0 0\n')
            d12_file.write(sg_name.upper() + '\n')
            for coord in cellParam:
                d12_file.write(str(round(coord, numDig)) + '  ')
            d12_file.write('\n')
            d12_file.write(str(len(fracCoordBlock)) + '\n')
            for coordLine in fracCoordBlock:
                for coordNum, coord in enumerate(coordLine):  #  This is enumerated and packed in an if/else, to properly align stuff in the input, not that anyone would care, but copy-pasting is easier this way...
                    if coordNum == 0:
                        d12_file.write(str(coord) + '\t')
                    else:
                        d12_file.write('{: .{numDig}f}'.format(coord, numDig=numDig) + '\t')
                d12_file.write('\n')
                
        
        if template == None:
            
            d12_file.write('ENDGEOM\n')
            d12_file.write('\n****************************\n')
            d12_file.write('add additional commands here\n')
            d12_file.write('****************************\n\n')
            
            d12_file.write('END')
            
            if verbosity > 0:
                print('Saved .d12 input as %s.\n'%d12_filename)    
            
            return
            
        
        elif template.upper() == 'FINDSYM':
            
            d12_file.write('TESTGEOM\n')
            d12_file.write('FINDSYM\n')
            d12_file.write('ENDGEOM\n')
            d12_file.write('END')
            
            if verbosity > 0:
                print('Saved .d12 input as %s.\n'%d12_filename)    

            return
        
        
        elif template.upper() == 'OPTGEOM':
            
            d12_file.write('OPTGEOM\n')
            d12_file.write('ENDOPT\n')
            d12_file.write('ENDGEOM\n')

        elif template.upper() == 'FREQUENCY':

            d12_file.write('FREQCALC\n')
            d12_file.write('NUMDERIV\n')
            d12_file.write('2\n')
            d12_file.write('ENDFREQ\n')
            d12_file.write('ENDGEOM\n')
        
        elif template.upper() == 'IR' or template.upper() == 'RAMAN':
            
            d12_file.write('FREQCALC\n')
            d12_file.write('RESTART\n')
            d12_file.write('NUMDERIV\n')
            d12_file.write('2\n')
            d12_file.write('INTENS\n')
            d12_file.write('INTRAMAN\n')
            d12_file.write('INTCPHF\n')
            d12_file.write('ENDINTCPHF\n')
            d12_file.write('IRSPEC\n')
            d12_file.write('DAMPFAC\n')
            d12_file.write('8\n')
            d12_file.write('ENDIRSPEC\n')
            d12_file.write('RAMANEXP\n')
            d12_file.write('298.15 632.8\n')
            d12_file.write('RAMSPEC\n')
            d12_file.write('DAMPFAC\n')
            d12_file.write('8\n')
            d12_file.write('VOIGT\n')
            d12_file.write('0.5\n')
            d12_file.write('ENDRAMSPEC\n')
            d12_file.write('ENDFREQ\n')
            d12_file.write('ENDGEOM\n')
        
        elif template.upper() == 'GRADCAL':
            d12_file.write('ENDGEOM\n')      
        
        else:
            try:
                with open(template, 'r') as inFile:
                    lines = [l.split() for l in (line.strip() for line in inFile) if l]  # gets rid of blank lines and "\n" at the end of each line like magic
            except:
                raise ValueError('Invalid TEMPLATE option or filename!\n')
            else:
                for line in lines:
                    for el in line:
                        d12_file.write(el + ' ')
                    d12_file.write('\n')
            
        
        if add_basisSets != []:
            for basis_set in add_basisSets:
                for file in os.listdir(path+'/BasisSets/'):
                    if file.startswith(basis_set[0].lower() + '_') == True and file.find(basis_set[1].lower()) != -1:
                        with open(path+'/BasisSets/' + file, 'r') as basis:
                            for line in basis:
                                if line.split() != []:
                                    d12_file.write(line)
                        
                        d12_file.write('\n')
                        continue

                        
            d12_file.write('99 0 \n')
            d12_file.write('ENDBAS\n')

        
        if add_standard == True:
            d12_file.write('DFT\n')
            d12_file.write('PBE0\n')
            d12_file.write('ENDDFT\n')
            d12_file.write('SHRINK\n')
            
            shrink_a = int(np.ceil(30/cell[0]))
            shrink_b = int(np.ceil(30/cell[1]))
            shrink_c = int(np.ceil(30/cell[2]))
            
            d12_file.write(f'0 {shrink_a}\n')            
            d12_file.write(f'{shrink_a} {shrink_b} {shrink_c}\n')            

            d12_file.write('TOLINTEG\n')
            d12_file.write('8 8 8 8 16\n')
            d12_file.write('MAXCYCLE\n')
            d12_file.write('100\n')
            d12_file.write('FMIXING\n')
            d12_file.write('80\n')
            d12_file.write('EXCHSIZE\n')
            d12_file.write('30000000\n')
            d12_file.write('BIPOSIZE\n')
            d12_file.write('30000000\n')


        if template.upper() in ['FREQUENCY', 'IR', 'RAMAN']:
            d12_file.write('TOLDEE\n')
            d12_file.write('9\n')
            
        if template.upper() == 'GRADCAL':
            d12_file.write('GRADCAL\n')

        
        d12_file.write('END')
        d12_file.close()
    
    
    
    if verbosity > 0:
        print('Saved .d12 input as %s.\n'%d12_filename)    
            
    
    return





if __name__ == '__main__':
    
    print(r'/\___/\   Script to create CRYSTAL .d12input-files from')
    print(r':O: :O:   ase atoms objects.')
    print(r' ¨   ¨')
    
    """
    This script can be used to generate .d12 input files for CRYSTAL from cif files
    crystal structures will be standardized by spglib or processed within spacegroup P1
    
    variables to be set by yhe user:
        cif  -> cif-file with crystal structure
        p1   -> set True if cif should be handled as space group P1
        step -> set the calculation step or give a template to add additional commands to the input
        fn   -> file name of the .d12
    
    """
    
    parser = argparse.ArgumentParser(description='Wizard to create .d12 inputs from cif files', formatter_class=argparse.RawTextHelpFormatter, epilog='This script can handle creating the standard calculation steps of optimization, frequency and IR/RAMAN. Basis sets can be loaded from folder "basisSets".')

    
    parser.add_argument("-cif", help="CIF-file from which the structure is read", nargs='?', type=str)
    
    parser.add_argument("-fn", "--FNAME", help="set the filename of the created .d12 input", nargs='?', default='', type=str)
    
    parser.add_argument("-s", "--STEP", help='set the step of the calculation\n' + 
                        'None: just create the geometry part\n' + 
                        'FINDSYM: to test if the symmetry is right\n' + 
                        'OPTGEOM: geometry optimization\n' + 
                        'FREQUENCY: harmonic frequencies\n' +
                        'IR/RAMAN: theoretical IR and RAMAN sprectra\n' +
                        'GRADCAL: forces and potential energy calculation', type=str, nargs='?', default=None)
    
    parser.add_argument("-p1", "--P1", help='option to handle the crystal structure in space group P1', type=bool, default=False)
    
    args, unknown = parser.parse_known_args()
    
    
    try:
        compound = cif.read_cif(args.cif)
    except:
        new_cif = correct_annoying_CIFs(args.cif)
        compound = cif.read_cif(new_cif)
    
    elements = set(compound.get_chemical_symbols())
    
    basis_set_input = []
    
    for el in elements:
        if el in ['Li', 'Na', 'K', 'Rb', 'Cs', 'Mg', 'Ca', 'Sr', 'Ba']:
            basis_set_input.append([el, 'SVP'])
            
        else:
            basis_set_input.append([el, 'TZVP'])
                
    if args.STEP == None:
        add_std = False
    else:
        add_std = True

        
    if args.STEP == None:
        
        write_d12_GEOM_from_aseAtoms(compound, d12_filename=args.FNAME, template=args.STEP)
        
    else:
    
        write_d12_GEOM_from_aseAtoms(compound, d12_filename=args.FNAME, template=args.STEP, add_basisSets=basis_set_input, add_standard=add_std, P1=args.P1)
    


    
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

    
    
















    

    

