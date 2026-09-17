#!/usr/bin/env python 

"""
Script to handle the workflow of getting pyroelectricity constants at different temperatures

The script clusters the steps in the workflow such that in between CRYSTAL, alm or anphon commands need to be run on the cluster or console. Each step is called with the -s attribute and the number of the step (plotting steps are handled sepparately). File names for the CRYSTAL inputs of the displacements are handled automatically, for ALM and ANPHON a file_name might be set. For the second part that handles the modes at different temperatures the file names need to incorporate both of those so that they can be found. If no prefix is set they will always consist of the compound name (given by the unit cell read via read_cif()) and a step specific suffix. Files that are essental for each step, that were created by prior steps do not need to be set as attributes, since they will be found in the directory. Those include:
    .xml files
    .pattern_HARMONIC
    CRY_TEMP or CRY_TEMP_PYRO templates
    .evec files

The following will give a rundown on what happens in each step and what attributes need to be set by the user

Step 0:
    load the stucture form a cif file (ideally an optimized cif)
    make a supercell (either custom dimension set with -dim or default 2x2x2 sized)
    save cell and supercell ase.Atoms() objects as json
    set to True to standardize the structure via spg-lib before use
    
    -> if the structure is not available via a cif file, save it via json from an aseAtoms() object as 'cell'

    Arguments:
    -cif:   cif-file containing the structure
    -dim:   set custom supercell dimension
    -std:   standardize the crystal structure before creating the cell/supercell 
    
    
Step:1 
    write the ALM0.in file, which needs to be run with ALM to get the Harmonic pattern file

    Arguments:
    -p:     prefix for the .pattern_HARMONIC
    -fn:    file name for the ALM input
    
-> in the console run 'alm ALM0.in > ALM0.log'    
    
    

Step 2:
    create the harmonic displacements according to the -pattern_HARMONIC
    calls displace() from ASE_Alamode_interface for disp_mode = pattern_file and set magnitude
    looks for .pattern_HARMONIC automatically in directory
    generates the geometry block of a .d12 file for each displacement (additional commands may be read from a CRY_TEMP file within the same directory)
    generated aseAtoms() objects of the harmonic displacements are saved via json
    
    Arguments:
    -mag:   magnitude of the displacement (default 0.01 Å)
    
-> run all GRADCAL calculations with CRYSTAL to get the forces
    

Step 3:
    extract the forces from the harmonic displacements via extract()
    writes the DFSET_harmonic file that alm can read
    write the ALM1.in file to calculate the harmonic forces

    Arguments:
    -p:     prefix for the .xml file containing the harmonic forces
    -fn:    file name for the ALM input
    
-> in the console run 'alm ALM1.in > ALM1.log'


    
Step 31: 
    write the phband.in to get the .band to plot the phonon dispersion

    Arguments:
    -ks:    k-point source (ase (default) or seekpath)
    -cont:  continuity of the k-path (default is True)
    -p:     prefix for the .bands file containing the phonon dispersion
    -fn:    file name for the ANPHON input
    
-> in the console run 'anphon phband.in > phband.log'



Step 32:
    plot the harmonic phonon dispersion
    No arguments to set here, file names were either stored or will be found by the script
    
    
    
    
Step 4:
    create the random displacements
    calls displace() function for disp_mode = random and set magnitude
    generates the geometry block of a .d12 file for each displacement (additional commands may be read from a CRY_TEMP file within the same directory)
    generated aseAtoms() objects of the random displacements are saved via json


    Arguments:
    -mag:   magnitude of the displacement (default 0.01 Å)
    -nd:    number of random displacements that should be created (default 1)
    
-> run all GRADCAL calculations with CRYSTAL to get the forces
  
  
    
Step 5:
    extract forces from the random displacements via extract()
    writes the DFSET_random that alm can read
    create ALM2.in and ALM3.in to fit the second and third order forces

    Arguments:
    -p:     prefix for the .xml file containing the second and third order forces (two should be given)
    -fn:    file name for the ALM input (two should be given)
    
-> in the console run 'alm ALM2.in > ALM2.log'
-> on the cluster run 'alm ALM3.in > ALM3.log'
    

    
Step 51:
    write the scph.in to run the self-consistent phonon calculations

    Arguments:
    -ks:    k-point source (ase (default) or seekpath)
    -cont:  continuity of the k-path (default is True)
    -p:     prefix for the quartic .xml file
    
-> on the cluster run 'anphon scph.in > scph.log'
-> check if the scph converged for each temperature
    -> here sometimes additional steps might be needed if the self consistent phonon calculations do not converge -> e.g with CV tag set

-> if scph converges nicely this will generate the .scph_dfc2 file that is needed for the next steps
-> for each temperature that should be considered further on do the following:
    -> call dfc2 and go through the wizard (set correct xml file with harmonic forces and temperature)
    -> the name of the created new .xml needs to include the temperature for which it was set!!



Step 52:
    plot the scph fitting of the phonon modes at all temperatures
    the harmonic phonon dispersion will be added as dashed lines in the plot
    file names were either stored or will be found automatically
    
Step 53:
    plot the scph fitting for 300 K and the harmonic dispersion (added as dashed lines in the plot)
    file names were either stored or will be found automatically   
    
Step 54:
    plot scph fitting of the phonon modes at different temperatures
    the harmonic phonon dispersion will be added as dashed lines in the plot
    file names were either stored or will be found automatically

    Arguments:
    -temp:  temperature(s) for which the anharmonic phonon dispersion(s) should be plotted
    


Step 6:
    write evec(temp).in files to get the eigenvector files at different temperatures
    if run with anphon they provide the .evec files for the next step

    Arguments:
    -temp:  temperature(s) for which the calculation(s) should be done
    
-> in the console run 'anphon evec(temp).in > evec(temp).log'
    


Step 7: (improved speed compared to the original code by just reading in the eigenvector at gamma)
    create random_normalcoord_pyro displacements for each mode that would lead to polarization in the direction given by the space group (look up in table -> automatize later?)
    generates the geometry block of a .d12 file for each displacement (additional commands may be read from a CRY_TEMP_PYRO file within the same directory (basically CRY_TEMP without the GRADCAL command))
    if the -flip option is set, the sign of the eigenvector coordinates will be flipped
    if the -dc option is set the displacement direction along the (given) polarization axis is checked and displacements are corrected such that they occurr in the "more common" direction.
    
    -> TBC: within this step the displace function should (after creating the displacement) should check if the displacement vector of one specific mode points in + or - direction for each atom for the different temperatures
    -> they need to point in the same direction and be flipped for that temperature if not (QUESTION: change sign of only polarization direction or all directions -> test if there is a difference)
    -> idea on how to implement:
        get a vector mapping un-displaced atom to displaced atom for adding or substracting the displacement (which should have an equal effect to flipping the signs of the eigenvectors, if not just manipulate a second set of thos eigenvectors that have opposite signs)
        see if that vecor is positive or negative (along the polarization direction) for each temperature
        build two sets of displacements, one with positive one with negative displacement along that direction
        do the following steps for both sets (slight re-ajustment for the code to loop to two sets of files (normal and flip))
    
    for each displacement from the this step one temp_modeXX.d3 file needs to be created to calculate the polarisation. The files contain this block (hard-coded). Those can only be run, after the CRYSTAL calculation has been finished (.w file needed)
        NEWK
        4 4   #   this needs to be adjusted depending on the size of the k-mesh
        1 0
        POLARI
        END    
    
        
    Arguments:
    -temp:  temperature(s) for which the calculation(s) should be done
    -mode:  phonon mode along which polarization should happen
    -flip:  flips all eigenvectors (default: False)
    -dc:    performs a direction check and correction for the displacements (default: False)
    -po:    polarization orientation (default: z)

-> run all SCP calculations with CRYSTAL to get the wave-functions
-> run all .d3 calculations with CRYSTAL afterwards

Step 7.0: Displace along the potential energy surface (non supercell) for each specified mode. A crystallographic supercell is cut according to the crystal system. The q range can be given by -q with 'start stop step'.

    Arguments:
    -temp:
    -mode:
    -crySys:  crystallographic system (ortho, rhomb, etc.) to select the right matrix for cutting out the crystallographic cell
    -q:       q range for the PES scan according to 'start stop step'. So please give three values


Step 7.1: Displace the non-supercell primitvive cell with the eigenvector. Then make a new "supercell" that is just a crystallographic cell by cutting along the matrix that is provided by CRYSTAL in the appendix; which cutting matrix is used will be determined by the input symmetry (aka ortho, rhomb, etc. -> implement as it goes)

    Arguments:
    -temp:
    -mode:
    -crySys:  crystallographic system (ortho, rhomb, etc.) to select the right matrix for cutting out the crystallographic cell


Step 8: --> obsolete, combined with prev. step (get rid of it in the code)
    for each displacement from the previous step one temp_modeXX.d3 file needs to be created to calculate the polarisation. The files contain this block (hard-coded)
        NEWK
        4 4   #   this needs to be adjusted depending on the size of the k-mesh
        1 0
        POLARI
        END
    
    Arguments:
    -temp:  temperature(s) for which the calculation(s) should be done
    -mode:  phonon mode along which polarization should happen

-> run all .d3 calculations with CRYSTAL afterwards


    
    
Step 9:
    for each temperature and mode the following needs to be done: 
        copy the wavefuction from the -20 K to modeX_temp.w (with temp being the X00 temp)
        copy the -20K .polari to modeX_temp.polari0
        copy the +20K .polari to modeX_temp.polari1
    this is taken care of via a subprocess.run() function working only if the script is run in a console
    for each temperature make a d3 file modeX_temp.d3 with:
            SPOLBP
            END

    Arguments:
    -temp:  temperature(s) for which the calculation(s) should be done
    -mode:  phonon mode along which polarization should happen

-> run all .d3 calculations with runcrys (if they take longer just submitt them to the cluster)



Step10:
    read out the P_s values along the direction for which the polarisation is noticable
    plot those values P_s (y) against temperature (x)
    -> here it would be visible which modes at which temperature would need to be flipped since their curves in the diagram show a downward trend
    -> also some corrections using the correcting factor P might have to be applied
    -> idea about handling this:
        if bith sets of displacemnets would be made and calculated from the beginning one could compare the differences between each two points and take the one that has a higher P_s for plotting
        also checking for the possibility that the P_s is actually one from a higher/lower branch
        the other (lower) values would then be used for a second figure so that both can be compared
        QUESTION: since I'm comparing every P_s with the previous, how do I know that the direction of the displacement for e.g. the 300 K was the right one and that I don't need the flipped eigenvector for that one?
        
    Arguments:
    -temp:  temperature(s) for which the calculation(s) should be done
    -mode:  phonon mode along which polarization should happen

"""


import matplotlib.pyplot as plt



import argparse  #  to be able to set variables when calling the script
import sys
import os
import subprocess
import numpy as np


from ase.io import cif
from ase.io.crystal import write_crystal
from ase.io.jsonio import write_json, read_json
from ase.build import make_supercell
from ase.spacegroup.symmetrize import check_symmetry

from ASE_Alamode_Interface import write_ALM, write_ANPHON, displace, extract, get_polarization_direction_from_sg, add_tags_to_crystal, get_evec_prop, get_pol_dir

from cif2D12 import write_d12_GEOM_from_aseAtoms, standardize_atomsAtom_crstalStructure, correct_sorting_of_stand_to_prev_ase

from plottingFunctions import plot_phonon_bands

from cif2ALM import get_fcsXML


#  imports to use the class:
    
import json


#  -> implement better way for the eigenvector flipping (maybe set different set of mode and tempearture tuples that will be treated as pairs)
#  -> rethink some naming so that it can be more flexible with filesnames including the _flip extension

class pyroCalculation():
    
    def __init__(self, **kwargs):
            
        #  step 0
        self.cif = kwargs.get('cif', None)
        self.dimension = kwargs.get('dimension', None)
        self.standardized = kwargs.get('standardize', None)
        # self.cell = kwargs.get('cell', None)
        # self.supercell = kwargs.get('supercell', None)
        
        #  step 1
        self.prefix_pattern = kwargs.get('prefix_pattern', None)
        
        #  step 2
        self.magnitude_harm = kwargs.get('magnitude_harm', None)
        # self.harm_disp = kwargs.get('harm_disp', [])
        
        #  step 3
        self.prefix_xml_harm = kwargs.get('prefix_xml_harm', None)
        
        #  step 31
        self.prefix_phonon_bands = kwargs.get('prefix_phonon_bands', None)
        
        #  step 4
        self.magnitude_rdm = kwargs.get('magnitude_rdm', None)
        # self.rdm_disp = kwargs.get('rdm_disp', '')
        
        #  step 5
        self.prefix_xlm_cubic = kwargs.get('prefix_xlm_cubic', None)
        self.prefix_xlm_quartic = kwargs.get('prefix_xlm_quartic', None)
        
        #  step 51
        self.prefix_scph_bands = kwargs.get('prefix_scph_bands', None)
        self.prefix_dfc2 = kwargs.get('prefix_dfc2', None)
        
        #  step 6
        self.temperatures = kwargs.get('temperatures', None)
        
        #  step 7
        self.pyro_modes = kwargs.get('pyro_modes', None)
        self.pol_orient = kwargs.get('pol_orient', None)
        # self.pyro_disp = kwargs.get('pyro_disp', [])
        
        self.space_group = kwargs.get('space_group', None)
        
        
        
        
    def save_pyroCalculation(self, pyro_Calc_file):
        
        with open(pyro_Calc_file, 'w') as out_file:
            json.dump(vars(self), out_file)
        
        out_file.close()
        
        
        return
            
def chem(formula, charge=''):
    sumForm = ''
    numFound = False
    
    for i, l in enumerate(formula):
        
        #print('Data-Type of', l, ':', type(l))
                
        if l.isdigit() == False and numFound == False:
            sumForm = sumForm + l
        
        if l.isdigit() == True and numFound == True:
            sumForm = sumForm + l
        
        if l.isdigit() == True and numFound == False:
            sumForm = sumForm + '$_{' + l
            numFound = True
    
        if l.isdigit() == False and numFound == True:
            sumForm = sumForm + '}$' + l
            numFound = False          
    
    if numFound == True:
        sumForm = sumForm + '}$'
    
    if charge != '':
        sumForm = sumForm + '$^{' + charge + '}$'
    
    return sumForm




def eos_raccoon():
    
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

    
    return


if __name__ == '__main__':
    
    
    print(r'/\___/\   ')
    print(r':O: :O:   Script for the pyroelectricity workflow')
    print(r' ¨   ¨')
    print()
    
    

    parser = argparse.ArgumentParser(description='Script to walk through the pyroelectricity wokflow step-by-step.', formatter_class=argparse.RawTextHelpFormatter, epilog='This script will require additional arguments depending on the step that is selected. See code for a more detailed description of each step.')
    
    
    parser.add_argument("-s", "--STEP", help='Set the step of the calculation. See code for more details on each step.\n', type=str, nargs='?', default='')
    
    parser.add_argument("-opt", "--OPTION", help='Set the sub-step option of the calculation. See code for more details on each step.\n', type=str, nargs='?', default='')

    parser.add_argument("-p", "--PREFIX", nargs='+', type=str, help="Set the PREFIX attribute of the ALM or ANPHON &general input block.", default=['', '', ''])  #  sets the ALM['PREFIX] attribute
    
    parser.add_argument("-cif", "--CIF", help="CIF-file from which the structure is read; If not set cell and supercell are loaded if available.", nargs='?', default='', type=str)

    parser.add_argument("-d", "-dim", "--DIM", help="Set the dimension of the supercell; If not set a 3x3x3 supercell is created.", nargs='+', type=int, default=[])
    
    parser.add_argument("-std", "--STANDARDIZE", help='Set to prim or crys to standardize the unit cell by spg-lib to the primitive or crystallographic cell', nargs='?', type=str, default='')
    
    parser.add_argument("-fn", "--FNAME", help="Set the filename for the created ALM/ANPHON input(s).", nargs='+', default=['', '', ''], type=str)

    parser.add_argument("-mag", "--MAGNITUDE", help="Set the magnitude of the displacement in angström; default is 0.02 Å.", nargs='?', type=float, default=0.02)
    
    parser.add_argument("-nd", "--NO_DISP", help="Set the number of random displacements that will be created.", nargs='?', type=int, default=40)
    
    parser.add_argument("-temp", "--TEMPERATURE", help="Set (a list of) temperature(s).", nargs='+', default=[], type=int)
    
    parser.add_argument("-mode", "--PYRO_MODE", help="Set (a list of) frequencies along which the pyro displacement(s) are created.", nargs='+', default=[], type=int)

    parser.add_argument("-flip", "--FLIP", help="Set to 'True' if the sign of the eigenvectors should be flipped.", nargs='?', default=False, type=bool)

    parser.add_argument("-crySys", "--CRYSYS", help="The crystal class to select the right transformation matrix.", nargs='?', default='tetra', type=str)
    
    parser.add_argument("-q", "--Q_VALUES", help="Set the q_list for which the PES scan is done as 'start stop step' ", nargs='+', default=[-0.5, 0.5, 0.1], type=float)


    parser.add_argument("-dc", "--DIRECTION_CHECK", help="Set to 'True' if the direction of the displacements should be checked and corrected for the more common direction.", nargs='?', default=False, type=bool)

    parser.add_argument("-po", "--POLARIZATION_ORIENTATION", help="Set the main axis of polarization.", nargs='?', type=str, default=None)

      
    parser.add_argument("-ks", "--KPOINT_SOURCE", help="Set the source of the k-points (default ase)", nargs='?', type=str, default='ase')
    
    parser.add_argument("-cont", "--CONTINOUS", help="Set this to True so that the retrived k-path is evaluated as if it was continous.", nargs='?', type=bool, default=True)
            
    parser.add_argument("-ovwr", "--OVERWRITE", help="Set this to false, if variables of the calculation should not be overwritten after each step.", nargs='?', type=bool, default=True)



    args, unknown = parser.parse_known_args()
    cwd = os.getcwd()

    try:
        with open('pyroCalc_progress', 'r') as in_file:
            loaded_calc = json.load(in_file)
        
        in_file.close()
                
        calculation = pyroCalculation(**loaded_calc)
        
    except:
        calculation = pyroCalculation()
    
    
    # if calculation.temperatures != None and args.TEMPERATURE != []:        
        
    #     if all(a in calculation.temperatures for a in args.TEMPERATURE) == False:
    #         print('Warning: Temperature(s) of this step and the last step do not align. New temperature(s) were added to the calculation process documentation.\n')
    #         for temp in args.TEMPERATURE:
    #             if temp not in calculation.temperatures:
    #                 calculation.temperatures.append(temp)
    
                    
    # if calculation.pyro_modes != None and args.PYRO_MODE != []:        
        
    #     if all(a in calculation.pyro_modes for a in args.PYRO_MODE) == False:
    #         print('Warning: pyro_mode(s) of this step and the last step do not align. New pyro_mode(s) were added to the calculation process documantation.\n')
    #         for mode in args.PYRO_MODE:
    #             if mode not in calculation.pyro_modes:
    #                 calculation.pyro_modes.append(mode)
 
            
    # MAking sure the right variables are present
    if args.TEMPERATURE == [] and args.OPTION.lower() in ['pes', 'pyro']:
        print('')
        
        if calculation.temperatures != None:
            args.TEMPERATURE = calculation.temperatures
        else:
            print('Please specify the temperature(s) for which the step should be calculated!')
            sys.exit()
            
    if args.PYRO_MODE == [] and args.OPTION.lower() in ['pes', 'pyro']:
        if calculation.pyro_modes != None:
            args.PYRO_MODE = calculation.pyro_modes
        else:
            print('Please specify the pyro mode(s) for which the step should be calculated!')
            sys.exit()
             
                
    
    """
    Step BUILD: build cell, supercell from cif and generate ALM0.in
    
    Options:
        -std / standardize: standardizes crystal structure to either primitive or crystallographic cell (default: '', options: prim/crys)
        
        -dim / dimension: set dimension of the supercell; 3 or 9 numbers need to be given (default: [3 3 3])
        
        -fn / filename: filenalme for ALM0.in (default: '')
        
        -pf / prefix: prefix for the .pattern_file (default: '')
    
    """
    
    if args.STEP.lower() == 'build':
        
        # failsafe, to check that one does not redefine the cell accidentally
        
        try:
            with open('cell', 'r') as in_cell:
                load_cell = json.load(in_cell)
            
            in_cell.close()
            
            print('For this calculation cell and supercell have already been defined. Do you want to build a new cell/supercell?')
            
            answer = ''
            
            while answer.lower() not in ['yes', 'y', 'no', 'n']:
                answer = input('yes/no:  ')
                
            if answer.lower() in ['yes', 'y']:
                pass
            else:
                exit()
            
        except FileNotFoundError:
            pass
        
        
        # read cif 
        # standardize to primitive or crystallographic
        
        if args.CIF != '':
            cifIn = args.CIF
            compound_in = cif.read_cif(cifIn)    
            
            if args.STANDARDIZE == 'prim':
                
                print('Standardizing the crystal structure to the primitive cell.\n')
                
                std_compound = standardize_atomsAtom_crstalStructure(compound_in, to_primitive=True)
                
                compound = correct_sorting_of_stand_to_prev_ase(compound_in, std_compound)
                
                calculation.standardized = args.STANDARDIZE
                
            elif args.STANDARDIZE == 'crys':
                
                print('Standardizing the crystal structure to the crystallographic cell.\n')
                
                std_compound = standardize_atomsAtom_crstalStructure(compound_in, to_primitive=False)
                
                compound = correct_sorting_of_stand_to_prev_ase(compound_in, std_compound)
                
                calculation.standardized = args.STANDARDIZE
                
            else:
                
                print('Crystal structure used without standardizing.\n')
                
                compound = compound_in
                
                calculation.standardized = 'no'
                
                
            write_json(open('cell', 'w'), compound)

            calculation.cif = args.CIF
            
        else:
            try:
                compound = read_json(open('cell', 'r'))
                print('Cell loaded from previous step.\n')
            
            except FileNotFoundError:
                print('No saved cell and/or supercell found.\n')
                print('Please rerun script with -cif option set!\n')
                sys.exit()
        
        
        
        # build supercell from cif 
        
            #use dimension variable to create matrix for the supercell, otherwise it is set as a x,y,z = 3 3 3 supercell
        if args.DIM != None and len(args.DIM) == 3:
            matrix = [[args.DIM[0], 0, 0], [0, args.DIM[1], 0], [0, 0, args.DIM[2]]]
            print('Supercell dimension set to custom %s.\n'%args.DIM)
            
            calculation.dimension = args.DIM
        
        elif args.DIM != None and len(args.DIM) == 9:
            matrix = [[args.DIM[0], args.DIM[1], args.DIM[2]], [args.DIM[3], args.DIM[4], args.DIM[5]], [args.DIM[6], args.DIM[7], args.DIM[8]]]
            print('Supercell dimension set to custom %s.\n'%args.DIM)
            
            calculation.dimension = args.DIM
        
        elif args.DIM != None and len(args.DIM) not in [3,9]:
            print('Supercell dimension has to be of length 3 or 9.\n')

        else:
            matrix = [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
            print('Supercell dimension set by default as [3,3,3].\n')
            
            calculation.dimension = [3, 3, 3]

        
        #  create supercell
        compound_sc = make_supercell(compound, matrix, order='atom-major')
        write_json(open('supercell', 'w'), compound_sc)
        
        # make ALM0 input for harmonic pattern
        
        print('Writing input for Alamode to generate the .pattern_HARMONIC file.\n')
        
        if args.FNAME[0] == '':
            args.FNAME[0] = 'ALM0.in'
        
        if args.PREFIX[0] == '':
            args.PREFIX[0] = str(compound.get_chemical_formula(empirical=True))

        calculation.prefix_pattern = args.PREFIX[0] + '.pattern_HARMONIC'

        ALM0_in = {'PREFIX': [args.PREFIX[0]], 'MODE': ['suggest'], 'NORDER': [1]}

        write_ALM(compound_sc, file_name=args.FNAME[0], command_list=ALM0_in, get_XML_files=[False, False, False])

        print('Please run %s with alm.\n'%args.FNAME[0])


    else:  #  for all following steps it is assumed that the cell and supercell json files are available
        compound = read_json(open('cell', 'r'))
        compound_sc = read_json(open('supercell', 'r'))
        print('Cell and supercell loaded from previous step.\n')

        
        
    """
    STEP displace: displce the cell or supercell with different options:
        
        harmoic: harmonic phonon displacement
        
        anharmonic: random phonon displacement
        
        pes: along the potential energy surface (normal cell by default)
        
        pyro: along different pyro modes and different temperatures
    """
        
    if args.STEP.lower() in ['displace', 'disp']:
        
        if args.OPTION.lower() in ['harmonic', 'harm']:
            
            print('Writing harmonic displacement inputs.\n')
            
            calculation.magnitude_harm = args.MAGNITUDE
            
            
            found_harm_pattern = False
            
            if calculation.prefix_pattern != None:
                
                harm_disp = displace(compound_sc, 'pattern_file', file_pattern=calculation.prefix_pattern, mag=args.MAGNITUDE)
            
            else:

                for file in os.listdir(cwd):
                    if file.endswith('.pattern_HARMONIC'):
                        harm_pattern_file = file
                        found_harm_pattern = True
                        print('Loading harmonic displacement pattern from %s.\n'%harm_pattern_file)
                        break
                
                if found_harm_pattern == False:
                    print('No file with the harmonic displacement pattern found.\n')
                    sys.exit()
                    
                harm_disp = displace(compound_sc, 'pattern_file', file_pattern=harm_pattern_file, mag=args.MAGNITUDE)
            
            write_json(open('harmonic_displacements', 'w'), harm_disp)
            
            print('Number of created haromic displacements: %s\n'%len(harm_disp))
            
            counter = 0
            
            #  Look for a TEMPLATE file to add to the CRYSTAL input
            found_cry_temp = False
            for file in os.listdir(cwd):
                if file == 'CRY_TEMP' or file == 'TEMP':
                    found_cry_temp = True
                    crystal_template = file
                
            if found_cry_temp == False:
                crystal_template = None
                    
            for disp in harm_disp:
                counter += 1
                file_name_d12 = 'h_disp%s.d12'%counter

                add_tags_to_crystal(disp)
                    
                write_d12_GEOM_from_aseAtoms(disp, d12_filename=file_name_d12, P1=True, verbosity=0, template=crystal_template, numDig=8, external=True)
                    
                #write_crystal('h_disp%s.ext'%counter, disp)


        elif args.OPTION.lower() in ['anharmonic', 'anharm']:
            
            print('Writing random displacement inputs.\n')
            
            calculation.magnitude_rdm = args.MAGNITUDE
            
            try:
                rdm_disp = read_json(open('random_displacements', 'r'))
                print('Random displacements have already been creted for this compound in this folder!\n\nWARNING: Once overwritten, random displacements cannot be retrived!\n\nDo you still want to overwrite the random displacements?\n')
                answer = ''
                
                while answer.lower() not in ['yes', 'y', 'no', 'n']:
                    answer = input('yes/no:  ')
                    
                if answer.lower() in ['yes', 'y']:
                    print('\nNew displacements will replace existing list.\n')
                    
                    rdm_disp = displace(compound_sc, 'random', mag=args.MAGNITUDE, ndata=args.NO_DISP)
                    
                    counter = 0
                    
                    with open('calc_documentation.txt', 'a') as doc_file:
                        
                        doc_file.write('WARNING: Existing random displacements have been replaced!\n\n')
                        
                        doc_file.close()
                    
                    
                elif answer.lower() in ['no', 'n']:
                    print('\nDisplacements will be added to the existing list.\n')
                    
                    counter = len(rdm_disp)
                    
                    add_rdm_disp = displace(compound_sc, 'random', mag=args.MAGNITUDE, ndata=args.NO_DISP)
                    
                    for disp in add_rdm_disp:
                        rdm_disp.append(disp)

            except:
                counter = 0
                
                rdm_disp = displace(compound_sc, 'random', mag=args.MAGNITUDE, ndata=args.NO_DISP)
                
            
            # rdm_disp = displace(compound_sc, 'random', mag=args.MAGNITUDE, ndata=args.NO_DISP)
            
            write_json(open('random_displacements', 'w'), rdm_disp)

            
            print('Number of created random displacements: %s\n'%len(rdm_disp))
            
            
            #  Look for a TEMPLATE file to add to the CRYSTAL input
            found_cry_temp = False
            for file in os.listdir(cwd):
                if file == 'CRY_TEMP' or file == 'TEMP':
                    found_cry_temp = True
                    crystal_template = file
                    break
                
            if found_cry_temp == False:
                crystal_template = None
                    
            for disp in rdm_disp[counter::]:
                counter += 1
                file_name_d12 = 'rdm_disp%s.d12'%counter
                

                add_tags_to_crystal(disp)
                    
                write_d12_GEOM_from_aseAtoms(disp, d12_filename=file_name_d12, P1=True, verbosity=0, template=crystal_template, numDig=8, external=True)
                
                #write_crystal('rdm_disp%s.ext'%counter, disp)
            
        
            """ This requires new definition of the cell / supercell, maybe this can be handled better with the given key words..."""
        elif args.OPTION.lower() == 'pes':
            ########################################################################
            
            print('Create the crystallographic cell pes displacements and their .d12 files for:')
            
            q_list = [str(np.round(q, 5)) for q in np.array(np.arange(args.Q_VALUES[0], args.Q_VALUES[1]+args.Q_VALUES[2], args.Q_VALUES[2]))]
            
            
            print('q_list: %s'%q_list)
            print('Modes: %s\n'%args.PYRO_MODE)
            
            
            if calculation.pyro_modes == None:
                calculation.pyro_modes = args.PYRO_MODE
            
            ndata = int(((args.Q_VALUES[1] - args.Q_VALUES[0]) / args.Q_VALUES[2]) + 1)
            
            print(ndata)
            print()
            
            evec_files = []
            pes_disp = [[] for i in args.PYRO_MODE]
            
            # if args.FLIP == True:
            #     flip_ext = '_flip'
            # else:
            #     flip_ext = ''


            temp = args.TEMPERATURE[0]
            temp = args.TEMPERATURE[0]

            for file in os.listdir(cwd):
                if file.find('_%s'%temp)!=-1 and file.find('.evec') !=-1:
                    evec_file = file
                    break
                
                
            for file in os.listdir(cwd):
                if file == 'CRY_TEMP_PYRO':
                    crystal_template = file
                    break
                else:
                    crystal_template = None
            
            
            for modeNum, mode in enumerate(args.PYRO_MODE):

                cur_disp_pes = displace(compound, 'pes', evec_file = evec_file, temp=temp, pyro_mode=mode, primitive_aseAtom=compound, qrange=[args.Q_VALUES[0], args.Q_VALUES[1]], ndata=ndata)
                
                    
                if args.CRYSYS == 'orthoA':
                    print('Crystal system orthorhombic A centered selected.\n')
                        
                    matrix = [[1,0,0], [0,1,-1], [0,1,1]] #ortho A centered
                        
                elif args.CRYSYS == 'orthoC':
                    print('Crystal system orthorhombic C centered selected.\n')
                        
                    matrix = [[1,1,0], [-1,1,0], [0,0,1]] #ortho C centered
                    
                elif args.CRYSYS == 'rhomb':
                    print('Crystal system rhombohedral selected.\n')
                    
                    matrix = [[1,-1,0], [0,1,-1], [1,1,1]]
                    
                elif args.CRYSYS == 'tetra':
                    
                    matrix = [[1,0,0], [0,1,0], [0,0,1]]
                    
                else:
                    print('The crystal class is not recognized by the script, please add the right matrix for it.', '\n')
                    matrix = [[1,0,0], [0,1,0], [0,0,1]]
                
                
                for dispNum, disp in enumerate(cur_disp_pes):
                    
                    add_tags_to_crystal(disp)
                    
                    pes_sc = make_supercell(disp, matrix, order='atom-major')
                    
                    if args.DIM != None and len(args.DIM) == 3:
                        matrixSC = [[args.DIM[0], 0, 0], [0, args.DIM[1], 0], [0, 0, args.DIM[2]]]
                        
                        pes_sc = make_supercell(pes_sc, matrixSC, order='atom-major')
                        
                    elif args.DIM != None and len(args.DIM) == 9:
                        matrixSC = [[args.DIM[0], args.DIM[1], args.DIM[2]], [args.DIM[3], args.DIM[4], args.DIM[5]], [args.DIM[6], args.DIM[7], args.DIM[8]]]
                        
                        pes_sc = make_supercell(pes_sc, matrixSC, order='atom-major')
                    
                    pes_disp[modeNum].append(pes_sc)
                    
                    cry_file_name = f'pes_disp_{q_list[dispNum]}_mode{mode}.d12'
                    prop_file_name = f'pes_disp_{q_list[dispNum]}_mode{mode}.d3'
                
                
                    write_d12_GEOM_from_aseAtoms(pes_sc, d12_filename=cry_file_name, P1=True, verbosity=0, template=crystal_template, numDig=8, external=True)
                
                    #write_crystal(f'{temp}_mode{mode}{flip_ext}.ext', pyro_sc)
                    
                    with open(prop_file_name, 'w') as polari_file:
                        polari_file.write('NEWK\n')
                        polari_file.write('4 4\n')
                        polari_file.write('1 0\n')
                        polari_file.write('POLARI\n')
                        polari_file.write('END')
                        polari_file.close()
                    
                    print(f'Displacement for mode {mode} at q_value {q_list[dispNum]} created.')
                print()
            
            print()
            write_json(open('pes_displacements', 'w'), pes_disp)

            
            
        elif args.OPTION.lower() == 'pyro':
            
            print('Create the crystallographic cell pyro displacements and their .d12 files for:')
            print('Temperatures: %s'%args.TEMPERATURE)
            print('Modes: %s\n'%args.PYRO_MODE)
            
            if calculation.pyro_modes == None:
                calculation.pyro_modes = args.PYRO_MODE
            
            temps = []
            evec_files = []
            pyro_disp = [[] for i in args.PYRO_MODE]
            
            # print(pyro_disp)
            
            
            # Flipping should be handled automaticall as well as making sure the mode stays the same over the temperature range
            
            # Idea: mode is identified at the lowest temperature and then all others are aligned to that one and the right number is found for other temperatures -> check out that one method that was written for the plotting evec thing...
            if args.FLIP == True:
                flip_ext = '_flip'
            else:
                flip_ext = ''
            
            #  create matrix with -20/+20 temperatures        
            for temp in args.TEMPERATURE:
                temps.append(temp-20)
                temps.append(temp+20)
                
                for file in os.listdir(cwd):
                    if file.find('_%s'%temp)!=-1 and file.find('.evec') !=-1:
                        evec_file = file
                        break
                
                if len(evec_files) > 0 and evec_file == evec_files[-1]:
                    raise FileNotFoundError(f'The evec file for {temp} was not found!')
                
                evec_files.append(evec_file)
                evec_files.append(evec_file)
            
            
            for file in os.listdir(cwd):
                if file == 'CRY_TEMP_PYRO':
                    crystal_template = file
                    break
                else:
                    crystal_template = None

            
            #  This step also needs to loop trough modes for each of the +/-temperatures
            for modeNum, mode in enumerate(args.PYRO_MODE):        
                for tempNum, temp in enumerate(temps):
                    
                    cur_pyro_disp = displace(compound, 'random_normalcoord_pyro', evec_file = evec_files[tempNum], temp=temp, pyro_mode=mode, primitive_aseAtom=compound, flip=args.FLIP)
                    
                    if args.CRYSYS == 'orthoA':
                        print('Crystal system orthorhombic A centered selected.\n')
                        
                        matrix = [[1,0,0], [0,1,-1], [0,1,1]] #ortho A centered
                        
                    elif args.CRYSYS == 'orthoC':
                        print('Crystal system orthorhombic C centered selected.\n')
                        
                        matrix = [[1,1,0], [-1,1,0], [0,0,1]] #ortho C centered
                    
                    elif args.CRYSYS == 'rhomb':
                        print('Crystal system rhombohedral selected.\n')
                    
                        matrix = [[1,-1,0], [0,1,-1], [1,1,1]]
                    
                    elif args.CRYSYS == 'tetra':
                    
                        matrix = [[1,0,0], [0,1,0], [0,0,1]]
                    
                    else:
                        print('The crystal class is not recognized by the script, please add the right matrix for it.', '\n')
                        matrix = [[1,0,0], [0,1,0], [0,0,1]]
                        
                    
                    pyro_sc = make_supercell(cur_pyro_disp[0], matrix, order='atom-major')
                    
                    
                    if args.DIM != None and len(args.DIM) == 3:
                        matrixSC = [[args.DIM[0], 0, 0], [0, args.DIM[1], 0], [0, 0, args.DIM[2]]]
                        
                        pyro_sc = make_supercell(pyro_sc, matrixSC, order='atom-major')
                        
                    elif args.DIM != None and len(args.DIM) == 9:
                        matrixSC = [[args.DIM[0], args.DIM[1], args.DIM[2]], [args.DIM[3], args.DIM[4], args.DIM[5]], [args.DIM[6], args.DIM[7], args.DIM[8]]]
                        
                        pyro_sc = make_supercell(pyro_sc, matrixSC, order='atom-major')

                    
                    pyro_disp[modeNum].append(pyro_sc)
                    
                    cry_file_name = f'{temp}_mode{mode}{flip_ext}.d12'
                    prop_file_name = f'{temp}_mode{mode}{flip_ext}.d3'
                    
                    add_tags_to_crystal(pyro_sc)
                    
                    write_d12_GEOM_from_aseAtoms(pyro_sc, d12_filename=cry_file_name, P1=True, verbosity=0, template=crystal_template, numDig=8, external=True)
                
                    #write_crystal(f'{temp}_mode{mode}{flip_ext}.ext', pyro_sc)
                    
                    with open(prop_file_name, 'w') as polari_file:
                        polari_file.write('NEWK\n')
                        polari_file.write('4 4\n')
                        polari_file.write('1 0\n')
                        polari_file.write('POLARI\n')
                        polari_file.write('END')
                        polari_file.close()
                    
                    print(f'Displacement for mode {mode} at {temp} K created.')
            
            print()
            write_json(open('pyro_displacements', 'w'), pyro_disp)

        
        elif args.OPTION.lower() in ['pyroa', 'pyro auto', 'pyro automatic']:
            
            print('Create the crystallographic cell pyro displacements and their .d12 files atomatically for:')
            print('Temperatures: %s'%args.TEMPERATURE)
            
            get_evec_prop(compound, args.TEMPERATURE)
            
            pyroDispPattern = get_pol_dir(compound, args.TEMPERATURE)
            
            autoPyroDisplacements = []
            
            for pyroDisp in pyroDispPattern:
                
                temps = [pyroDisp[0]-20, pyroDisp[0]+20]
                
                for file in os.listdir(cwd):
                    if file.find('_%s'%pyroDisp[0])!=-1 and file.find('.evec') !=-1:
                        evec_file = file
                        break
                
                
                for file in os.listdir(cwd):
                    if file == 'CRY_TEMP_PYRO':
                        crystal_template = file
                        break
                    else:
                        crystal_template = None


                for temp in temps:
                    
                    cur_pyro_disp = displace(compound, 'random_normalcoord_pyro', evec_file = evec_file, temp=temp, pyro_mode=pyroDisp[1], primitive_aseAtom=compound, flip=pyroDisp[2])
                    
                    if args.CRYSYS == 'orthoA':
                        print('Crystal system orthorhombic A centered selected.\n')
                        
                        matrix = [[1,0,0], [0,1,-1], [0,1,1]] #ortho A centered
                        
                    elif args.CRYSYS == 'orthoC':
                        print('Crystal system orthorhombic C centered selected.\n')
                        
                        matrix = [[1,1,0], [-1,1,0], [0,0,1]] #ortho C centered
                    
                    elif args.CRYSYS == 'rhomb':
                        print('Crystal system rhombohedral selected.\n')
                    
                        matrix = [[1,-1,0], [0,1,-1], [1,1,1]]
                    
                    elif args.CRYSYS == 'tetra':
                    
                        matrix = [[1,0,0], [0,1,0], [0,0,1]]
                    
                    else:
                        print('The crystal class is not recognized by the script, please add the right matrix for it.', '\n')
                        matrix = [[1,0,0], [0,1,0], [0,0,1]]
                        
                    
                    pyro_sc = make_supercell(cur_pyro_disp[0], matrix, order='atom-major')
                    
                    
                    if args.DIM != None and len(args.DIM) == 3:
                        matrixSC = [[args.DIM[0], 0, 0], [0, args.DIM[1], 0], [0, 0, args.DIM[2]]]
                        
                        pyro_sc = make_supercell(pyro_sc, matrixSC, order='atom-major')
                        
                    elif args.DIM != None and len(args.DIM) == 9:
                        matrixSC = [[args.DIM[0], args.DIM[1], args.DIM[2]], [args.DIM[3], args.DIM[4], args.DIM[5]], [args.DIM[6], args.DIM[7], args.DIM[8]]]
                        
                        pyro_sc = make_supercell(pyro_sc, matrixSC, order='atom-major')

                    
                    autoPyroDisplacements.append(pyro_sc)
                    
                    cry_file_name = f'{temp}_mode{pyroDisp[1]}.d12'
                    prop_file_name = f'{temp}_mode{pyroDisp[1]}.d3'
                    
                    add_tags_to_crystal(pyro_sc)
                    
                    write_d12_GEOM_from_aseAtoms(pyro_sc, d12_filename=cry_file_name, P1=True, verbosity=0, template=crystal_template, numDig=8, external=True)
                
                    #write_crystal(f'{temp}_mode{mode}{flip_ext}.ext', pyro_sc)
                    
                    with open(prop_file_name, 'w') as polari_file:
                        polari_file.write('NEWK\n')
                        polari_file.write('4 4\n')
                        polari_file.write('1 0\n')
                        polari_file.write('POLARI\n')
                        polari_file.write('END')
                        polari_file.close()
                    
                    print(f'Displacement for mode {pyroDisp[1]} at {temp} K created.')
            
            print()
            write_json(open('pyro_displacements_auto', 'w'), autoPyroDisplacements)

            
            

            
        else:
            
            print('Please choose one of the following options: harmonic, anharmonic, pes or pyro, pyroa!\n')
        
        
    
        """
    STEP extract: extract forces from harmonic/anharmonic force calculations and write DFSET files for the options:
        
        harmonic: also creates ALM1.in and phband.in
        
        anharmonic: creates ALM2.in, ALM3.in and scph.in
    
    
    for pyro and PES scan prepare the files for Barry phase calculations such that the last step can be run immediately:
        
        copies and renames the right .polari and .w functions, writes inputs for Berry phase calculation
    """
    elif args.STEP.lower() in ['extract', 'extr']:
        
        if args.OPTION.lower() in ['harmonic', 'harm']:
            
            print('Extracting forces and epot of harmonic displacements.\n')
            
            try:
                harm_disp = read_json(open('harmonic_displacements', 'r'))
                print('Loading harmonic displacements from previous step.\n')
            except:
                raise ValueError('Failed to load harmonic displacements. Please check if step 2 was executed properly!')
            
            # harm_disp_out_files = []
            
            # for file in os.listdir(cwd):
            #     if file.find('h_disp')!=-1 and file.find('.out') !=-1:
            #         harm_disp_out_files.append(file)
            
            # print(len(harm_disp))
            # print()
            
            harm_disp_out_files = [f'h_disp{aNum+1}.out' for aNum, a in enumerate(harm_disp)]
            
            # print(harm_disp_out_files)
            # print()

            extract(compound_sc, harm_disp, file_name='DFSET_harmonic', in_unit_energy='hartree', CRYSTAL=True, out_files = harm_disp_out_files)

            print('Crating ALM1.in input for Alamode to fit the harmonic forces.\n')

            if args.FNAME[0] == '':
                args.FNAME[0] = 'ALM1.in'

            if args.PREFIX[0] == '':
                args.PREFIX[0] = str(compound.get_chemical_formula(empirical=True))
                
            calculation.prefix_xml_harm = args.PREFIX[0] + '.xml'

            ALM1_in = {'PREFIX': [args.PREFIX[0]], 'MODE': ['opt'], 'NORDER': [1], 'DFSET': ['DFSET_harmonic']}

            write_ALM(compound_sc, file_name=args.FNAME[0], command_list=ALM1_in, get_XML_files=[False, False, False])


            print('Creating phband.in file to run with anphon.\n')
            
            if args.FNAME[1] == '':
                args.FNAME[1] = 'phband.in'
    
            if args.PREFIX[1] == '':
                args.PREFIX[1] = str(compound.get_chemical_formula(empirical=True)) + '_bands'
                    
            calculation.prefix_phonon_bands = args.PREFIX[1] + '.bands'
            
            
            if calculation.prefix_xml_harm != None:
                phband_in = {'PREFIX': [args.PREFIX[1]],'MODE': ['phonons'], 'FCSXML': [calculation.prefix_xml_harm]}
                            
                if args.STANDARDIZE == 'prim':
                    write_ANPHON(compound, file_name=args.FNAME[1], command_list=phband_in, get_XML_files=False, contPath=args.CONTINOUS, source=args.KPOINT_SOURCE)
                    
                elif args.STANDARDIZE == 'crys':
                    write_ANPHON(compound, file_name=args.FNAME[1], command_list=phband_in, get_XML_files=False, contPath=args.CONTINOUS, source=args.KPOINT_SOURCE, stdIn = 'crys')

    
            else:
    
                phband_in = {'PREFIX': [args.PREFIX[1]],'MODE': ['phonons']}
                
                if args.STANDARDIZE == 'prim':

                    write_ANPHON(compound, file_name=args.FNAME[1], command_list=phband_in, get_XML_files=True, contPath=args.CONTINOUS, source=args.KPOINT_SOURCE)
                    
                elif args.STANDARDIZE == 'crys':
                    write_ANPHON(compound, file_name=args.FNAME[1], command_list=phband_in, get_XML_files=True, contPath=args.CONTINOUS, source=args.KPOINT_SOURCE, stdIn = 'crys')
                    


            
            
        elif args.OPTION.lower() in ['anharmonic', 'anharm']:
            
            print('Extracting forces and epot of random displacements.\n')
            
            try:
                rdm_disp = read_json(open('random_displacements', 'r'))
                print('Loading random displacements from previous step.\n')
            except:
                raise ValueError('Failed to load random displacements. Please check if step 4 was executed properly!')
            
            # rdm_disp_out_files = []
            
            # for file in os.listdir(cwd):
            #     if file.find('rdm_disp')!=-1 and file.find('.out') !=-1:
            #         rdm_disp_out_files.append(file)
            
            rdm_disp_out_files = [f'rdm_disp{aNum+1}.out' for aNum, a in enumerate(rdm_disp)]

            extract(compound_sc, rdm_disp, file_name='DFSET_random', in_unit_energy='hartree', CRYSTAL=True, out_files = rdm_disp_out_files)

            print('Crating ALM2.in and ALM3.in input for Alamode to optimize the random forces with second and third order.\n')  #  Is that really what is happening?? -> ask Kim

            if args.FNAME[0] == '':
                args.FNAME[0] = 'ALM2.in'
                
            if args.FNAME[1] == '':
                args.FNAME[1] = 'ALM3.in'

            if args.PREFIX[0] == '':
                args.PREFIX[0] = str(compound.get_chemical_formula(empirical=True)) + '_cubic'
                
            if args.PREFIX[1] == '':
                args.PREFIX[1] = str(compound.get_chemical_formula(empirical=True)) + '_quartic'

            calculation.prefix_xlm_cubic = args.PREFIX[0] + '.xml'
            calculation.prefix_xlm_quartic = args.PREFIX[1] + '.xml'


            if calculation.prefix_xml_harm != None:
                ALM2_in = {'PREFIX': [args.PREFIX[0]], 'MODE': ['opt'], 'NORDER': [2], 'DFSET': ['DFSET_random'], 'FC2XML': [calculation.prefix_xml_harm]}
                
                ALM3_in = {'PREFIX': [args.PREFIX[1]], 'MODE': ['opt'], 'NORDER': [3], 'DFSET': ['DFSET_random'], 'FC2XML': [calculation.prefix_xml_harm], 'FC3XML': [calculation.prefix_xlm_cubic]}
                
                write_ALM(compound_sc, file_name=args.FNAME[0], command_list=ALM2_in, get_XML_files=[False, False, False])
                
                write_ALM(compound_sc, file_name=args.FNAME[1], command_list=ALM3_in, get_XML_files=[False, False, False])
               
                
            else:
                ALM2_in = {'PREFIX': [args.PREFIX[0]], 'MODE': ['opt'], 'NORDER': [2], 'DFSET': ['DFSET_random']}
                
                ALM3_in = {'PREFIX': [args.PREFIX[1]], 'MODE': ['opt'], 'NORDER': [3], 'DFSET': ['DFSET_random'], 'FC3XML': [args.PREFIX[0] + '.xml']}

                write_ALM(compound_sc, file_name=args.FNAME[0], command_list=ALM2_in, get_XML_files=[True, False, False])
                
                write_ALM(compound_sc, file_name=args.FNAME[1], command_list=ALM3_in, get_XML_files=[True, False, False])

            
            print('Creating scph.in file to run with anphon.\n')
            
            if args.FNAME[2] == '':
                args.FNAME[2] = 'scph.in'

            if args.PREFIX[2] == '':
                args.PREFIX[2] = str(compound.get_chemical_formula(empirical=True)) + '_scph'

                
            calculation.prefix_scph_bands = args.PREFIX[2] + '.scph_bands'
            calculation.prefix_dfc2 = args.PREFIX[2] + '.scph_dfc2'
            
            
            if calculation.prefix_xlm_quartic != None:
                scph_in = {'PREFIX': [args.PREFIX[2]],'MODE': ['SCPH'], 'SELF_OFFDIAG': [0], 'MAXITER': [1000], 'MIXALPHA': [0.1], 'KMESH_INTERPOLATE': [4, 4, 4], 'KMESH_SCPH': [8, 8, 8], 'LOWER_TEMP': 0, 'DT': [50], 'FCSXML': [calculation.prefix_xlm_quartic]}
                
                if args.STANDARDIZE == 'prim':
                
                    write_ANPHON(compound, file_name=args.FNAME[2], command_list=scph_in, get_XML_files=False, contPath=args.CONTINOUS, source=args.KPOINT_SOURCE)
                    
                elif args.STANDARDIZE == 'crys':
                    
                    write_ANPHON(compound, file_name=args.FNAME[2], command_list=scph_in, get_XML_files=False, contPath=args.CONTINOUS, source=args.KPOINT_SOURCE, stdIn = 'crys')
                    
                
            else:
            
                scph_in = {'PREFIX': [args.PREFIX[2]],'MODE': ['SCPH'], 'SELF_OFFDIAG': [0], 'MAXITER': [1000], 'MIXALPHA': [0.1], 'KMESH_INTERPOLATE': [4, 4, 4], 'KMESH_SCPH': [8, 8, 8], 'LOWER_TEMP': 0, 'DT': [50]}
                
                if args.STANDARDIZE == 'prim':
                
                    write_ANPHON(compound, file_name=args.FNAME[2], command_list=scph_in, get_XML_files=True, contPath=args.CONTINOUS, source=args.KPOINT_SOURCE)
                    
                elif args.STANDARDIZE == 'crys':
                    
                    write_ANPHON(compound, file_name=args.FNAME[2], command_list=scph_in, get_XML_files=True, contPath=args.CONTINOUS, source=args.KPOINT_SOURCE, stdIn = 'crys')

            
            
            
            
        elif args.OPTION.lower() == 'pes':
            
            #modeList = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
            modeList = args.PYRO_MODE
            
            
    
            #q_list_eval = ['-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0.1', '0.2', '0.3', '0.4', '0.5']
            
            q_list_eval = [str(round(q, 5)) for q in np.array(np.arange(args.Q_VALUES[0], args.Q_VALUES[1]+args.Q_VALUES[2], args.Q_VALUES[2]))]

            
            for mode in modeList:
            
                for q in q_list_eval:
                
                    with open(f'pes_mode{mode}_{q}.d3', 'w') as polari_file:
                        polari_file.write('SPOLBP\n')
                        polari_file.write('END')
                        polari_file.close()
                    
                
                    subprocess.run(f'cp pes_disp_-0.0_mode{mode}.polari pes_mode{mode}_{q}.polari0', shell=True)
                # pes_disp_0_mode1.polari    pes_mode1_01.polari0
                
                    subprocess.run(f'cp pes_disp_{q}_mode{mode}.polari pes_mode{mode}_{q}.polari1', shell=True)
                # pes_disp_01_mode1.polari    pes_mode1_01.polari1
                
                    subprocess.run(f'cp pes_disp_{q}_mode{mode}.w pes_mode{mode}_{q}.w', shell=True)
                # pes_disp_01_mode1.w    pes_mode1_01.w
    
                    subprocess.run(f'runcrys pes_mode{mode}_{q}.d3', shell=True)

            
            
        elif args.OPTION.lower() == 'pyro':
            
            print('Create the SPOLB .d3 files for each temperature for each mode:')
            print('Temperatures: %s'%args.TEMPERATURE)
            print('Modes: %s\n'%args.PYRO_MODE)
            
            
            for mode in args.PYRO_MODE:
                for temp in args.TEMPERATURE:
                    
                    temp_lower = temp - 20
                    temp_higher = temp + 20
                        
                    SPOLBP_name = f'mode{mode}_{temp}.d3'
                    
                    with open (SPOLBP_name, 'w') as SPOLBP_file:
                        SPOLBP_file.write('SPOLBP\n')
                        SPOLBP_file.write('END')
                        SPOLBP_file.close()
                        
                    # this should run subprocesses of copying and renaming wavefunctions and .polari files -> works only if script is called from a console
                    #  use some method from os/sys instead of subprocess (or shutil)
                    if args.FLIP == False:
                        subprocess.run(f'cp {temp_lower}_mode{mode}.w mode{mode}_{temp}.w', shell=True)
                        subprocess.run(f'cp {temp_lower}_mode{mode}.polari mode{mode}_{temp}.polari0', shell=True)
                        subprocess.run(f'cp {temp_higher}_mode{mode}.polari mode{mode}_{temp}.polari1', shell=True)
                        subprocess.run(f'runcrys mode{mode}_{temp}.d3', shell=True)

                        
                    elif args.FLIP == True:
                        
                        subprocess.run(f'cp {temp_lower}_mode{mode}_flip.w mode{mode}_{temp}_flip.w', shell=True)
                        subprocess.run(f'cp {temp_lower}_mode{mode}_flip.polari mode{mode}_{temp}_flip.polari0', shell=True)
                        subprocess.run(f'cp {temp_higher}_mode{mode}_flip.polari mode{mode}_{temp}_flip.polari1', shell=True)
                        
                        subprocess.run(f'runcrys mode{mode}_{temp}.d3', shell=True)

            

        elif args.OPTION.lower() in ['pyroa', 'pyro auto', 'pyro automatic']:

            print('Create the SPOLB .d3 files for the modes with polarization along the x-direction automatically for:')
            print('Temperatures: %s'%args.TEMPERATURE)
            
            get_evec_prop(compound, args.TEMPERATURE)
            
            pyroDispPattern = get_pol_dir(compound, args.TEMPERATURE)
            
            
            for pyroDisp in pyroDispPattern:
                
                temp_lower = pyroDisp[0] - 20
                temp_higher = pyroDisp[0] + 20
                
                SPOLBP_name = f'mode{pyroDisp[1]}_{pyroDisp[0]}.d3'

                with open (SPOLBP_name, 'w') as SPOLBP_file:
                    SPOLBP_file.write('SPOLBP\n')
                    SPOLBP_file.write('END')
                    SPOLBP_file.close()
            
                subprocess.run(f'cp {temp_lower}_mode{pyroDisp[1]}.w mode{pyroDisp[1]}_{pyroDisp[0]}.w', shell=True)
                subprocess.run(f'cp {temp_lower}_mode{pyroDisp[1]}.polari mode{pyroDisp[1]}_{pyroDisp[0]}.polari0', shell=True)
                subprocess.run(f'cp {temp_higher}_mode{pyroDisp[1]}.polari mode{pyroDisp[1]}_{pyroDisp[0]}.polari1', shell=True)
                
                subprocess.run(f'runcrys mode{pyroDisp[1]}_{pyroDisp[0]}.d3', shell=True)
                
                
                
                
    
        """
    STEP plot: handling all plots with options for:
        
        harmonic: plot the harmonic phonon dispersion
        
        anharmonic: plots the scph phonons (different options dependin on whether temperature is given or not)
        
        pes: gives the modes plot for the pes scan
        
        pyro: plots the pyro coefficicent along different modes; if more than one temperature is given, maybe have it add the summed up line?
    """
    
    
    elif args.STEP.lower() in ['plot', 'plt']:
        
        if args.OPTION.lower() in ['harmonic', 'harm']:
            
            
            print('Plotting the harmonic phonon dispersion.\n')
            
            if args.PREFIX[0] == '':
                label = chem(str(compound.get_chemical_formula(empirical=True)))
            else:
                label = chem(args.PREFIX[0])
            
            
            if calculation.prefix_phonon_bands != None:
                plot_phonon_bands(calculation.prefix_phonon_bands, labels=[label], save_path=str(compound.get_chemical_formula(empirical=True)) + '_harm.png')
            else:
                
                for file in os.listdir(cwd):
                    if file.find('.bands') !=-1:
                        bands_file = file
                        break
                
                plot_phonon_bands(bands_file, labels=[label], save_path=str(compound.get_chemical_formula(empirical=True)) + '_harm.png')
                

            
            
        elif args.OPTION.lower() in ['anharmonic', 'anharm']:
            
            if args.TEMPERATURE == []:
            
                print('Plotting the anharmonic phonon dispersion.\n')
                
                if args.PREFIX == []:
                    label1 = chem(str(compound.get_chemical_formula(empirical=True)))+' (anharmonic)'
                    label2 = chem(str(compound.get_chemical_formula(empirical=True)))+' (harmonic)'
                else:
                    label1 = chem(args.PREFIX[0]) +' (anharmonic)'
                    label2 = chem(args.PREFIX[0])+' (harmonic)'
    
    
                if calculation.prefix_phonon_bands != None and calculation.prefix_scph_bands != None:
                    plot_phonon_bands(calculation.prefix_scph_bands, harmData=calculation.prefix_phonon_bands, labels=[label1, label2], save_path=str(compound.get_chemical_formula(empirical=True)) + '_anharm.png')
                else:
                    
                    found_harm = False
                    found_anharm = False
                    
                    for file in os.listdir(cwd):
                        if file.find('.bands') !=-1:
                            harm_bands_file = file
                            found_harm = True
                        elif file.find('.scph_bands') !=-1:
                            anharm_bands_file = file
                            found_anharm = True
                        
                        if found_harm == True and found_anharm == True:
                            break
                    
                    plot_phonon_bands(anharm_bands_file, harmData=harm_bands_file, labels=[label1, label2], save_path=str(compound.get_chemical_formula(empirical=True)) + '_anharm.png')
            
            else:
                
                print('Plotting the harmonic and anharmonic dispersion at different K.\n')
                
                labels = []
                
                for temp in args.TEMPERATURE:
                    if args.PREFIX == []:
                        labels.append(chem(str(compound.get_chemical_formula(empirical=True))) + f' {temp} K')
                    else:
                        labels.append(chem(args.PREFIX[0]) +f' {temp} K')


                    
                    #labels.append(chem(str(compound.get_chemical_formula(empirical=True))) + f' {temp} K')

                if args.PREFIX == []:
                    labels.append(chem(str(compound.get_chemical_formula(empirical=True)))+' (harmonic)')
                else:
                    labels.append(chem(args.PREFIX[0])+' (harmonic)')

                #labels.append(chem(str(compound.get_chemical_formula(empirical=True)))+' (harmonic)')
                

                if calculation.prefix_phonon_bands != None and calculation.prefix_scph_bands != None:
                    
                    plot_phonon_bands(calculation.prefix_phonon_bands, harmData=calculation.prefix_phonon_bands, labels=labels, save_path=str(compound.get_chemical_formula(empirical=True)) + '_anharmTemp.png', getTemp=args.TEMPERATURE)
                else:
                    
                    found_harm = False
                    found_anharm = False
                    
                    for file in os.listdir(cwd):
                        if file.find('.bands') !=-1:
                            harm_bands_file = file
                            found_harm = True
                        elif file.find('.scph_bands') !=-1:
                            anharm_bands_file = file
                            found_anharm = True
                        
                        if found_harm == True and found_anharm == True:
                            break

                    plot_phonon_bands(anharm_bands_file, harmData=harm_bands_file, labels=labels, save_path=str(compound.get_chemical_formula(empirical=True)) + '_anharmTemp.png', getTemp=args.TEMPERATURE)

            
            
        elif args.OPTION.lower() == 'pes':
            
            #x_values = [-0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5]
            #q_list_eval = ['-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0.1', '0.2', '0.3', '0.4', '0.5']
            
            q_list_eval = [str(round(q, 5)) for q in np.array(np.arange(args.Q_VALUES[0], args.Q_VALUES[1]+args.Q_VALUES[2], args.Q_VALUES[2])) if round(q, 2) not in [0.0, -0.0]]

            x_values = [round(q, 5) for q in np.array(np.arange(args.Q_VALUES[0], args.Q_VALUES[1]+args.Q_VALUES[2], args.Q_VALUES[2])) if round(q, 2) not in [0.0, -0.0]]

            
            fig = plt.figure(figsize=(10,5))
            ax = fig.add_subplot()
            
            cmap = plt.get_cmap('gist_rainbow')
            
            colors = [cmap(value) for value in np.linspace(0, 1, len(args.PYRO_MODE))]
            
            for modeNum, mode in enumerate(args.PYRO_MODE):
                
                cur_mode_p = []
                
                for q in q_list_eval:
                    
                    for file in os.listdir():
                        
                        if file.find(f'mode{mode}_{q}.prop.out') != -1:
                            cur_out_file = file
                            # print(file)
                            # # print()
                            break
                        
                    found_polarization_values = False

                    with open(cur_out_file, 'r') as file:
                        for line in file:
                            if found_polarization_values == True and len(line.split()) > 3:
                                split_line = line.split()
                                cur_mode_p.append([float(split_line[1]), float(split_line[2]), float(split_line[3])])
                                break
                            
                            if line .find('POLARIZATION, P,') != -1:
                                found_polarization_values = True
            
                # x_values = [item[0] for item in cur_mode_p]
                # y_values = [item[1] for item in cur_mode_p]
                z_values = [item[2] for item in cur_mode_p]
                
                
                ax.plot(x_values, z_values, marker='x', label=f'mode {mode}', color=colors[modeNum])
                
                #ax.annotate(f'mode {mode})', xy=(x_values[-1], z_values[-1]), xytext=(x_values[-1]+0.05, z_values[-1]))


                
                
            ax.set_xlabel(r'Displacement (amu$^{1/2} \AA$)', fontsize=12)
            ax.set_ylabel(r'P (Cm$^{-2}$)', fontsize=12)
            
            # plt.title('orthorhombic BaTiO$_{3}$ (crystallographic supercell, noTranspose matrix)')



            

            plt.legend(bbox_to_anchor=(1.5, 1), loc='upper right')
            plt.tight_layout()

            plt.savefig(f'{str(compound.get_chemical_formula(empirical=True))}_pes.png')

            
            
        elif args.OPTION.lower() == 'pyro':
            
            fig = plt.figure(figsize=(8,5))
            
            ax = fig.add_subplot()
            
            pyro_all = []

            
        
            for mode in args.PYRO_MODE:
                cur_mode_Ps = []
                
                for temp in args.TEMPERATURE:
                    
                    for file in os.listdir(cwd):
                        if file.find('_%s'%temp)!= -1 and file.find('mode'+str(mode)+'_') != -1 and file.find('.out') != -1:
                            cur_out_file = file                        
                            break

                    found_polarization_values = False

                    with open(cur_out_file, 'r') as file:
                        for line in file:
                            if found_polarization_values == True and len(line.split()) > 3:
                                split_line = line.split()
                                cur_mode_Ps.append([float(split_line[1]), float(split_line[2]), float(split_line[3])])
                                break
                            
                            if line .find('POLARIZATION, P,') != -1:
                                found_polarization_values = True
            
                x_values = [item[0] for item in cur_mode_Ps]
                y_values = [item[1] for item in cur_mode_Ps]            
                z_values = [item[2]/40 * 10**6 for item in cur_mode_Ps]
                
                ax.plot(args.TEMPERATURE, z_values, marker='x', label=f'mode {mode}')
                
                pyro_all.append(z_values)
                
                #ax.annotate(f'mode {mode}', xy=(args.TEMPERATURE[0], z_values[0]), xytext=(args.TEMPERATURE[0]+1, z_values[0]))
                
            pyroTotal = []
            
            for tempNum, temp in enumerate(args.TEMPERATURE):
                
                cuPyro = 0
                
                for modeNum, mode in enumerate(args.PYRO_MODE):
                    
                    cuPyro += pyro_all[modeNum][tempNum]
                
                pyroTotal.append(cuPyro)
            
            ax.plot(args.TEMPERATURE, pyroTotal, marker='x', label=f'total p$^{(1)}$', color='black')
            
            ax.set_xlabel('Temperature (K)', fontsize=12)
            ax.set_ylabel(r'p ($\mu$Cm$^{-2}$K$^{-1}$)', fontsize=12)

            plt.legend(bbox_to_anchor=(1.25, 1), loc='upper right')
            
            plt.tight_layout()
            plt.savefig(f'{str(compound.get_chemical_formula(empirical=True))}_pyro.png')

            
            plt.show()

        
    
        """
    STEP eigenvector: Is supposed to handle the dfc2 script and create the inputs for evec generation for anphon
    """
    
    elif args.STEP.lower() in ['eigenvector', 'evec']:
        
        print('Creating the evecXXX.in file to run with anphon for temperatures %s.\n'%args.TEMPERATURE)
        
        if calculation.temperatures == None:
            calculation.temperatures = args.TEMPERATURE
        
        #  make sure that the length of the FNAME and temperature listst are the same -> add empty elements that will be renamed according to temp list        
        if len(args.FNAME) != len(args.TEMPERATURE):
            for temp in args.TEMPERATURE[len(args.FNAME)::]:
                args.FNAME.append('')
        
        #  same for prefix
        if len(args.PREFIX) != len(args.TEMPERATURE):
            for temp in args.TEMPERATURE[len(args.PREFIX)::]:
                args.PREFIX.append('')
                
                
        XML_harm = get_fcsXML(fcs2=True, fcs3=False, fcs4=False)

        for tempNum, temp in enumerate(args.TEMPERATURE):
            
            if args.FNAME[tempNum] == '':
                args.FNAME[tempNum] = f'evec{temp}.in'
                
            if args.PREFIX[tempNum] == '':
                args.PREFIX[tempNum] = str(compound.get_chemical_formula(empirical=True)) + f'_{temp}'
                
            if calculation.prefix_dfc2 != None:
                scph_dfc2_file = calculation.prefix_dfc2
            else:
                for file in os.listdir(cwd):
                    if file.find('.scph_dfc2') !=-1:
                        scph_dfc2_file = file
                        break
                        
            
            subprocess.run('dfc2', shell=True, input=f'{XML_harm[0]}\n{args.PREFIX[tempNum]}.xml\n{scph_dfc2_file}\n{temp}'.encode('utf-8'))
            
            print()
            
            
            for file in os.listdir(cwd):
                if file.find(f'_{temp}')!=-1 and file.find('.xml') !=-1:
                    xml_file = file
                    break

            evec_in = {'PREFIX': [args.PREFIX[tempNum]],'MODE': ['phonons'], 'PRINTEVEC': [1], 'FCSXML': [xml_file]}
            
            if args.STANDARDIZE == 'prim':

                write_ANPHON(compound, file_name=args.FNAME[tempNum], command_list=evec_in, get_XML_files=False, contPath=True, source=args.KPOINT_SOURCE)
                
            elif args.STANDARDIZE == 'crys':
                
                write_ANPHON(compound, file_name=args.FNAME[tempNum], command_list=evec_in, get_XML_files=False, contPath=True, source=args.KPOINT_SOURCE, stdIn = 'crys')
                


        
        #  this part of the script will loop throuh the set temperatures, since it has to be done for each temperature
        # for tempNum, temp in enumerate(args.TEMPERATURE):
        #     if args.FNAME[tempNum] == '':
        #         args.FNAME[tempNum] = f'evec{temp}.in'
                
        #     if args.PREFIX[tempNum] == '':
        #         args.PREFIX[tempNum] = str(compound.get_chemical_formula(empirical=True)) + f'_{temp}'
            
            
        #     #  look for the right FCSXML files that have the temp in their name (not via the get_FCSXML() function, just by name)
        #     #  -> seems to be some problem with finding the right files if not named "comp_temp.xml"...
        #     for file in os.listdir(cwd):
        #         if file.find(str(temp))!=-1 and file.find('.xml') !=-1:
        #             xml_file = file
        #             break

        #     evec_in = {'PREFIX': [args.PREFIX[tempNum]],'MODE': ['phonons'], 'PRINTEVEC': [1], 'FCSXML': [xml_file]}
            
        #     write_ANPHON(compound, file_name=args.FNAME[tempNum], command_list=evec_in, get_XML_files=False, contPath=True, source=args.KPOINT_SOURCE)
            

    elif args.STEP.lower() in ['test', 't']:
        
        evec0 = []
        
        #mode = args.PYRO_MODE[0]
        mode = 6
        
        for temp in args.TEMPERATURE:
            
            cur_evec_diff = []
            
            #get the corredsponding evec files
            for file in os.listdir(cwd):
                if file.find(str(temp))!=-1 and file.find('.evec') !=-1:
                    evec_file = file
                    break
                
            from ASE_Alamode_Interface import Alamode
            
            ALM = Alamode(compound)
            
            ALM.read_evec_gamma(evec_file)
            
            print('Length of Alamode evecs: %s'%len(ALM._evec))
            
            print('Mode 6: \n %s'%ALM._evec[0][5])
            print()
            
            

            
            # set the initial evec as the basic one
            
            if temp == args.TEMPERATURE[0]:
                evec0 = ALM._evec[0]
                print('saved the first temp evec.\n')
                
            else:
                
                for num, ev in enumerate(ALM._evec[0]):
                    print(evec0[num])
                    print(ev)
                    print()
                    cur_evec_diff = evec0[num]-ev
                    
                
                
        print('Difference:')
        print(cur_evec_diff)
        print()
        print('Sum of difference at evec 6:')
        print(sum(cur_evec_diff[0][6]))            
        print('---')
                
            
            
    
        
        
        
        
        
        




        """
        This is left of the old script...
        """




    #  add step 53: plotting only the harmonic phonons and the 300K anharmonic phonons
    elif args.STEP == 53:
        print('Plotting the harmonic and anharmonic dispersion at 300 K.\n')
        
        if args.PREFIX == []:
            label1 = chem(str(compound.get_chemical_formula(empirical=True)))+' 300 K (anharmonic)'
            label2 = chem(str(compound.get_chemical_formula(empirical=True)))+' (harmonic)'
        else:
            label1 = chem(args.PREFIX[0]) +' 300 K (anharmonic)'
            label2 = chem(args.PREFIX[0])+' (harmonic)'


        if calculation.prefix_phonon_bands != None and calculation.prefix_scph_bands != None:
            
            plot_phonon_bands(calculation.prefix_phonon_bands, harmData=calculation.prefix_phonon_bands, labels=[label1, label2], save_path=str(compound.get_chemical_formula(empirical=True)) + '_anharm300K.png', get300=True)
        else:
            
            found_harm = False
            found_anharm = False
            
            for file in os.listdir(cwd):
                if file.find('.bands') !=-1:
                    harm_bands_file = file
                    found_harm = True
                elif file.find('.scph_bands') !=-1:
                    anharm_bands_file = file
                    found_anharm = True
                
                if found_harm == True and found_anharm == True:
                    break

            plot_phonon_bands(anharm_bands_file, harmData=harm_bands_file, labels=[label1, label2], save_path=str(compound.get_chemical_formula(empirical=True)) + '_anharm.png', get300=True)



    #  add some routine to set pairs/sets of temp and mode for compounds where modes switch for temperatures
    #  maybe add some function that can anlyse the modes based on the frequqnecies given by the frequqncy calculation
    #  -> read in the frequqncy output from crystal
    #  -> get all frequqncies with A1(g) symmetry and the corresponding cm^-1
    #  -> get the eigenvectors and their cm^-1 from the .log file
    #  -> compare things by weeding out the potentially E frequqncies -> list of remaining non-degenerate ones; cut first 3 cause of rotation
    #  -> assign them via user input or by creating an excel table
    #  -> get list of eigenvectors with similar +/- distribution
    #  -> let user select modes for maybe 300K or after input that shows similar pattern to the frequqncy distribution of the CRYSTAL output
    #  -> signs of the eigenvectors from that input will be tracked along the different temperatures (and flipped right), such that -temp and -mode are distinctive pair lists for each mode (naming will be handeled automatically for the right 300K mode, so if it changes the mode for another temp it will still be named after the mode at 300K; maybe documented in an output)
    
    elif args.STEP == 7:  #  add routine here that automatically gets the polarization direction based on the space group (probably saved in cell) and relevant modes for each .evec file
        
        print('Create the random_normalcoord_pyro displacements and their .d12 files for:')
        print('Temperatures: %s'%args.TEMPERATURE)
        print('Modes: %s\n'%args.PYRO_MODE)
        
        if calculation.pyro_modes == None:
            calculation.pyro_modes = args.PYRO_MODE
        
        temps = []
        evec_files = []
        pyro_disp = [[] for i in args.PYRO_MODE]
        
        # print(pyro_disp)
        
        if args.FLIP == True:
            flip_ext = '_flip'
        else:
            flip_ext = ''
        
        #  create matrix with -20/+20 temperatures        
        for temp in args.TEMPERATURE:
            temps.append(temp-20)
            temps.append(temp+20)
            
            for file in os.listdir(cwd):
                if file.find(str(temp))!=-1 and file.find('.evec') !=-1:
                    evec_file = file
                    break
            
            if len(evec_files) > 0 and evec_file == evec_files[-1]:
                raise FileNotFoundError(f'The evec file for {temp} was not found!')
            
            evec_files.append(evec_file)
            evec_files.append(evec_file)
        
        
        for file in os.listdir(cwd):
            if file == 'CRY_TEMP_PYRO':
                crystal_template = file
                break
            else:
                crystal_template = None

        
        #  This step also needs to loop trough modes for each of the +/-temperatures
        for modeNum, mode in enumerate(args.PYRO_MODE):        
            for tempNum, temp in enumerate(temps):
                
                cur_pyro_disp = displace(compound_sc, 'random_normalcoord_pyro', evec_file = evec_files[tempNum], temp=temp, pyro_mode=mode, primitive_aseAtom=compound, flip=args.FLIP)
                
                pyro_disp[modeNum].append(cur_pyro_disp[0])
                
                cry_file_name = f'{temp}_mode{mode}{flip_ext}.d12'
                
                if calculation.is_rhomb == True:
                    add_tags_to_crystal(cur_pyro_disp[0])
                    
                    write_d12_GEOM_from_aseAtoms(cur_pyro_disp[0], d12_filename=cry_file_name, P1=True, verbosity=0, template=crystal_template, numDig=8, external=True)
                
                    write_crystal(f'{temp}_mode{mode}{flip_ext}.ext', cur_pyro_disp[0])
                
                else:
                    write_d12_GEOM_from_aseAtoms(cur_pyro_disp[0], d12_filename=cry_file_name, P1=True, verbosity=0, template=crystal_template, numDig=8)
                
                print(f'Displacement for mode {mode} at {temp} K created.')
        
        print()
        write_json(open('pyro_displacements', 'w'), pyro_disp)
        
        # maybe take this step as extra step?
        if args.DIRECTION_CHECK == True:
            
            sgNum = check_symmetry(compound).number
            calculation.space_group = sgNum
            
            
            if args.POLARIZATION_ORIENTATION == None:
                args.POLARIZATION_ORIENTATION = get_polarization_direction_from_sg(sgNum)

            
            
            print('Performing a direction check for the created displacements.\n')
            
            if args.POLARIZATION_ORIENTATION.lower() == 'x':
                coord_pos = 0
            elif args.POLARIZATION_ORIENTATION.lower() == 'y':
                coord_pos = 1
            elif args.POLARIZATION_ORIENTATION.lower() == 'z':
                coord_pos = 2
            else:
                raise ValueError("Preferred orientation has to be 'x', 'y' or 'z'!")
                
            calculation.pol_orient = args.POLARIZATION_ORIENTATION
            
        
            # from cell and supercell get the dimension of the supercell
            dim = compound_sc.cell.cellpar()/compound.cell.cellpar()
            
            first_positions = [int(i*dim[0]*dim[1]*dim[2]) for i in range(0, len(compound.positions), 1)]
            disp_sign = [[] for i in args.PYRO_MODE]
            
            elements = [el for el in set(compound.numbers)]
            
            #  sort the coordinates according to the Alamode way (at some point that shoul be made unnecessary)
            sc_positions = []
            
            for atom in elements:
                for numNum, num in enumerate(compound_sc.numbers):
                    if atom == num:
                        sc_positions.append(compound_sc.positions[numNum])
                        
                        
            #  difference is always origPos - dispPos -> get +/-1 as values, counts which one occurs more often -> is set as "right direction" for that mode -> the displacements that deviate are flipped
            for modeNum, mode in enumerate(args.PYRO_MODE):
                
                for disp in pyro_disp[modeNum]:
                    
                    diff = sc_positions[first_positions[-1]][coord_pos] - disp.positions[first_positions[-1]][coord_pos]
                    
                    if diff > 0:
                        disp_sign[modeNum].append(1)
                    elif diff < 0:
                        disp_sign[modeNum].append(-1)
                    else:
                        disp_sign[modeNum].append(0)
                        print('Displacement is 0 or something went wrong.')
                                    
                        
                more_of_sign = max(set(disp_sign[modeNum]), key=disp_sign[modeNum].count)
                
                
                for signNum, sign in enumerate(disp_sign[modeNum]):
                    if sign != more_of_sign:
                        # print(f'For mode {mode} at {temps[signNum]} evec will be flipped!\n')
                        # print(f'Reading {evec_files[signNum]} as evec file.\n')
                        
                        new_pyro_disp = displace(compound_sc, 'random_normalcoord_pyro', evec_file = evec_files[signNum], temp=temps[signNum], pyro_mode=mode, primitive_aseAtom=compound, flip=True)
                        
                        pyro_disp[modeNum][signNum] = new_pyro_disp
                        
                        cry_file_name = f'{temp}_mode{mode}_flip.d12'
                        
                        
                        add_tags_to_crystal(new_pyro_disp[0])
                        
                        write_d12_GEOM_from_aseAtoms(new_pyro_disp[0], d12_filename=cry_file_name, P1=True, verbosity=0, template=crystal_template, numDig=8, external=True)

                            
                        
                        print(f'New displacement for mode {mode} at {temp} K with flipped eigenvector created.')
                

                
                    
            # dump the new list of siplacements again
            write_json(open('pyro_displacements_cor', 'w'), pyro_disp)
            print()
            
            


    elif args.STEP == 8:
        print('Create the POLARI .d3 files for each random_normalcoord_pyro displacement for:')
        print('Temperatures: %s'%args.TEMPERATURE)
        print('Modes: %s\n'%args.PYRO_MODE)
        
        
        temps = []
        
        #  create matrix with -20/+20 temperatures        
        for temp in args.TEMPERATURE:
            temps.append(temp-20)
            temps.append(temp+20)
            
            
        for mode in args.PYRO_MODE:
            for temp in temps:        
        
                polari_file_name = f'{temp}_mode{mode}.d3'
                
                with open(polari_file_name, 'w') as polari_file:
                    polari_file.write('NEWK\n')
                    polari_file.write('4 4\n')
                    polari_file.write('1 0\n')
                    polari_file.write('POLARI\n')
                    polari_file.write('END')
                    polari_file.close()

    elif args.STEP == 70:
        print('Create the crystallographic cell pes displacements and their .d12 files for:')
        
        q_list = [str(np.round(q, 2)) for q in np.array(np.arange(args.Q_VALUES[0], args.Q_VALUES[1]+args.Q_VALUES[2], args.Q_VALUES[2]))]
        
        
        print('q_list: %s'%q_list)
        print('Modes: %s\n'%args.PYRO_MODE)
        
        if calculation.pyro_modes == None:
            calculation.pyro_modes = args.PYRO_MODE
        
        ndata = int(((args.Q_VALUES[1] - args.Q_VALUES[0]) / args.Q_VALUES[2]) + 1)

        print(ndata)
        print()
        
        
        evec_files = []
        pes_disp = [[] for i in args.PYRO_MODE]

        if args.FLIP == True:
            flip_ext = '_flip'
        else:
            flip_ext = ''

        temp = args.TEMPERATURE[0]

        for file in os.listdir(cwd):
            if file.find(str(temp))!=-1 and file.find('.evec') !=-1:
                evec_file = file
                break

        for file in os.listdir(cwd):
            if file == 'CRY_TEMP_PYRO':
                crystal_template = file
                break
            else:
                crystal_template = None

        for modeNum, mode in enumerate(args.PYRO_MODE):

            cur_disp_pes = displace(compound, 'pes', evec_file = evec_file, temp=temp, pyro_mode=mode, primitive_aseAtom=compound, qrange=[args.Q_VALUES[0], args.Q_VALUES[1]], ndata=ndata)
                
            if args.CRYSYS == 'ortho':
                print('Crystal system orthorhombic selected.\n')
                
                matrix = [[1,0,0], [0,1,-1], [0,1,1]]
                
            elif args.CRYSYS == 'rhomb':
                print('Crystal system rhombohedral selected.\n')
                
                matrix = [[1,-1,0], [0,1,-1], [1,1,1]]

            elif args.CRYSYS == 'tetra':
                
                matrix = [[1,0,0], [0,1,0], [0,0,1]]
                
            else:
                print('The crystal class is not recognized by the script, please add the right matrix for it.', '\n')
                matrix = [[1,0,0], [0,1,0], [0,0,1]]
                    
            
            for dispNum, disp in enumerate(cur_disp_pes):
                
                add_tags_to_crystal(disp)
                
                pes_sc = make_supercell(disp, matrix, order='atom-major')
                
                pes_disp[modeNum].append(pes_sc)
                    
                cry_file_name = f'pes_disp_{q_list[dispNum]}_mode{modeNum+1}.d12'
                prop_file_name = f'pes_disp_{q_list[dispNum]}_mode{modeNum+1}.d3'
            
            
                write_d12_GEOM_from_aseAtoms(pes_sc, d12_filename=cry_file_name, P1=True, verbosity=0, template=crystal_template, numDig=8, external=True)
            
                #write_crystal(f'{temp}_mode{mode}{flip_ext}.ext', pyro_sc)
                
                with open(prop_file_name, 'w') as polari_file:
                    polari_file.write('NEWK\n')
                    polari_file.write('4 4\n')
                    polari_file.write('1 0\n')
                    polari_file.write('POLARI\n')
                    polari_file.write('END')
                    polari_file.close()
                
                print(f'Displacement for mode {mode} at q_vale {q_list[dispNum]} created.')
        
        print()
        write_json(open('pes_displacements', 'w'), pes_disp)








    calculation.save_pyroCalculation('pyroCalc_progress')
    
    obj_dict = vars(calculation)
    
    # for attribute in vars(calculation):
    #     print(attribute)
    #     print(obj_dict[attribute])


    #  write a file to document the progress of calculations:
    with open('calc_documentation.txt', 'a') as doc_file:
        doc_file.write('Calculation step: %s\n'%args.STEP)
        doc_file.write('Saved parameters:\n')

        for attribute in vars(calculation):
            if obj_dict[attribute] != None:
        
                doc_file.write('%s: %s\n'%(attribute, obj_dict[attribute]))
                
        if args.TEMPERATURE != [] or args.PYRO_MODE != []:
            doc_file.write('-----------------------------------\n')
            doc_file.write('Calculations done for:\n')
        
        if args.TEMPERATURE != []:
            doc_file.write('Temperatures: %s\n'%args.TEMPERATURE)
        
        if args.PYRO_MODE != []:
            doc_file.write('Pyro modes: %s\n'%args.PYRO_MODE)
            
        doc_file.write('\n')



    """
    TBC with this:
        - go through whole workflow with ortho BTO to see where possible bugs are (BTO maybe 2 temp & 2 mode)
        
        - discuss and decide how that flipping of eigenvectors should be handeled and try to implement it into the script
            -> two sets that are handled accordingly in step 10 (where the numbers are plotted)
            -> or define (at least for compounds with only one polarisation axis possible) some plane by the other two and get relative position of the remaining component
            -> depending on +/- add to the right set
            -> with this steps 7 and onward need to be handled for both cases sepperately
            
        - maybe think about implementing something that picks the right modes according to the eigenvectors and the polarization direction based on the spacegroup
            -> sets direction(s)
            -> chooses appropriate modes which vibrate along the right direction(s) (maybe look at them relatively if e.g. the eigenvector along pol. direction is much larger than the others)
            -> set them as modes to take into account
            
        - implement some way to handle the case, where the actual number of the mode changes for the different temperatures (guess this could happen with energy degenerate modes at higher temperatures) -> how likely is that -> QUESTION?
        
        - implement some check for too high P_s (if one would end up at the higher branches) so that a correction with the quantum would be needed
        
        - make a table/file/list/smthg similar that has a list of the possible directions for the right spacegroups (probably smthg. to hardcode...)
        
    """








    #  simple code to plot the P_s values -> incorporate into script or make own function that reads the values from the prop.o output files...

    # temps = [300, 400, 500]
    
    # mode6 = [0.011567, 0.010352, 0.009377]
    # mode12 = [0.001828, 0.001878, 0.001829]
    # mode15 = [-0.00347, -0.00423, -0.00450]
    
    # mode6_flip = [-0.011496, -0.010327, -0.009377]
    # mode12_flip = [-0.001901, -0.001970, -0.001960]
    # mode15_flip = [0.002850, 0.003426, 0.003602]
    
    
    # fig = plt.figure()
    
    # ax = fig.add_subplot()
    
    # ax.plot(temps, mode6, marker='x', label='mode 6')
    # ax.plot(temps, mode12, marker='x', label='mode 12') 
    # ax.plot(temps, mode15, marker='x', label='mode 15')
    
    # plt.gca().set_prop_cycle(None)

    # ax.plot(temps, mode6_flip, marker='x')
    # ax.plot(temps, mode12_flip, marker='x')
    # ax.plot(temps, mode15_flip, marker='x')

    # ax.set_xlabel('temperature', fontsize=12)
    # ax.set_ylabel(r'P$_s$', fontsize=12)
    
    # ax.legend()
    
    # plt.tight_layout()
    # plt.savefig('BTO_pyro_all.png')
    
    # plt.show()
    
    
    eos_raccoon()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
