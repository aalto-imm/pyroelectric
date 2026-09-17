"""
This interface will replace the displace.py and extract.py scripts that are normally used to interact with the Alamode program

It will contain the main class 'ALAMODE' which will load upon calling all relevant information from an aseAtoms object (numbers positions, etc.) -> For now coordinates and co will be sorted to handle them more human readable, later on those options will prob. be omitted

it will include two functions displace() and extract() that will handle creating displacements and extracting forces depending on the calculator set for the aseAtoms object

this will be a parent class for the two classes ALM and ANPHON, that handle input creation for the two Alamode inputs


"""
import numpy as np
import random
import math
import cmath

from ase import Atoms
from ase.data import chemical_symbols, atomic_numbers

from ase.units import Bohr, Rydberg, kJ, Hartree, nm, mol


from cif2ALM import get_fcsXML, get_keywords_from_template, ALMgeneral, interaction, optimize, write_ALM_input

from cif2ANPHON import ANPHONgeneral, scph, qha, relax, analysis, write_ANPHON_input

#from TestAlamodeInput.ASE_CRTYSTAL_interface import get_energy_CRYSTAL, get_forces_CRYSTAL
from ASE_CRTYSTAL_interface import get_energy_CRYSTAL, get_forces_CRYSTAL

from cif2D12 import standardize_atomsAtom_crstalStructure, correct_sorting_of_stand_to_prev_ase, get_spacegroup_from_primitive

## Needed to load for testing:
from ase.io import cif
from ase.build import make_supercell

#from tblite.ase import TBLite





class Alamode():
    
    def __init__(self, aseAtom, verbosity=0, length_unit='angström', energy_unit='hartree', force_unit='hartree/bohr', primitive_aseAtom=None, **kwargs):
        #self.lattice_vector = aseAtom.get_cell().transpose()
        self.lattice_vector = np.reshape([vec for vec in aseAtom.get_cell()], (3,3))
        #self.cell = aseAtom.cell.cellpar()
        self.inverse_lattice_vector = np.linalg.inv(self.lattice_vector)
#        self.elements = [chemical_symbols[el] for el in set(aseAtom.numbers)]
        self.elements = [chemical_symbols[el] for el in list(dict.fromkeys(aseAtom.numbers))]

        #  [chemical_symbols[el] for el in list(dict.fromkeys(LNO.numbers))] -> chnage to this, is always ordered according to the first occurence of each (otherwise set() could mess up stuff)
        self.nat = len(aseAtom.positions)
        self.x_fractional = self.get_sorted_coordinates(aseAtom)
        self.x_cartesian = self.get_sorted_coordinates(aseAtom, cartesian=True)
        self.numbers = self.get_sorted_numbers(aseAtom)
        self.nat_elem = [self.numbers.count(atomic_numbers[el]) for el in self.elements]
        self.kd = np.array([self.elements.index(chemical_symbols[i]) for i in self.numbers], dtype=int)
        self.symbols = [chemical_symbols[i] for i in self.numbers]
        self.calculator = aseAtom.calc
        
        self.verbosity = verbosity  ##  should handle all print() commands -> nothing if 0, everything if 1
        
        #  at some point go over this whole unit conversion thing and make it easier...
        self.length_unit = length_unit
        self.energy_unit = energy_unit
        
        self._legth_conversion_factor = 1.0
        self._energy_conversion_factor = 1.0
        self._force_conversion_factor = 1.0
        
        #  attributes that are set if an evec file is read in:
        self._qpoints = None
        self._omega2 = None
        self._evec = None
        self._qlist_real = None
        self._qlist_uniq = None
        self._mass = None
        self._nmode = None
        
        #  attributes that are needed for the random_normalcoord_pyro option
        self._mapping_s2p = None
        self._classical = False
        self._RYDBERG_TO_JOULE = 4.35974394e-18 / 2.0
        self._mapping_shift = None
        amu = 1.660538782e-27
        electron_mass = 9.10938215e-31
        self._AMU_RYD = amu / electron_mass / 2.0
        self._K_BOLTZMANN = 1.3806488e-23


        #  attributes that get loaded if a primitve aseAtom is given:
        if primitive_aseAtom != None:
            self._primitive_lattice_vector = np.reshape([vec for vec in primitive_aseAtom.get_cell()], (3,3))
            self._inverse_primitive_lattice_vector = np.linalg.inv(self._primitive_lattice_vector)
            self._nat_primitive = len(primitive_aseAtom.positions)
            
            self._xp_fractional = self.get_sorted_coordinates(primitive_aseAtom, primitive=True)
            self._numbers_primitive = self.get_sorted_numbers(primitive_aseAtom)
            self._primitive_kd = np.array([self.elements.index(chemical_symbols[i]) for i in self._numbers_primitive], dtype=int)
            self._generate_mapping_s2p()
            self._find_commensurate_q()
            

        
        #self._BOHR_TO_ANGSTROM = 0.5291772108
        #self._RYDBERG_TO_EV = 13.60569253
        #self._HARTREE_TO_EV = 27.211386245988

        
        
    def refold(x):
        if x >= 0.5:
            # print('doing smthg with x>=0.5')
            return x - 1.0
        elif x < -0.5:
            # print('doing smthg with x<=-0.5')
            return x + 1.0
        else:
            # print('doing nothing...')
            return x
        
    def _get_atomic_forces(self, aseAtom, CRYSTAL=False, **kwargs):  #  uses the set calculator to calculate the atomic forces (gives error message if no calculator is set); the CRYSTAL option is a small workaraound to call our inbuild function that reads crystal outputs **kwargs can be passed e.g. as additional arguments for the calculator; try and find out which units this is...
        
        if type(aseAtom) == str or aseAtom.calc == None:
            if CRYSTAL == True:
                forces = get_forces_CRYSTAL(aseAtom)
                return forces
            else:
                print('Please set a calculator for your Atoms object!\n')
                return None
        else:
            forces = aseAtom.get_forces(**kwargs)
                       
            return forces


        
    def _get_total_energy(self, aseAtom, CRYSTAL=False, **kwargs):  #  gives the total energy of the system (; try and find out which units this is...); gives back None if no calculator is set
        
        if type(aseAtom) == str or aseAtom.calc == None:
            if CRYSTAL == True:
                etot = get_energy_CRYSTAL(aseAtom)
                return np.array(etot, dtype=np.float64)
            else:
                print('Please set a calculator for your Atoms object!\n')
                return None
        else:
            etot = aseAtom.get_total_energy(**kwargs)
            
            return np.array(etot, dtype=np.float64)

                      
        
    
    #  that two functions should take care of any kind of unit comversion, one should convert anything to bohr, the other should convert all energies to rydberg

    def set_length_conversion_factor(self, in_unit):  #  as length Alamode always uses the units for bohr, therefore possible lengths are converted to bohr
                
        if in_unit.lower() == 'bohr':
            
            self._legth_conversion_factor = 1.0
        
        elif in_unit.lower() in ['angström', 'ang', 'angstrom']:
        
            self._legth_conversion_factor = 1 / Bohr
            
        elif in_unit.lower() in ['nm', 'nanometer', 'nano meter']:  #  convert from nm to ang then to bohr
            self._legth_conversion_factor = 1 / Bohr * nm

        return
    
    
    def set_energy_conversion_factor(self, in_unit):  #  energies in Alamode are always in rydberg
        
        if in_unit.lower() == 'rydberg':
            self._energy_conversion_factor = 1.0
        
        elif in_unit.lower() == 'hartree':  #  conversion factor hartree to rydberg is 2
            self._energy_conversion_factor = 1 / Rydberg * Hartree  #  go from Hartree to eV then to Rydberg
            
        elif in_unit.lower() in ['ev', 'electronvolt', 'electron volt']:
            self._energy_conversion_factor = 1 / Rydberg
            
        elif in_unit.lower() in ['kj', 'kilo-joule', 'kilo joule', 'kilojoule']:
            self._energy_conversion_factor = 1 / Rydberg * kJ
            
        elif in_unit.lower() in ['kj/mol', 'kilo joule per mole', 'kj mol', 'kilojoule mole']:
            self._energy_conversion_factor = 1 / (Rydberg * mol) * kJ
            
        
        return
    
    #  simplified version of the read_evec function
    def read_evec_gamma(self, evec_file, flip=False):
        
        
        found_gamma = False
        
        with open(evec_file, 'r') as evec_file:
            for line in evec_file:
                
                if line.find('Number of phonon modes:') != -1:
                    nmode = int(line.split()[-1])
                    gamma_evecs = np.zeros((1, nmode, nmode), dtype=np.complex128)
                    omega2 = np.zeros((1, nmode))
                  
                #  this might be unnecessary, since masses are set anyway when the Alamode object is created
                if line.find('# Atomic masses') != -1:
                    mass = [float(line.split()[4]), float(line.split()[5]), float(line.split()[6])]
                
                if line.find('## kpoint       2 :') != -1:                    
                    found_gamma = False
                    break
                    
                if found_gamma == True:
                    if line.split() == []:
                        continue
                    elif line.find('###') != -1:
                        mode_num = int(line.split()[2]) - 1
                        omega2[0][mode_num] = float(line.split()[4])
                        counter = 0
                    
                    else:
                        if flip == True:
                            gamma_evecs[0][mode_num][counter] = complex(float(line.split()[0])*(-1), float(line.split()[1])*(-1))
                        else:
                            gamma_evecs[0][mode_num][counter] = complex(float(line.split()[0]), float(line.split()[1]))
                        counter += 1
                            
                                            
                if line.find('## kpoint       1 :') != -1:
                    xq = np.array([[float(line.split()[4]), float(line.split()[5]), float(line.split()[6])]])
                    found_gamma = True
                             
            evec_file.close
            
            
            #  The second part of the read evec_file function is mostly used in the random_noramlcoord displacement, since the q_lists are only used for calculations there
            #  Most of this stuff becomes redundant for the gamma point since nothing really happens since the number of q points there is only 1
            #  For the displace_random_normalcoord thing, the q_listst are sorted in such a way, that points that could be there twice (e.g. if a k-point is repeated in the path) would get sorted out from the list.

        # qlist_real = []
        # qlist_uniq = []
        
        # nq = 1

        # # Prepare q point list -> this is the most time consuming step...
        # for iq in range(nq):
        #     xq_tmp = xq[iq, :]  # list of 3 integs
        #     xq_minus = -xq_tmp
        #     xdiff = (xq_tmp - xq_minus) % 1.0
        #     norm = math.sqrt(np.dot(xdiff, xdiff))
        #     if norm < tol_zero:
        #         qlist_real.append(iq)

        #     flag_uniq = True

        #     for jq in qlist_uniq:
        #         xq_tmp2 = xq[jq, :]
        #         xdiff = (xq_tmp2 - xq_tmp) % 1.0
        #         xdiff2 = (xq_tmp2 + xq_tmp) % 1.0
        #         norm = math.sqrt(np.dot(xdiff, xdiff))
        #         norm2 = math.sqrt(np.dot(xdiff2, xdiff2))

        #         if norm < tol_zero or norm2 < tol_zero:
        #             flag_uniq = False
        #             break
        #     if flag_uniq:
        #         qlist_uniq.append(iq)
                
        # qlist_uniq = list(set(qlist_uniq) - set(qlist_real))

        # self._qpoints = xq
        # self._omega2 = omega2  #  is actually used somewhere
        # self._evec = gamma_evecs
        # self._qlist_real = qlist_real
        # self._qlist_uniq = qlist_uniq
        # self._mass = mass
        # self._nmode = nmode

        self._qpoints = xq
        self._omega2 = omega2  #  is actually used somewhere
        self._evec = gamma_evecs
        self._nmode = nmode
        self._mass = mass       
        
        return
    
    
    
    def read_evec_file(self, evec_file):
                
        tol_zero = 1.0e-3

        f = open(evec_file, 'r')

        # skip 10 lines
        for i in range(10):
            f.readline()

        nmode = int(f.readline().split(':')[1])  #  number of phonon modes (15 BTO)
        nq = int(f.readline().split(':')[1])     #  number of k-points -> a lot (2750 BTO)
        nkd = int(f.readline().split(':')[1])    #  number of atomic kinds
        mass = [float(t) for t in f.readline().split(':')[1].split()]  #  atomic masses
                
        # skip 3 lines
        for i in range(3):
            f.readline()

        omega2 = np.zeros((nq, nmode))

        evec = np.zeros((nq, nmode, nmode), dtype=np.complex128)
        xq = np.zeros((nq, 3))

        for iq in range(nq):
            xq_tmp = [float(a) for a in (f.readline().split(':')[1]).split()]
            xq[iq, :] = xq_tmp[:]
            for imode in range(nmode):
                omega2[iq, imode] = float(f.readline().split(':')[1])
                for jmode in range(nmode):
                    line = f.readline().split()
                    evec[iq, imode, jmode] = complex(
                        float(line[0]), float(line[1]))
                f.readline()
            f.readline()
        
        
        
        qlist_real = []
        qlist_uniq = []

        # Prepare q point list -> this is the most time consuming step...
        for iq in range(nq):
            xq_tmp = xq[iq, :]  # list of 3 integs
            xq_minus = -xq_tmp
            xdiff = (xq_tmp - xq_minus) % 1.0
            norm = math.sqrt(np.dot(xdiff, xdiff))
            if norm < tol_zero:
                qlist_real.append(iq)

            flag_uniq = True

            for jq in qlist_uniq:
                xq_tmp2 = xq[jq, :]
                xdiff = (xq_tmp2 - xq_tmp) % 1.0
                xdiff2 = (xq_tmp2 + xq_tmp) % 1.0
                norm = math.sqrt(np.dot(xdiff, xdiff))
                norm2 = math.sqrt(np.dot(xdiff2, xdiff2))

                if norm < tol_zero or norm2 < tol_zero:
                    flag_uniq = False
                    break
            if flag_uniq:
                qlist_uniq.append(iq)
                
        qlist_uniq = list(set(qlist_uniq) - set(qlist_real))

        f.close()

        self._qpoints = xq
        self._omega2 = omega2  #  is actually used somewhere
        self._evec = evec
        self._qlist_real = qlist_real
        self._qlist_uniq = qlist_uniq
        self._mass = mass
        self._nmode = nmode

        return


    def _n_bose(self, omega, temperature):
        if abs(temperature) < 1.0e-15 or omega < 1.0e-10:
            return 0.0
        else:
            temperature_au =  self._K_BOLTZMANN * temperature / self._RYDBERG_TO_JOULE
            x = omega / temperature_au

            return 1.0 / (math.exp(x) - 1.0)


    def _n_classical(self, omega, temperature):
        if abs(temperature) < 1.0e-15 or omega < 1.0e-10:
            return 0.0
        else:
            temperature_au = self._K_BOLTZMANN * temperature / self._RYDBERG_TO_JOULE
            return temperature_au / omega



    def _get_gaussian_sigma(self, temp, ignore_imag):
        """
         Computes the deviation of Q, i.e., sqrt{<Q^2>} from the phonon frequency
         and input temperature. Since omega is defined in the Rydberg atomic unit,
         the return value (sigma) is also in the Rydberg atomic unit (bohr*amu_ry^(1/2)).
        """

        nq = len(self._qpoints)
        nmode = self._nmode
        omega = np.zeros((nq, nmode))
        sigma = np.zeros((nq, nmode))

        for iq in range(nq):
            for imode in range(nmode):
                if self._omega2[iq, imode] < 0.0:
                    if ignore_imag:
                        omega[iq, imode] = 0.0
                        if self.verbosity > 0:
                            print("Warning: Detected imaginary mode at iq = %d, imode = %d.\n"
                                  "This more will be ignored.\n" % (iq + 1, imode + 1))
                    else:
                        omega[iq, imode] = math.sqrt(-self._omega2[iq, imode])
                        if self.verbosity > 0:
                            print("Warning: Detected imaginary mode at iq = %d, imode = %d.\n"
                                  "Use the absolute frequency for this mode.\n" % (iq + 1, imode + 1))
                else:
                    omega[iq, imode] = math.sqrt(self._omega2[iq, imode])

                if omega[iq, imode] > 1.0e-6:
                    if self._classical:
                        sigma[iq, imode] = math.sqrt(self._n_classical(
                            omega[iq, imode], temp) / omega[iq, imode])
                    else:
                        sigma[iq, imode] = math.sqrt(
                            (1.0 + 2.0 * self._n_bose(omega[iq, imode], temp)) / (2.0 * omega[iq, imode]))
                
        return sigma


    def _find_commensurate_q(self):

        tol_zero = 1.0e-3

        nqmax = self.nat // self._nat_primitive
        convertor = np.dot(self.inverse_lattice_vector,
                           self._primitive_lattice_vector)
        #print('\nConvertor in function:\n', convertor, '\n')
    
        nmax = 10
        qlist = []

        for i in range(3):
            for j in range(3):
                frac = abs(convertor[i, j])
                
                # print()
                # print('Vale of variable frac: ', frac)

                if frac < tol_zero:
                    convertor[i, j] = 0.0
                else:
                    found_nnp = False
                    for nnp in range(1, 1000):
                        if abs(frac * float(nnp) - 1.0) < tol_zero:
                            found_nnp = True
                            break
                    if found_nnp:
                        convertor[i, j] = np.sign(convertor[i, j]) / float(nnp)
                    else:
                        raise RuntimeError("Failed to express the inverse transformation matrix"
                                           "by using fractional numbers.\n\n"
                                           "Please make sure that the lattice parameters of \n"
                                           "the supercell and primitive cell are consistent.\n")

        comb = []
        for Lx in range(nmax):
            for Ly in range(nmax):
                for Lz in range(nmax):
                    for sx in (1, -1):
                        for sy in (1, -1):
                            for sz in (1, -1):
                                comb.append([Lx * sx, Ly * sy, Lz * sz])

        for entry in comb:
            vec = np.array([entry[0], entry[1], entry[2]])
            qvec = np.dot(vec, convertor) % 1.0

            for i in range(3):
                if qvec[i] >= 0.5:
                    qvec[i] -= 1.0
                elif qvec[i] < -0.5:
                    qvec[i] += 1.0

            new_entry = True

            for elem in qlist:
                diff = (elem - qvec) % 1.0

                for i in range(3):
                    if diff[i] >= 0.5:
                        diff[i] -= 1.0
                    elif diff[i] < -0.5:
                        diff[i] += 1.0

                norm = math.sqrt(np.dot(diff, diff))

                if norm < tol_zero:
                    new_entry = False
                    break

            if new_entry:
                qlist.append(qvec)

            if len(qlist) == nqmax:
                break

        self._commensurate_qpoints = qlist

        return


    def _generate_mapping_s2p(self):

        tol_zero = 1.0e-3
        convertor = np.dot(self.lattice_vector,
                           self._inverse_primitive_lattice_vector)

        # print('Convertor: ', convertor)

        for i in range(3):
            for j in range(3):
                convertor[i, j] = float(round(convertor[i, j]))

        shift = np.zeros((self.nat, 3))
        map_s2p = np.zeros(self.nat, dtype=int)

        for iat in range(self.nat):
            xtmp = self.x_fractional[iat, :]
            xnew = np.dot(xtmp, convertor)
            iloc = -1


            for jat in range(self._nat_primitive):
                xp = self._xp_fractional[jat, :]
                xdiff = np.array((xnew - xp) % 1.0)
                for i in range(3):
                    if xdiff[i] >= 0.5:
                        xdiff[i] -= 1.0
                diff = math.sqrt(np.dot(xdiff[:], xdiff[:]))
                # print('diff: ', diff)
                if diff < tol_zero:
                    iloc = jat

                    # print('\n', 'iloc: ', iloc)

                    break
            if iloc == -1:
                raise RuntimeError("Equivalent atom not found")

            map_s2p[iat] = iloc
            shift[iat, :] = [float(round(xnew[i] - self._xp_fractional[iloc, i]))
                             for i in range(3)]

        self._mapping_shift = shift
        self._mapping_s2p = map_s2p


        return


    def get_sorted_coordinates(self, aseAtom, cartesian=False, primitive=False):  #  returns the fractional coordinates of the atoms object, sorted makes stuff more human readable...
        
        if cartesian == False:
            positions = aseAtom.get_scaled_positions()
        else:
            positions = aseAtom.get_positions()
        
        xf_sort = []
        
        for atom in self.elements:
            for numNum, num in enumerate(aseAtom.numbers):
                if atomic_numbers[atom] == num:
                    xf_sort.append(positions[numNum])
        
        if primitive == False:
            return np.reshape(xf_sort, (self.nat, 3))
        else:
            return np.reshape(xf_sort, (self._nat_primitive, 3))
    
    
    def get_sorted_numbers(self, aseAtom):
        
        num_sort = []
        
        for atom in self.elements:
            for num in aseAtom.numbers:
                if atomic_numbers[atom] == num:
                    num_sort.append(num)
                    
        return num_sort
        
    
    
    #  disp_mode: sets the displacement mode; options are [pattern_file, random, md, random_normalcoord, random_normalcoord_pyro]
    #  file_pattern needed for disp_mode == 'pattern_file' -> ALM generated file that holds the information on harmonic displacements
    #  mag: magnitude of the displacement given in angström
    #  ndata: number of displacements that will be created
    #  mode: mode for the generation of the random displacements currently gauss (default) and uniform are supported
    #  evec: file with eigenvalue vectors provided by Alamode
    # temp: temperature
    
def displace(aseAtom, disp_mode, file_pattern=None, mag=0.02, ndata=1, mode='gauss', evec_file=None, temp=None, ignore_imag=False, primitive_aseAtom=None, classical=False, pyro_mode=1, in_unit='ang', flip=False, qrange=[-0.5,0.5], altALM=None, **kwargs):
    
    """ the altALM thing does not do anything worthy right now (pes scan looks way worse) """
              
    ALM = Alamode(aseAtom, primitive_aseAtom=primitive_aseAtom)
    
    
    ALM.set_length_conversion_factor(in_unit)
    
    displacements = []
        
    if disp_mode == 'pattern_file':
        if not file_pattern:
            raise RuntimeError('ALM generated pattern file must be given with "file_pattern" variable')
                    
        disps = []
        
        with open(file_pattern, 'r') as inFile:
            for line in inFile:
                line_split = line.split()

                if line_split[0].find(':') != -1 or line_split[0] == 'Basis':
                    continue
                else:
                    disps.append([float(i) for i in line_split])
                                        
        for disp in disps:
            #print(disp)
            
            cur_disp_x_frac = []
            
            for num, pos in enumerate(ALM.x_fractional):
                if num + 1 == int(disp[0]):
                    
                    disp_vec = [i*mag for i in disp[1::]]
                    
                    disp_x_frac = pos + np.dot(disp_vec, ALM.inverse_lattice_vector.T) ## !!transpose!!
                    #print(ALM.inverse_lattice_vector)
                    # print(num)
                    # print(pos)
                    # print(disp_x_frac)
                    # print()

                    # disp_x_frac = pos + np.dot(ALM.inverse_lattice_vector, disp_vec)
                                        
                else:
                    disp_x_frac = [i for i in pos]
                
                
                cur_disp_x_frac.append(disp_x_frac)

            disp_frac = np.reshape(cur_disp_x_frac, (ALM.nat, 3))
            
            #print(disp_frac)
        
            newAseDisp = Atoms(numbers=ALM.numbers, scaled_positions=disp_frac, cell=ALM.lattice_vector, pbc=True)
        
            displacements.append(newAseDisp)
        

    elif disp_mode == 'random':
        
        disp_random = np.zeros((ndata, ALM.nat, 3))
        
        if mode == 'gauss':
            for idata in range(ndata):
                for i in range(ALM.nat):
                    disp_xyz = [random.gauss(0.0, 1.0), random.gauss(0.0, 1.0), random.gauss(0.0, 1.0)]
                    norm = np.linalg.norm(disp_xyz)
                    disp_random[idata, i, :] = disp_xyz[:] / norm * mag
                    disp_random[idata, i] = np.dot(disp_random[idata, i],                          ALM.inverse_lattice_vector.T) ##!!transpose!!
                    
        elif mode == "uniform":
            for idata in range(ndata):
                for i in range(ALM.nat):
                    # Generate a random number following the Gaussian distribution
                    disp_xyz = [random.uniform(-mag, mag), random.uniform(-mag, mag), random.uniform(-mag, mag)]
                    # Transform to the fractional coordinate
                    disp_random[idata, i] = np.dot(disp_xyz[:],
                                                   ALM.inverse_lattice_vector.T) ##!!transpose!!
        else:
            raise RuntimeError("Invalid option for the random number distribution types.")


        for disp in disp_random:
            cur_disp_x_frac = []
            
            for coordNum, coords in enumerate(disp):
                cur_disp_x_frac.append(ALM.x_fractional[coordNum] + coords)
                
            disp_frac = np.reshape(cur_disp_x_frac, (ALM.nat, 3))
        
            newAseDisp = Atoms(numbers=ALM.numbers, scaled_positions=disp_frac, cell=ALM.lattice_vector, pbc=True)
        
            displacements.append(newAseDisp)

    
    #  add maybe later
    elif disp_mode == 'md':        
        print('Method currently not available')
        
    #  add later, since Kim said it could be useful
    elif disp_mode == 'random_normalcoord':
        print('Method currently not available')
        
        
    #  useful for the scan thing to see which modes do contribute
    elif disp_mode == 'pes':
        
        #######################################################################
        
        
        
        #ALM.read_evec_file(evec_file)
        ALM.read_evec_gamma(evec_file, flip=flip)
        
        
        #  this will always target the gamma mode at the beginning of the list (needs to be changed in case one wants to target another k-point or if Gamma is not the first k-point on the list)
        target_q = 0
        
        target_mode = pyro_mode-1

        Qmin = qrange[0]
        Qmax = qrange[1]
        
        

        Qlist = np.array(np.linspace(float(Qmin), float(Qmax), ndata))  # in units of u^{1/2} Angstrom
        
        # if self._verbosity > 0:
        #     print(" Displacement mode              : Along specific normal coordinate\n")
        #     print(" %d configurations are generated from\n"
        #           " the original supercell structure" % number_of_displacements)
        #     print(" The normal coordinate of the following phonon mode is excited:\n\n"
        #           " xk = %f %f %f" %
        #           (self._qpoints[target_q, 0],
        #            self._qpoints[target_q, 1],
        #            self._qpoints[target_q, 2]))
        #     print(" branch : %d" % (target_mode + 1))
        #     print("")


        use_imaginary_part = False
        
        
        tol_zero = 1.0e-3

        ndata = len(Qlist)
        xq_tmp = ALM._qpoints[target_q, :]

        is_commensurate = False
        
        for xq in ALM._commensurate_qpoints:
            xtmp = (xq - xq_tmp) % 1.0
            if math.sqrt(np.dot(xtmp, xtmp)) < tol_zero:
                is_commensurate = True
                break
        if not is_commensurate:
            raise RuntimeWarning("The q point specified by --pes is "
                                 "not commensurate with the supercell.")

        disp = np.zeros((ALM.nat, 3, ndata))

        for iat in range(ALM.nat):
            xshift = ALM._mapping_shift[iat, :]
            jat = ALM._mapping_s2p[iat]
            phase_base = 2.0 * math.pi * np.dot(xq_tmp, xshift)
            cexp_phase = cmath.exp(1.0j * phase_base)

            if use_imaginary_part:
                for icrd in range(3):
                    disp[iat, icrd, :] += Qlist[:] * \
                                          (ALM._evec[target_q, target_mode, 3 * jat + icrd] * cexp_phase).imag
            else:
                for icrd in range(3):
                    disp[iat, icrd, :] += Qlist[:] * \
                                          (ALM._evec[target_q, target_mode, 3 * jat + icrd] * cexp_phase).real

        factor = np.zeros(ALM.nat)
        
        
        for iat in range(ALM.nat):
            # factor[iat] = 1.0 / math.sqrt(self._mass[self._supercell.atomic_kinds[iat]] * float(nq))
            factor[iat] = 1.0 / math.sqrt(ALM._mass[ALM.kd[iat]])


        for idata in range(ndata):
            for i in range(3):
                disp[:, i, idata] = factor[:] * disp[:, i, idata]
                
                # print('first disp: ', disp[:, i, idata])

            # Transform to the fractional coordinate
            # The lattice vector is defined in units of Angstrom, so this operation is fine.
            for iat in range(ALM.nat):
                disp[iat, :, idata] = np.dot(disp[iat, :, idata],
                                             ALM.inverse_lattice_vector)#.T) # !! TRANSPOSE !!
                # print('sce disp: ', disp[iat, :, idata])
                # print(ALM.inverse_lattice_vector.T, '\n')

        if use_imaginary_part:
            if np.max(np.abs(disp)) < 1.0e-10:
                print("Warning. All displacements are zero for the phonon modes specified by --pes.")

        
        disp_list = []


        for i in range(ndata):
            disp_list.append(disp[:, :, i])

        
        for dispNum, disp in enumerate(disp_list):
            cur_disp_x_frac = []

            
            for coordNum, coord in enumerate(disp):
                
                cur_disp_x_frac.append(np.add(ALM.x_fractional[coordNum], coord))
                

            newAseDisp = Atoms(numbers=ALM.numbers, scaled_positions=cur_disp_x_frac, cell=ALM.lattice_vector, pbc=True)
            
            displacements.append(newAseDisp)
            
        
        #######################################################################
        
        
    

    elif disp_mode == 'random_normalcoord_pyro':
        if temp == None or evec_file == None or primitive_aseAtom == None:
            raise ValueError('Please set options for temp (temperature), evec_file (file containing the eigenvectors) and primitive_aseAtom (ase Atoms() object containing the primitive cell)')
        
        
        ALM.read_evec_gamma(evec_file, flip=flip)
        #ALM.read_evec_file(evec_file)
                                
        if altALM != None:
            ALM2 = Alamode(altALM)

        
        nq = 1

        Q_R = np.zeros((nq, ALM._nmode, ndata))
        Q_I = np.zeros((nq, ALM._nmode, ndata))

        # get sigma in units of bohr*amu_ry^(1/2)
        sigma = ALM._get_gaussian_sigma(temp, ignore_imag)

        #Gamma point iq=0, imode given as input
        
        #  -> since iq is only used for the gamma probably no one cares about other k-points; see if it is possible to just skip all the other kpoints, to make the function faster...

        iq = 0
        imode = int(pyro_mode - 1)
        

        if sigma[iq, imode] < 1.0e-10:
            Q_R[iq, imode, :] = 0.0
            Q_I[iq, imode, :] = 0.0
        else:
            Q_R[iq, imode, :] = sigma[iq, imode]
            Q_I[iq, imode, :] = sigma[iq, imode]

        disp_pyro = np.zeros((ALM.nat, 3, ndata))

        for iat in range(ALM.nat):
            xshift = ALM._mapping_shift[iat]
            jat = ALM._mapping_s2p[iat]

            xq_tmp = ALM._qpoints[iq, :]
            phase_real = math.cos(2.0 * math.pi * np.dot(xq_tmp, xshift))


            for icrd in range(3):
                disp_pyro[iat, icrd, :] += Q_R[iq, imode, :] * \
                                      ALM._evec[iq, imode, 3 * jat + icrd].real \
                                      * phase_real

            xq_tmp = ALM._qpoints[iq, :]
            phase = cmath.exp(complex(0.0, 1.0) * 2.0 *
                              math.pi * np.dot(xq_tmp, xshift))

            for icrd in range(3):
                ctmp = ALM._evec[iq, imode, 3 * jat + icrd] * phase
                disp_pyro[iat, icrd, :] += math.sqrt(2.0) * (
                        Q_R[iq, imode, :] * ctmp.real - Q_I[iq, imode, :] * ctmp.imag)

        factor = np.zeros(ALM.nat)
        for iat in range(ALM.nat):
            factor[iat] = 1.0 / math.sqrt(ALM._mass[ALM.kd[iat]]
                                          * ALM._AMU_RYD * float(nq))

        for idata in range(ndata):
            for i in range(3):
                # convert the unit of disp from bohr*amu_ry^(1/2) to angstrom
                disp_pyro[:, i, idata] = factor[:] * disp_pyro[:, i, idata] * (1/ALM._legth_conversion_factor)#ALM._BOHR_TO_ANGSTROM
                #disp[:, i, idata] = factor[:] * disp[:, i, idata]






            #Transform to the fractional coordinate -> double check whether this line needs to be the normal vector or really the transposed version...
            for iat in range(ALM.nat):
                
                disp_pyro[iat, :, idata] = np.dot(disp_pyro[iat, :, idata], ALM.inverse_lattice_vector)#.T) ##!!transpose!!

        # print(disp_pyro)
        
        # print('**********************')
        # for dis in disp_pyro:
        #     print(dis[2])
        #     print()
        # #print(disp_pyro[::][2][::])
        # print(np.shape(disp_pyro))
        # # print(1/ALM._legth_conversion_factor)
        # print('**********************')

        
        cur_disp_x_frac = []
        for coordNum, coords in enumerate(disp_pyro):
            coord_np = np.reshape(coords, (3))
            # print(ALM.x_fractional[coordNum])
            # print(coord_np)   
            # print(np.add(ALM.x_fractional[coordNum], coord_np))
            # print()
            
            cur_disp_x_frac.append(np.add(ALM.x_fractional[coordNum], coord_np))
        
        
        disp_frac = np.reshape(cur_disp_x_frac, (ALM.nat, 3))
    
        newAseDisp = Atoms(numbers=ALM.numbers, scaled_positions=disp_frac, cell=ALM.lattice_vector, pbc=True)
        
    
        displacements.append(newAseDisp)

    else:
        raise ValueError('disp_mode has to be one of ["pattern_file", "random", "pes", "random_normalcoord", "random_normalcoord_pyro"]')
    
    
    return displacements



def extract(aseAtom, dispAtoms, file_name='extractFile', get='disp-force', in_unit_energy='Hartree', in_unit_legth='Angström', offset=None, emin=None, emax=None, CRYSTAL=False, out_files=None, **kwargs):
    
    ALM = Alamode(aseAtom)
    ALM.set_energy_conversion_factor(in_unit_energy)
    ALM.set_length_conversion_factor(in_unit_legth)
    
    if CRYSTAL == False:
    
        ALM._force_conversion_factor = ALM._energy_conversion_factor / ALM._legth_conversion_factor
    
    else:
        
        ALM._force_conversion_factor = ALM._energy_conversion_factor

    #  This makes sure that the file is created as empty file, later on it is appended displacement by displacement
    file = open(file_name, 'w')
    file.close()
    
    
    #  Check what should be printed (-> add new stuff if needed)
    print_disp = False
    print_force = False

    if get.lower() == "disp-force":
        print_disp = True
        print_force = True
    elif get.lower() == "disp":
        print_disp = True
    elif get.lower() == "force":
        print_force = True
    else:
        raise ValueError(
            "Error: Please specify which quantity to extract by the 'get' option.")


    if print_disp == True or print_force == True:
        x0 = np.round(ALM.x_fractional, 8)  #  just rounds the fractional coordinates... -> just use x_frac
        #lavec_transpose = ALM.lattice_vector#*ALM._legth_conversion_factor  # get the lavec in bohr
        
        vec_refold = np.vectorize(Alamode.refold)  #  results in some shift of the origin of the displacements from the corner of the cell to the centre of the cell; weird but maybe the C++ code needs that...


    if offset is None:
        displacement_offset = np.zeros((ALM.nat, 3))
        force_offset = np.zeros((ALM.nat, 3))
        epot_offset = 0.0
    else:
        ALM_offset = Alamode(offset)
        x_offset = ALM_offset.x_fractional
        
        displacement_offset = x_offset - x0
        force_offset = ALM_offset._get_atomic_forces(offset, **kwargs)
        epot_offset = ALM_offset._get_total_energy(offset, **kwargs)
            
#  Probably rework that part by not calling the functions to get energies/forces, just directly call ase function (and then reshape properly)      
#  double check if reshaping every time is necessary... cause ase should have proper shapes...
###############################################################################
    counter = 1
    for disp_num, disp in enumerate(dispAtoms):
        
        if disp != 0:
        
            #  create an Alaomde object to load coordinates etc. in the right sorted way
            ALM_disp = Alamode(disp)
            
            x = ALM_disp.x_fractional
            
            if CRYSTAL == False:
                force = ALM_disp._get_atomic_forces(disp, **kwargs)  #  get forces as [len(x), 3] matrix
                epot = ALM_disp._get_total_energy(disp, **kwargs)  #  get them as a single value
    
            else:
                force = ALM_disp._get_atomic_forces(out_files[disp_num], CRYSTAL=True, **kwargs)  #  get forces as [len(x), 3] matrix
                epot = ALM_disp._get_total_energy(out_files[disp_num], CRYSTAL=True, **kwargs)  #  get them as a single value
                            
            
    
    
            if x is None or force is None or epot is None:
                continue
            
            
            # shape_of_force = np.shape(force)  #  (40,3)
                    
            # num_elements_in_force = shape_of_force[0]*shape_of_force[1]  #  120
            
            # num_data_force = num_elements_in_force // 3  #  1
                    
            # force = np.reshape(force, (num_data_force, ALM_disp.nat, 3))
            
            # num_data_disp = np.shape(x)  #  (40,3)
                            
    
            # if num_data_disp[0] != num_data_force and print_disp and print_force:
            #     print(
            #         "Error: The number of entries of displacement and force is inconsistent.")
            #     print("Ndata disp : %d, Ndata force : %d" %
            #           (num_data_disp[0], num_data_force))
            #     exit(1)
    
            # ndata_energy = len(epot)
            
            # if ndata_energy != num_data_disp:
    
                
                # raise RuntimeError("The numbers of displacement and energy entries are different.")
    
            
            epot -= epot_offset
            
    
            if emin is not None:
                if emin > epot:
                    continue
    
            if emax is not None:
                if emax < epot:
                    continue
    
            if print_disp:
                
                displ = x - x0 - displacement_offset
                
                # print(displ)
                
                displ = np.dot(vec_refold(displ), ALM.lattice_vector)#.transpose())
    
                displ *= ALM._legth_conversion_factor
    
            if print_force:
                f = force - force_offset
                f *= ALM._force_conversion_factor  
                
                # print(ALM._force_conversion_factor)
            
            
            with open(file_name, 'a') as file:
                
                file.write("# Displacement: %s,  E_pot (%s): %s\n" %
                      (counter, in_unit_energy, epot))
    
    
                if print_disp and print_force:
                    for i in range(ALM_disp.nat):
                        file.write("%15.7F %15.7F %15.7F %20.8E %15.8E %15.8E\n"%(displ[i, 0],                                                                         displ[i, 1], displ[i, 2], f[i, 0], f[i, 1], f[i, 2]))
                        
                elif print_disp:
                    for i in range(ALM_disp.nat):
                        file.write("%15.7F %15.7F %15.7F\n" % (displ[i, 0], displ[i, 1], displ[i, 2]))
                        
                elif print_force:
                    for i in range(ALM_disp.nat):
                        file.write("%15.8E %15.8E %15.8E\n" % (f[i, 0], f[i, 1], f[i, 2]))
                
                counter += 1
    return





#  command_list is a dictionary of different keywords and their values used in the ALM input -> easiest to handle as **kwargs
def write_ALM(aseAtom, file_name='ALMinput.in', command_list = None, cutoffs=None, get_XML_files=[False, False, False], **kwargs):
    
    if file_name != '' and file_name.endswith('.in') == False:
        file_name += '.in'
        
    if type(command_list) is not dict:
        print('Loading variables from a template file.\n')        
        ALM_keys = get_keywords_from_template(command_list)

    else:
        ALM_keys = command_list
        
    if any(get_XML_files) == True:
        
        xml_keys = ['FC2XML', 'FC3XML', 'FC4XML']
                
        fc2xml, fc3xml, fc4xml = get_fcsXML(fcs2=get_XML_files[0], fcs3=get_XML_files[1], fcs4=get_XML_files[2])
        
        for xmlNum, xml in enumerate([fc2xml, fc3xml, fc4xml]):
            if xml != '':
                ALM_keys[xml_keys[xmlNum]] = xml
                
    if cutoffs != None:
        cutoffs_harmonic = cutoffs[0]
        cutoffs_anharmonic = cutoffs[1]
    else:
        cutoffs_harmonic = []
        cutoffs_anharmonic = []
        
    
    general_In = ALMgeneral(aseAtom, **ALM_keys)
    interaction_In = interaction(**ALM_keys)
    optimize_In = optimize(**ALM_keys)
    
    #  This block checks if attributes are set for the different blocks, so that only the blocks with 
    blocks = ['generalInput', 'interactionInput', 'optimizeInput']
    in_list = [general_In, interaction_In, optimize_In]
    set_check = [False for i in blocks]
    
    for obj_num, obj in enumerate(in_list):
        cur_block = obj.get_all_att()
        for key in cur_block:
            if cur_block.get(key) != None:
                set_check[obj_num] = True
            if set_check[obj_num] == True:
                ALM_keys[blocks[obj_num]] = obj
                continue
        
    write_ALM_input(aseAtom, Cutoffs_harm=cutoffs_harmonic, Cutoffs_anharm=cutoffs_anharmonic, fileName=file_name, **ALM_keys, **kwargs)
         
    return


def write_ANPHON(aseAtom, file_name='ANPHONinput.in', command_list = None, get_XML_files=False, strainInput=None, displaceInput=None, DISPMODE=None, source='ase', contPath=False, KPMODE=1, **kwargs):
    
    if file_name != '' and file_name.endswith('.in') == False:
        file_name += '.in'
        
    if type(command_list) is not dict:
        print('Loading variables from a template file.\n')        
        ANPHON_keys = get_keywords_from_template(command_list)

    else:
        ANPHON_keys = command_list
        
    if get_XML_files == True:
        
        if ANPHON_keys['MODE'] == 'phonons' or ANPHON_keys['MODE'] == ['phonons']:
            fc2xml, fc3xml, fc4xml = get_fcsXML(fcs2=True)
            ANPHON_keys['FCSXML'] = fc2xml
            
        elif ANPHON_keys['MODE'] == 'SCPH' or ANPHON_keys['MODE'] == ['SCPH']:
            fc2xml, fc3xml, fc4xml = get_fcsXML(fcs4=True)
            ANPHON_keys['FCSXML'] = fc4xml
        
    general_in = ANPHONgeneral(aseAtom, **ANPHON_keys)
    
    if ANPHON_keys['MODE'] == ['SCPH']:
        scph_in = scph(**ANPHON_keys)
    else:
        scph_in = scph()
    
    if ANPHON_keys['MODE'] == ['QHA']:
        qha_in = qha(**ANPHON_keys)
    else:
        qha_in = qha()
    
    relax_in = relax(**ANPHON_keys)
    analysis_in = analysis(**ANPHON_keys)
    
    blocks = ['generalInput', 'scphInput', 'qhaInput', 'relaxInput', 'analysisInput']
    in_list = [general_in, scph_in, qha_in, relax_in, analysis_in]
    set_check = [False for i in blocks]
    
    for obj_num, obj in enumerate(in_list):
        cur_block = obj.get_all_att()
        for key in cur_block:
            if cur_block.get(key) != None:
                set_check[obj_num] = True
            if set_check[obj_num] == True:
                ANPHON_keys[blocks[obj_num]] = obj
                continue
    
    write_ANPHON_input(aseAtom, strainInput=strainInput, displaceInput=displaceInput, DISPMODE=DISPMODE,  fileName=file_name, source=source, KPMODE=KPMODE, contPath=contPath, **ANPHON_keys, **kwargs)


    return






def get_polarization_direction_from_sg(space_group):
    
    z_pol_num = [6,7,8,9,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,75,76,77,78,79,80,99,100,101,102,103,104,105,106,107,108,109,110,143,144,145,146,156,157,158,159,160,161,168,169,170,171,172,173,183,184,185,186]
    
    # hex: 168,169,170,171,172,173,183,184,185,186
    # trigonal: 143,144,145,146 (rhomb),156,157,158,159,160 (rhomb),161 (rhomb)
    # tetragonal: 75,76,77,78,79,80,99,100,101,102,103,104,105,106,107,108,109,110 (a,a,c)
    # orthorhombic: 25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46 (a,b,c)
    # monoclinic: 3,4,5,6,7,8,9 (a,b,c,alpha)
    # triclinic: 1
    
    
    if space_group == 1:
        pol_dir = ['x', 'y', 'z']
        # print('Space group %s polarization along %s possible for.'%(space_group, pol_dir))
        
    elif space_group in [3,4,5]:
        pol_dir = ['x', 'y']
        # print('Space group %s polarization along %s possible for.'%(space_group, pol_dir))
        
    elif space_group in z_pol_num:
        pol_dir = ['z']
        # print('Space group %s polarization along %s possible for.'%(space_group, pol_dir))
    else:
        pol_dir = []
        # print('Space group %s has no ploarization.'%space_group)
    
    # print()
    
    return pol_dir





 







def add_tags_to_crystal(aseAtom, treshold=36):
    
    tag_list = []
    
    for number in aseAtom.numbers:
        if number > treshold:
            tag_list.append(200)
        else:
            tag_list.append(0)
    
    aseAtom.set_tags(tag_list)
    
    return


def select_displacements(disp_list, selected_list):
    
    pos_list = np.zeros(len(disp_list))
    
    
    for selected in selected_list:
        print(selected)
        
        if selected.find('-') != -1:
            parts = selected.split('-')

            
            num_range = np.arange(int(parts[0]), int(parts[1])+1, 1)

        else:
            num_range = [int(selected)]
        
        for num in num_range:
            
            pos_list[num-1] = 1
            
            
    print(pos_list)
    print()
    
    new_disps_list = []
    
    
    for disp_num, disp in enumerate(disp_list):
        if pos_list[disp_num] == 0:
            new_disps_list.append(0)
        else:
            new_disps_list.append(disp)

    return new_disps_list
    


    
    
import os

def get_evec_prop(cell, tempList=[100]):
    for tempNum, temp in enumerate(tempList):
        ALM = Alamode(cell)
        
        for file in os.listdir(os.getcwd()):
        
            if file.find('_%s.'%temp)!=-1 and file.find('.evec') !=-1:
                 
                ALM.read_evec_gamma(file)
                
                if temp == tempList[0]:
                    
                    ALM_ini = Alamode(cell)
                    ALM_ini.read_evec_gamma(file)
                    
                    evecMatchList = np.zeros((len(ALM_ini._evec[0][::]), len(tempList)))
                    
                    for evecNum, evec in enumerate(ALM_ini._evec[0][::]):
                        evecMatchList[evecNum][tempNum] = evecNum+1
                    
                    
                else:
                    
                    ALM_cur = Alamode(cell)
                    ALM_cur.read_evec_gamma(file)
                    
                    for evecIniNum, evecIni in enumerate(ALM_ini._evec[0][::]):
                        
                        curBestMatchNum = 50000000.0
                    
                        for evecNum, evec in enumerate(ALM._evec[0][::]):
                            
                            evecDiff = ALM_ini._evec[0][evecIniNum].real - ALM_cur._evec[0][evecNum].real
                            
                            evecDiff_sum = np.sum([np.abs(x) for x in evecDiff])
                            
                            evecSumm = ALM_ini._evec[0][evecIniNum].real + ALM_cur._evec[0][evecNum].real
                            evecSumm_sum = np.sum([np.abs(x) for x in evecSumm])
                            
                            if curBestMatchNum > np.abs(evecDiff_sum):
                                curBestMatchNum = np.abs(evecDiff_sum)
                                evecMatchList[evecIniNum][tempNum] = evecNum+1
                                
                            if curBestMatchNum > np.abs(evecSumm_sum):
                                curBestMatchNum = np.abs(evecSumm_sum)
                                
                                evecMatchList[evecIniNum][tempNum] = evecNum+1
                            
        
    with open('evecProp.txt', 'w') as infile:
        infile.write('Propagation of the eigenvectors: \n\n')
        
        for temp in tempList:
            infile.write(f'{temp:>6.0f}K ')
        
        infile.write('\n')

        for modeNum, mode in enumerate(evecMatchList):
            for modeT in mode:
                infile.write(f'{modeT:>7.0f} ')
            infile.write('\n')
            
            
    return

def get_pol_dir(cell, tempList=[100]):

    for tempNum, temp in enumerate(tempList):
        
        ALM = Alamode(cell)

        for file in os.listdir(os.getcwd()):
        
            if file.find('_%s.'%temp)!=-1 and file.find('.evec') !=-1:
                
                ALM.read_evec_gamma(file)
                
                if temp == tempList[0]:
                    
                    evecPolList = np.zeros((len(ALM._evec[0][::]), len(tempList)))
                    evecPolDirList = [[] for x in ALM._evec[0][::]]
                    
        
        for evecNum, evec in enumerate(ALM._evec[0][::]):
            #dispVec = [0,0,0]
            
            evecShape = np.reshape(evec, (len(cell.numbers), 3))
            
            dispX = 0
            dispY = 0
            dispZ = 0
            
            for pos in evecShape:
                dispX += pos[0].real
                dispY += pos[1].real
                dispZ += pos[2].real
            
            if np.abs(dispX) > max([np.abs(dispY), np.abs(dispZ)]):
                if dispX > 0:
                    evecPolDirList[evecNum].append('-x') 
                    evecPolList[evecNum][tempNum] = -dispX
                elif dispX < 0:
                    evecPolDirList[evecNum].append('x') 
                    evecPolList[evecNum][tempNum] = dispX
            
            elif np.abs(dispY) > max([np.abs(dispX), np.abs(dispZ)]):
                if dispY > 0:
                    evecPolDirList[evecNum].append('-y') 
                    evecPolList[evecNum][tempNum] = -dispY
                elif dispY < 0:
                    evecPolDirList[evecNum].append('y') 
                    evecPolList[evecNum][tempNum] = dispY
                    
            elif np.abs(dispZ) > max([np.abs(dispX), np.abs(dispY)]):
                if dispZ > 0:
                    evecPolDirList[evecNum].append('-z') 
                    evecPolList[evecNum][tempNum] = -dispZ
                elif dispZ < 0:
                    evecPolDirList[evecNum].append('z') 
                    evecPolList[evecNum][tempNum] = dispZ
            else:
                evecPolDirList[tempNum].append('prob 0') 
                evecPolList[evecNum][tempNum] = 0


    with open('evecPol.txt', 'w') as infile:
        infile.write('Possible polarisation direction based on displacement direction (direction/value): \n\n')
        
        infile.write('Mode ')
        
        for temp in tempList:
            infile.write(f'{temp:>13.0f}K')
        
        infile.write('\n')

        for polDirNum, polDir in enumerate(evecPolDirList):
            infile.write(f'{polDirNum+1:<5.0f} ')
            
            for polNum, pol in enumerate(polDir):
                
                if pol in ['z', '-z'] and np.abs(evecPolList[polDirNum][polNum]) > 1e-6:
                    infile.write(f'{pol:>6} ')
                    infile.write(f'{evecPolList[polDirNum][polNum]:>6.2f} ')
                elif pol in ['x', '-x'] and np.abs(evecPolList[polDirNum][polNum]) > 1e-6:
                    infile.write(f'{pol:>6} ')
                    infile.write(f'{evecPolList[polDirNum][polNum]:>6.2f} ')
                elif pol in ['y', '-y'] and np.abs(evecPolList[polDirNum][polNum]) > 1e-6:
                    infile.write(f'{pol:>6} ')
                    infile.write(f'{evecPolList[polDirNum][polNum]:>6.2f} ')
                    
                else:
                    s = ''
                    infile.write(f'{s:>14}')
            
            infile.write('\n')

    # list for the displacements that should be done since they probably are A1 modes in the form of:
    # [Temp, mode, flip]
    # [int, int, True/False]
    pyroDispList = []

    for polDirNum, polDir in enumerate(evecPolDirList):
        for polNum, pol in enumerate(polDir):

            if pol == 'z' and np.abs(evecPolList[polDirNum][polNum]) > 1e-6:
                pyroDispList.append([tempList[polNum], polDirNum+1, False])

            elif pol== '-z' and np.abs(evecPolList[polDirNum][polNum]) > 1e-6:
                pyroDispList.append([tempList[polNum], polDirNum+1, True])


    with open('dispPattern.txt', 'w') as infile:
        infile.write('Suggestion for displacements to be made for the temperature dependent pyro displacements: \n\n')
        
        for temp in tempList:
            infile.write(f'{temp:>6.0f}K')
        
        infile.write('\n')

        for polDirNum, polDir in enumerate(evecPolDirList):
            infile.write(f'{polDirNum+1:<5.0f} ')
            
            for polNum, pol in enumerate(polDir):
                if pol == 'z' and np.abs(evecPolList[polDirNum][polNum]) > 1e-6:
                    infile.write(f'{polDirNum+1:<6.0f} ')
                elif pol == '-z' and np.abs(evecPolList[polDirNum][polNum]) > 1e-6:
                    infile.write(f'{(polDirNum+1)*(-1):<6.0f} ')
                else:
                    s = ''
                    infile.write(f'{s:>6} ')
            
            infile.write('\n')

        
        
        return pyroDispList




if __name__ == '__main__':
    
    print(r'/\___/\   ')
    print(r':O: :O:   Testing Alamode class')
    print(r' ¨   ¨')
    print()
    
    
    
    print('Test')
    
    def get_ase_from_crystal_out(crystal_out_file, set_mode='prim', dump_object=True):
        
        with open(crystal_out_file, 'r') as in_file:
            found_converged = False
            found_prim = False
            found_prim_matrix = False
            found_conv = False
            counter_prim = 0
            counter_conv = 0
            numbers_prim = []
            positions_prim = []
            matrix_prim = []
            
            cell_prim = []
            cell_conv = []
            
            numbers_conv = []
            positions_conv = []


            for line in in_file:
                if line.find('OPT END - CONVERGED') != -1:
                    found_converged = True
                
                if found_prim == True and set_mode == 'prim':
                    if counter_prim == 1:
                        cell_prim = [float(x) for x in line.split()]
                    elif counter_prim > 5:
                        if line.split() != []:
                            number = int(line.split()[2])
                            if number > 55:                            
                                numbers_prim.append(number-200)
                            else:
                                numbers_prim.append(number)
                            
                            positions_prim.append([float(line.split()[4]), float(line.split()[5]), float(line.split()[6])])
                            
                            # print(line.split())
                        else:
                            found_prim = False
                        
                    counter_prim += 1
                        
                if found_prim_matrix == True and set_mode == 'prim':
                    if counter_prim > 0:
                        if line.split() != []:
                            matrix_prim.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
                        else:
                            found_prim_matrix = False
                    
                    # print(line.split())
                        
                    counter_prim += 1
                    
                if found_conv == True and set_mode == 'conv':
                    if counter_conv == 1:

                        cell_conv = [float(x) for x in line.split()]
                    elif counter_conv > 5:
                        if line.split() != []:                            
                            number = int(line.split()[2])
                            if number > 55:                            
                                numbers_conv.append(number-200)
                            else:
                                numbers_conv.append(number)
                            
                            positions_conv.append([float(line.split()[4]), float(line.split()[5]), float(line.split()[6])])
                            
                            
                            
                        else:
                            found_conv = False
                    
                    # print(line.split())
                    
                    
                
                    counter_conv += 1
                
                
                if found_converged == True:
                    if line.find('PRIMITIVE CELL - CENTRING') != -1:
                        found_prim = True
                    elif line.find('DIRECT LATTICE VECTORS') != -1:
                        found_prim_matrix = True
                        counter_prim = 0
                        # print(line.split())
                        #print()
                    elif line.find('CRYSTALLOGRAPHIC CELL (') != -1:
                        found_conv = True
                        counter_conv = 0

                #print(line)
            
        
        
        # print()
        # print(numbers_prim)
        
        # print()
        # print(positions_prim)
        
        # print()
        # print(matrix_prim)
        
        # print()
        # print(cell_conv)
        
#        TBC: use variables to build an ase Atom() object (either with matrix/scaled positions for prim cell or cell parametrs/scaled positions for conventional cell)

        if set_mode == 'prim':
            asePrim = Atoms(scaled_positions=positions_prim, numbers=numbers_prim, cell=cell_prim)
            return asePrim
            
        elif set_mode == 'conv':
            aseConv = Atoms(scaled_positions=positions_conv, numbers=numbers_conv, cell=cell_conv)
            return aseConv
        
        else:
            print('set_mode not supported')
            
    





























