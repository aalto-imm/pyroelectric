"""
This script should incorporate function to
- read forces from a crystal output
- read the total energy from a crystal output
- convert ase Atoms() object to proper .d12 files (with proper symmetry or in P1 (-> for the supercells))

maybe later on additional stuff to get the energies and born stuff (that's currently not really implemented in extract())
"""

def get_energy_CRYSTAL(CRYSTAL_outfile):    
    with open(CRYSTAL_outfile, 'r') as file:
        for line in file:
            if line.find('SCF ENDED - CONVERGENCE ON ENERGY') != -1:    
                epot = float(line.split()[8])
                break
                     
    return epot



def get_forces_CRYSTAL(CRYSTAL_outfile):
    
    found_forces = False
    forces = []
    
    with open(CRYSTAL_outfile, 'r') as file:
        for line in file:
            
            if found_forces == True:
                if line.split() == []:
                    break
                elif line.split()[0] == 'ATOM':
                    continue

                forces.append([float(line.split()[2]), float(line.split()[3]), float(line.split()[4])])
            
            
            if line.find('CARTESIAN FORCES IN HARTREE/BOHR (ANALYTICAL)') != -1:    
                found_forces = True
        
    return forces














if __name__ == '__main__':
    
    
    print(r'/\___/\   ')
    print(r':O: :O:   Testing ase CRYSTAL interface')
    print(r' ¨   ¨')
    print()
    
    
    energy = get_energy_CRYSTAL('disp1.out')
    
    print(energy)
    print()
    
    forces = get_forces_CRYSTAL('disp1.out')
    
    print(forces)
    print(len(forces))
    print()
    
    
    



