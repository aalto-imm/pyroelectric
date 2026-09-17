"""
This script contains the relevant functions for plotting the phonon bands and the bands from the scph calculation

"""

import sys
import matplotlib.pyplot as plt
import numpy as np

#  function to read in the output .bands for phonons with two options for the filetype (phonons and scph)

#  add option to get one (or more) specific temperatures

def read_bands_data(bands_file, get300=False, getTemp=[]):
    
    if bands_file.endswith('.bands') == True:
        file_type = 'phonons'
    elif bands_file.endswith('.scph_bands') == True:
        file_type = 'scph'
    else:
        print('unknown file type...')
        sys.exit()
    
    
    with open(bands_file, 'r') as infile:
        k_point = infile.readline().split()[1::]
        
        k_point_x = [float(x) for x in infile.readline().split()[1::]]

        if file_type == 'phonons':
                
            infile.readline()
            
            bands = [[[float(x)] for x in infile.readline().split()]]
                        
            for line in infile:
                
                for bandNum, band in enumerate(line.split()):
                    
                    bands[0][bandNum].append(float(band))
                    
        elif file_type == 'scph':
            if get300 == True:
                
                infile.readline()
                #firstLine = infile.readline()
                bands = [[[0.0] for x in infile.readline().split()[1::]]]
                
                
                for line in infile:
                    
                    
                    if line.split() != [] and float(line.split()[0]) == 300.:
                        #print(np.shape(bands))
                        
                        for bandNum, band in enumerate(line.split()[1::]):
                            
                            bands[0][bandNum].append(float(band))

                        

            elif getTemp != []:
                
                infile.readline()

                bands = [[[0.0] for x in infile.readline().split()[1::]] for temp in getTemp]
                
                for line in infile:
                    
                    
                    if line.split() != [] and float(line.split()[0]) in getTemp:                        
                        #print(np.shape(bands))
                        tempIndex = getTemp.index(float(line.split()[0]))
                        
                        for bandNum, band in enumerate(line.split()[1::]):
                            
                            bands[tempIndex][bandNum].append(float(band))
                
                
                
                
            else:
                infile.readline()
                
                first_line = infile.readline().split()
                
                curTemp = float(first_line[0])
                tempIndex = 0
                
                bands = [[[float(x)] for x in first_line[1::]]]
                
                for line in infile:
                    
                    cur_line = line.split()
                    
                    if cur_line == []:
                        continue
                    elif float(cur_line[0]) == curTemp:
                        for bandNum, band in enumerate(cur_line[1::]):
                            bands[tempIndex][bandNum].append(float(band))
                    elif float(cur_line[0]) != curTemp:
                        
                        curTemp = float(cur_line[0])
                        tempIndex += 1
                        bands.append([[float(x)] for x in cur_line[1::]])
        
        infile.close() 
    
    # print()              
    # print(bands)
    #print(np.shape(bands))
    
    return bands, k_point, k_point_x


def plot_phonon_bands(bands_file, labels=[], harmData='', save_path='plot.png', get300=False, getTemp=[]):
    
    bands, k_point, k_point_x = read_bands_data(bands_file, get300=get300, getTemp=getTemp)
        
    fig = plt.figure()
    
    ax = fig.add_subplot()
    
    # expand like needed to handle greek letters/indices
    for labelNum, label in enumerate(k_point):
        if label == 'G' or label.lower() == 'gamma':
            k_point[labelNum] = r'$\Gamma$'
        if label == 'C_0':
            k_point[labelNum] = r'C$_0$'
        if label == 'SIGMA_0':
            k_point[labelNum] = r'$\Sigma_0$'
        if label == 'H_0':
            k_point[labelNum] = r'H$_0$'
        if label == 'S_0':
            k_point[labelNum] = r'S$_0$'
        if label == 'A_0':
            k_point[labelNum] = r'A$_0$'
        if label == 'E_0':
            k_point[labelNum] = r'E$_0$'

    
    #  sort out breaks in k-path
    k_point_cor = [k_point[0]]
    k_point_x_cor = [k_point_x[0]]
    
    for kxNum, kx in enumerate(k_point_x):

        if kxNum > 0:
            if kx == k_point_x[kxNum-1]:
                k_point_cor[-1] += '|' + k_point[kxNum]
            else:
                k_point_cor.append(k_point[kxNum])
                k_point_x_cor.append(kx)
    
    #  plotting vertical lines for k-points
    for tick in k_point_x_cor:
        ax.axvline(tick, color='grey', lw=0.5)
    
    #  color map for plotting scph convergence plots
    cmap = plt.get_cmap('viridis')
    if len(bands) >1:
        colors = [cmap(value) for value in np.linspace(0, .9, len(bands))]
    else:
        colors = ['navy']
    
    
    for tempNum, temp in enumerate(bands):
        for band in temp[1::]:
            ax.plot(temp[0], band, color=colors[tempNum], ls='-', lw=0.7)
            


    if harmData != '':
        band_harm, k_point_harm, k_point_x_harm = read_bands_data(harmData)
        
        for tempNum, temp in enumerate(band_harm):
            for band in temp[1::]:
                ax.plot(temp[0], band, color='grey', ls=':', lw=0.7)
        


    ax.hlines(0, min(k_point_x_cor), max(k_point_x_cor), ls='-', color='black', lw=0.5)        
        
        
    if labels != []:
        if len(labels) == 1:
            ax.plot(0,0, marker='', color='navy', ls='-', lw=0.7, label=labels[0])
        
        else:
            
            for labelNum, label in enumerate(labels[:-1:]):
                ax.plot(0,0, marker='', color=colors[labelNum], ls='-', lw=0.7, label=label)
        
        if harmData != '':
            ax.plot(0,0, marker='', color='grey', ls=':', lw=0.7, label=labels[-1])

        
        if getTemp == []:
            plt.legend(loc='upper right')
        else:
            plt.legend(bbox_to_anchor=(1.01, 1.0), loc='upper left')


    for k_pointNum, k_point in enumerate(k_point_cor):
        if k_point.find('|') != -1:
            ax.axvline(k_point_x_cor[k_pointNum], color='grey', lw=0.7)
    
    
    ax.set_ylabel(r'Frequency (cm$^{-1}$)')
    
    ax.set_xticks(k_point_x_cor)
    ax.set_xticklabels(k_point_cor)
    
    ax.set_xlim(min(k_point_x_cor), max(k_point_x_cor))
    
    
    plt.tight_layout()
    
    plt.savefig(save_path, dpi=500)
    
    plt.show()
    
    
    return

import os


def plot_pyro_constants(pyro_modes, temperatures, folder='', pic_name=''):
    
    fig = plt.figure()
    
    ax = fig.add_subplot()
    
    if folder == '':
        cwd = os.get_cwd()
    else:
        cwd = folder

    all_modes = []

    for mode in pyro_modes:
        cur_mode_Ps = []
        
        for temp in temperatures:
            
            for file in os.listdir(cwd):
                #print(file)
                if file.find(str(temp))!= -1 and file.find('mode'+str(mode)+'_') != -1 and file.find('.out') != -1:
                    cur_out_file = file                        
                    break

            found_polarization_values = False

            with open(cwd + cur_out_file, 'r') as file:
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
        
        ax.plot(temperatures, z_values, marker='x', label=f'mode {mode}')
        
        #ax.annotate(f'mode {mode}', xy=(args.TEMPERATURE[0], z_values[0]), xytext=(args.TEMPERATURE[0]+1, z_values[0]))
        
        all_modes.append(z_values)
        
        sum_p = []
    
    for tempNum, temp in enumerate(all_modes[0][::]):
        
        sum_p.append(all_modes[0][tempNum] + all_modes[1][tempNum] + all_modes[2][tempNum])
        
    ax.plot(temperatures, sum_p, marker='x', color='black', label=r'total (p$^{(1)}$)')
        
        
    
    ax.set_xlabel('temperature (K)', fontsize=12)
    ax.set_ylabel(r'$p$ ($\mu$Cm$^{-2}$K$^{-1}$)', fontsize=12)

    plt.legend()
    
    plt.tight_layout()
    plt.savefig(cwd + f'{pic_name}_pyro.png')

    plt.show()
    
    
    
    
    return





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
    print(r':O: :O:   Script for plotting phonons')
    print(r' ¨   ¨')
    print()
    
    
    
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    