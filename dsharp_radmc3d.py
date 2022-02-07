import numpy as np
import pandas as pd
import os, glob, argparse, warnings, pickle

import dsharp_opac as opacity

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']


def create_tables(a, lam, fm_ice, filename):
    extrapolate_large_grains = True
    
    output_dir = 'data'
    if fm_ice == 0.2:
        name = 'default_opacities'
    elif fm_ice == 0.0:
        name = 'icefree_opacities'
    else:
        warnings.warn('fm_ice should be 0.0 or 0.2')
        name = 'other_opacities'
        
    if extrapolate_large_grains:
        name += '_extrapol'
        
    print('opacity name: {}'.format(name))
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    dc,rho_s = opacity.get_dsharp_mix(fm_ice=fm_ice)
    res_l = opacity.get_smooth_opacities(a, lam, rho_s, dc, smoothing='linear', n_angle=180, extrapolate_large_grains=extrapolate_large_grains)

    with open('{}.pickle'.format(filename), 'wb') as fid:
        pickle.dump(res_l, fid)


def dustkappa(a, amax, lam, fm_ice, filename):
    extrapolate_large_grains = True
    
    if fm_ice == 0.2:
        name = 'default_opacities'
    elif fm_ice == 0.0:
        name = 'icefree_opacities'
    else:
        warnings.warn('fm_ice should be 0.0 or 0.2')
        name = 'other_opacities'
        
    if extrapolate_large_grains:
        name += '_extrapol'
        
    dc,rho_s = opacity.get_dsharp_mix(fm_ice=fm_ice)
    res = opacity.get_smooth_opacities(a*1e-4, lam*1e-4, rho_s, dc, smoothing='linear', n_angle=180, extrapolate_large_grains=extrapolate_large_grains)

    k_abs = res['k_abs']
    k_sca = res['k_sca']
    g = res['g']
    S1 = res['S1']
    S2 = res['S2']

    theta = res['theta']

    s = a**0.5
    s[a > amax*1e-4] = 0
    s = s / s.sum()
    k_abs = (k_abs*s[:,None]).sum(0)
    k_sca = (k_sca*s[:,None]).sum(0)
    g = (g*s[:,None]).sum(0)

    frame = np.stack((lam, k_abs, k_sca, g))


    #---CONVERT INTO RADMC3D READABLE INPUT FILES
    f = open('dustkappa_dsharp-{}.inp'.format(filename),"w")
    f.write('# Opacity and scattering matrix file for bin{}\n'.format(str(bin+1).zfill(2)))
    f.write('# Please do not forget to cite in your publications the original paper of these optical constant measurements\n')
    f.write('# Made with the DSHARP_OPAC package by Cornelis Dullemond & Til Birnstiel\n')
    f.write('# using either the bhmie.py Mie code of Bohren and Huffman (python version by Cornelis Dullemond,\n')
    f.write('# or a F90 version by Til Birnstiel, both after the original bhmie.f code by Bruce Draine)\n')
    f.write('# Material density = {0:6.3f} g/cm^3\n'.format(rho_s))
    f.write('3\n')  # Format number
    f.write('{0:d}\n'.format(res["lam"].size))
    f.write('\n')

    for ilam in range(res['lam'].size):
        for i in range(0, 4, 1):
            f.write("%.6e  "%(frame[i][ilam]))
        f.write("\n")

    f.close()


def dustkapscatmat(a, amax, lam, fm_ice, filename):
    extrapolate_large_grains = True
    
    if fm_ice == 0.2:
        name = 'default_opacities'
    elif fm_ice == 0.0:
        name = 'icefree_opacities'
    else:
        warnings.warn('fm_ice should be 0.0 or 0.2')
        name = 'other_opacities'
        
    if extrapolate_large_grains:
        name += '_extrapol'
        

    dc,rho_s = opacity.get_dsharp_mix(fm_ice=fm_ice)
    res = opacity.get_smooth_opacities(a*1e-4, lam*1e-4, rho_s, dc, smoothing='linear', n_angle=180, extrapolate_large_grains=extrapolate_large_grains)
    print(res["a"].size)
    k_abs = res['k_abs']
    k_sca = res['k_sca']
    g = res['g']
    S1 = res['S1']
    S2 = res['S2']
    theta = res['theta']

    s = a**0.5
    s[a > amax*1e-4] = 0
    s = s / s.sum()
    k_abs = (k_abs*s[:,None]).sum(0)
    k_sca = (k_sca*s[:,None]).sum(0)
    g = (g*s[:,None]).sum(0)
    print(s)

    m  = 4 * np.pi / 3 * rho_s * res['a']**3
    MM = opacity.calculate_mueller_matrix(res['lam'], m, res['S1'], res['S2'], theta=theta, k_sca=k_sca)

    frame = np.stack((lam, k_abs, k_sca, g))

    #---CONVERT INTO RADMC3D READABLE INPUT FILES
    f = open('dustkapscatmat_dsharp-{}.inp'.format(filename),"w")
    f.write('# Opacity and scattering matrix file for bin{}\n'.format(str(bin+1).zfill(2)))
    f.write('# Please do not forget to cite in your publications theo riginal paper of these optical constant measurements\n')
    f.write('# Made with the DSHARP_OPAC package by Cornelis Dullemond & Til Birnstiel\n')
    f.write('# using either the bhmie.py Mie code of Bohren and Huffman (python version by Cornelis Dullemond,\n')
    f.write('# or a F90 version by Til Birnstiel, both after the original bhmie.f code by Bruce Draine)\n')
    f.write('# Material density = {0:6.3f} g/cm^3\n'.format(rho_s))
    f.write('1\n')  # Format number
    f.write('{0:d}\n'.format(res["lam"].size))
    f.write('{0:d}\n'.format(res["theta"].size))
    f.write('\n')

    for ilam in range(res['lam'].size):
        for i in range(0, 4, 1):
            f.write("%.6e  "%(frame[i][ilam]))
        f.write("\n")
    f.write("\n")

    for iang in range(len(theta)):
        f.write("%.6e  "%(theta[iang]))
        f.write("\n")
    f.write("\n")

    for ilam in range(res['lam'].size):
        for itheta in range(res['theta'].size):
            f.write('%13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n' %
                    (MM['zscat'][s, ilam, itheta, 0], MM['zscat'][s, ilam, itheta, 1],
                        MM['zscat'][s, ilam, itheta, 2], MM['zscat'][s, ilam, itheta, 3],
                        MM['zscat'][s, ilam, itheta, 4], MM['zscat'][s, ilam, itheta, 5]))
        f.write('\n')

    f.close()

def make_figure(lmin, lmax):
    filelist = sorted(glob.glob('dustkappa*'))

    fig = plt.figure(figsize=(9.6, 7.2)) #fig = plt.figure(figsize=(9.6, 7.2))
    ax = fig.add_subplot(111) 
    #ax2 = ax.twinx()
    ax.set_xlabel(r'$\lambda$ [$\mu$m]', fontsize=20)
    #ax.set_ylabel(r'$\omega$', fontsize=20)
    ax.set_ylabel(r'$\kappa$ [g/cm$^2$]', fontsize=20)
    ax.set_xlim(lmin, lmax)
    
    for i in filelist:
        idnb = i.split(".")[0].split("-")[-1]
        kappa = pd.read_table(i, sep="\s+", comment='#', header=None, skiprows=9)  
        ax.loglog(kappa[0], kappa[1], label=' absorption {}'.format(idnb), linewidth=2, color = 'black', linestyle='-')  
        ax.loglog(kappa[0], kappa[2], label='scattering {}'.format(idnb), linewidth=2, color = 'black', linestyle='--')
        
    ax.legend(fontsize=16)
    ax.tick_params(labelsize=16)
    plt.show()



#_________________________________________________#
#                      MAIN                       #
#_________________________________________________#
if __name__ == "__main__":
    #define parser
    parser = argparse.ArgumentParser(description='Create the dustkappa* files for RADMC3D using the dsharp_opac package.')


    parser.add_argument('-ice', '--fmice', type=float,
						help='ice fraction should be 0.0 or 0.2.')

    parser.add_argument('-amin', '--amin', type=float,
						help='min size value in microns.')

    parser.add_argument('-amax', '--amax', type=float,
						help='max size value in microns.')

    parser.add_argument('-nlam', '--nlam', type=int,
						help='number of wavelengths. Typically between 200 and 250.')

    parser.add_argument('-filename', '--filename', type=str,
						help='choose a filename. Please do not include an underscore (_) in the name.')

    parser.add_argument('-scatmat', '--scatmat', type=bool,
						help='True or False. If True, then the matrix elements will be written in the dustkapscatmat_*.inp file. Default value is False.')
						
    parser.add_argument('-fig', '--fig', type=bool,
						help='True or False. If True, then a figure is gerenated after the opacities are computed.')
						

    args = parser.parse_args()
    #------------

    if args.fmice is not None:
        fm_ice = args.fmice
    else:
        fm_ice = 0.2

    if args.amin is not None:
        amin = args.amin
    else:
        amin = 5e-3 #microns

    if args.amax is not None:
        amax = args.amax
    else:
        amax = 1e3 #microns

    if args.nlam is not None:
        nlam = args.nlam
    else:
        nlam = 250

    if args.filename is not None:
        filename = args.filename
    else:
        filename = 'standard'

    lamin = 1e-1 #microns
    lamax = 2e3 #microns
    a = np.logspace(np.log10(amin), np.log10(1e3), 200) #microns
    lam = np.logspace(np.log10(lamin), np.log10(lamax), nlam) #microns


    if args.scatmat is not None:
        if args.scatmat == True:
            dustkapscatmat(a, amax, lam, fm_ice, filename)
        if args.scatmat == False: 
            dustkappa(a, amax, lam, fm_ice, filename)
    else:
        dustkappa(a, amax, lam, fm_ice, filename)

    if args.fig is not None:
        if args.fig == True:
            make_figure(lamin, lamax)
        if args.fig == False: 
            pass
    else:
        make_figure(lamin, lamax)