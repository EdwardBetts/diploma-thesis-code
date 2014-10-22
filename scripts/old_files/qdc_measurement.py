from subprocess import check_call
import os
import shutil
from numpy import loadtxt, savetxt

def measure_pde(time, wavelength, dst_dir, channels_to_safe):
    process_string = '/home/astro/jammer/Programme_vme/'
    program_string = 'v965_1'

    filename = str(wavelength)
    #check if file already exists
    if not os.path.isfile(''.join((dst_dir, filename))):
        filename = ''.join((filename, '_2'))

    #start measrement one
    subprocess.check_call([process_string + program_string,str(time)])

    #open histogram-file
    data = loadtxt(process_string + 'data_v965.dat', unpack=True)
    #write important channels to new file
    with open(''.join((dst_dir, filename)), 'w') as f:
        savetxt(f, [data[channels_to_safe]])

    #wait for raw_input to change fibers
    print 'First measurement complete, switch fibers and continue'
    raw_input()

    #start measurement two
    subprocess.check_call([process_string + program_string,str(time)])

    #open histogram-file
    data = loadtxt(process_string + 'data_v965.dat', unpack=True)
    #write important channels to new file
    with open(''.join((dst_dir, filename)), 'a') as f:
        savetxt(f, [data[channels_to_safe]])
