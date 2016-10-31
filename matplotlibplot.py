try:
	import numpy as np 
except:
	print('Error: Unable to access the numpy module\n')
	sys.exit(1)

try:
	import matplotlib.pyplot as plt
except:
	print('Error: Unable to access the matplotlib module\n')
	sys.exit(1)

import argparse
import os
import sys


def cylinder():
	'''
	
    Returns the unit cylinder that corresponds to the curve r.
    INPUTS:  ro - a vector of radii (Default = 3.8)
             n - number of coordinates to return for each element in r (Default = 40)

    OUTPUTS: x,y,z - coordinates of points
    '''

    ro = 4.8
    n = 40

    points = np.linspace(0,2*np.pi,n+1)
    x = np.atleast_2d((ro-1) * np.cos(points))
    y = np.atleast_2d((ro-1) * np.sin(points))

    x = np.concactenate((x, x), axis = 1)
    y = np.concactenate((y, y), axis = 1)

    z = np.zeros((2,n+1), dtype = float)
    z[1,:] = z[1,:] + 1.0

    X = np.concactenate(((np.zeros(1,41, dtype = float)), x, (np.zeros(1,41, dtype = float))), axis = 1)
    Y = np.concactenate(((np.zeros(1,41, dtype = float)), y, (np.zeros(1,41, dtype = float))), axis = 1)
    Z = np.concactenate(((np.zeros(1,41, dtype = float)), z, (np.ones(1,41, dtype = float))), axis = 1)

    Z = (ro*Z) - (ro/2)

    return X,Y,Z

def wrapper_dna_coods(r, core_x, core_y, core_z):

	ro = 4.8
	d1 = 1.8

	theta_0 = 0.6 * np.pi
	thetas = np.linspace(theta_0, np.pi, num = 171, retstep = False, endpoint = True)

	Z = d1 - ((thetas - theta_0) * ((2*d1)/((4*np.pi) - theta_0)))
	X = ro * np.cos(thetas - (np.pi/2))
	Y = ro * np.sin(thetas - (np.pi/2))

	Xn = r[0] + (core_x[1]*X) + (core_y[1]*Y) + (core_z[1]*Z)
	Yn = r[1] + (core_x[2]*X) + (core_y[2]*Y) + (core_z[2]*Z)
	Zn = r[2] + (core_x[3]*X) + (core_y[3]*Y) + (core_z[3]*Z)

	return Xn,Yn,Zn


def nucleosome_cylinder(r,core_x,core_y,core_z):
	X, Y, Z = cylinder()

	Xc = r[0] + (core_x[1]*X) + (core_y[1]*Y) + (core_z[1]*Z)
	Yc = r[1] + (core_x[2]*X) + (core_y[2]*Y) + (core_z[2]*Z)
	Zc = r[2] + (core_x[3]*X) + (core_y[3]*Y) + (core_z[3]*Z)

	return Xc,Yc,Zc




parser = argparse.ArgumentParser(prog = 'CG Chromatin Visualizer', add_help = False,  formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
*************************************************************************
Chromatin CG Visualizer
author: Akshay Sridhar (as2637@cam.ac.uk)
git: https://github.com/akshay-sridhar/CGchromatin
*************************************************************************

Option 				Description                    
-------------------------------------------------------------------------
-f				: Co-ordinate file of bead positions
-p 				: File with details of histone cores
-o				: Output file (Format taken from this). DEFAULT = frame_xx.eps

Optional options
-------------------------------------------------------------------------
-fc 			: The index of the first core to be plotted. DEFAULT = 1
-lh				: Linker histone present (Y/N). DEFAULT = Y
-plh				: Plot linker histone (Y/N). DEFAULT = Y
-fr 				: File with frame numbers to be rendered 

Other options
-------------------------------------------------------------------------
-h, --help			: Show this menu and exit\n\n
''')


parser.add_argument('-f', nargs = 1, dest = 'xyzfile', help = argparse.SUPPRESS, required = True)
parser.add_argument('-p', nargs = 1, dest = 'histfile', help = argparse.SUPPRESS, required = True)

parser.add_argument('-fc', nargs = 1, dest = 'first_core', default = 1, type = int, help = argparse.SUPPRESS)
parser.add_argument('-lh', nargs = 1, dest = 'togglelinker', default = ['Y'], help = argparse.SUPPRESS)
parser.add_argument('-plh', nargs = 1, dest = 'plotlinker', default = ['Y'], help = argparse.SUPPRESS)
parser.add_argument('-fr', nargs = 1, dest = 'framefile', default = ['Auto'], help = argparse.SUPPRESS)

parser.add_argument('-o', nargs = 1, dest = 'outputfile', default = ['frame.eps'], help = argparse.SUPPRESS)
parser.add_argument('-h', '--help' , action = 'help', help = argparse.SUPPRESS)

args = parser.parse_args()

args.xyzfile = args.xyzfile[0]
args.histfile = args.histfile[0]
args.framefile = args.framefile[0]
args.outputfile = args.outputfile[0]
args.togglelinker = args.togglelinker[0]
args.plotlinker = args.plotlinker[0]
args.first_core = args.first_core[0]

if not os.path.isfile(args.xyzfile):
	print('Input XYZ coordinate file does not exist\n')
	sys.exit(1)

if not os.path.isfile(args.histfile):
	print('File containing details of histone cores does not exist\n')
	sys.exit(1)

if args.framefile != 'Auto':
	if not os.path.isfile(args.framefile):
		print('File containing the list of frames does not exist\n')
		sys.exit(1)

if args.togglelinker != 'Y' and args.togglelinker != 'N':
	print('Only values of "Y" and "N" are allowed with option -lh\n')
	sys.exit(1)

elif args.togglelinker == 'N':
	args.plotlinker = 'N'

if args.plotlinker != 'Y' and args.plotlinker != 'N':
	print('Only values of "Y" and "N" are allowed with option -plh\n')
	sys.exit(1)

histonedetails = np.genfromtxt(args.histfile, dtype = int, usecols = 0)
n_cores = histonedetails[0]
n_linkers = np.sum(histonedetails) - n_cores


if args.togglelinker == 'Y':
	lines_per_frame = (4*n_cores) + (4*n_linkers) + (53*n_cores)
elif args.togglelinker == 'N':
	lines_per_frame = (4*n_cores) + (4*n_linkers) + (50*n_cores)


xyz = np.genfromtxt(args.xyzfile, dtype = float, delimiter = '\t')
n_frames = int(np.floor(xyz.shape[0]/lines_per_frame))

if args.framefile == 'Auto':

	printframe = np.zeros((0), dtype = int)

	inputcheck = 'nok'
	while(inputcheck == 'nok'):

		tempframe = raw_input('Enter the frame number to be rendered:\n')
		try:
			tempframe = int(tempframe)
		except ValueError:
			print('Please Enter only a Number\n')
		else:
			if tempframe <= 0:
				print('Please Enter a positive value for frame number\n')

			elif tempframe > n_frames:
				ErrorOut = 'Frame number should be lesser than number of Frames in trajectory (<' + str(int(n_frames)) + ')\n'
				print(ErrorOut)

			else:
				printframe[0] = int(tempframe)
				inputcheck = 'ok'

else:
	printframe = np.genfromtxt(args.framefile, usecols = 0)

	if np.sum(np.int32(np.isnan(printframe))) > 0:
		ErrorOut = 'Check the frames in the input file ' + str(args.framefile) + '. Only Integers allowed.\n'
		print(ErrorOut)
		sys.exit(1)

	elif np.amin(printframe) <= 0:
		ErrorOut = 'Check the frames in the input file ' + str(args.framefile) + '. You have entered a negative frame number.\n'
		print(ErrorOut)
		sys.exit(1)

	elif np.amax(printframe) > n_frames:
		ErrorOut = 'Check the frames in the input file ' + str(args.framefile) + '. Entered a frame greater than the total number of frames' + str(n_frames) + '\n'
		print(ErrorOut)
		sys.exit(1)

for i in xrange(0,printframe.size):

	currentframe = printframe[i]
	starting_line = (currentframe - 1) * lines_per_frame
	ending_line = currentframe * lines_per_frame
	currentcoods = xyz[starting_line:ending_line, :]

	core_coods = np.zeros((n_cores,3), dtype = float)
	core_orientations = np.zeros((n_cores*3, 3), dtype = float)

	for core_index in xrange(0, n_cores):
		core_coods[core_index, :] = currentcoods[core_index*4, :]
		core_orientations[((core_index*3)+0),:] = currentcoods[((core_index*4)+1),:]
		core_orientations[((core_index*3)+1),:] = currentcoods[((core_index*4)+2),:]
		core_orientations[((core_index*3)+2),:] = currentcoods[((core_index*4)+3),:]


	dna_linker_coods = np.zeros((n_linkers,3), dtype = float)
	dna_linker_orientations = np.zeros((n_linkers*3,3), dtype = float)
	
	deltaZ = n_cores * 4 #deltaZ....This is to go down the array of co-ordinates. As the Core xyz(s) are at the top.

	for dna_linker_index in xrange(0,n_linkers):
		dna_linker_coods[dna_linker_index,:] = currentcoods[deltaZ+(dna_linker_index*4),:]
		dna_linker_orientations[((dna_linker_index*3)+0),:] = currentcoods[deltaZ+(dna_linker_index*4)+1,:]
		dna_linker_orientations[((dna_linker_index*3)+1),:] = currentcoods[deltaZ+(dna_linker_index*4)+2,:] 
		dna_linker_orientations[((dna_linker_index*3)+2),:] = currentcoods[deltaZ+(dna_linker_index*4)+3,:]

	n_histonetails = n_cores * 50
	histonetail_coods = np.zeros((n_histonetails,3), dtype = float)

	deltaZ = deltaZ + (n_linkers*4)

	for histonetail_index in xrange(0,n_histonetails):
		histonetail_coods[histonetail_index,:] = currentcoods[deltaZ+histonetail_index,:]


	if args.togglelinker == 'Y':

		deltaZ = deltaZ + n_histonetails
		linker_histone_coods1 = np.zeros((n_cores,3), dtype = float)
		linker_histone_coods2 = np.zeros((n_cores,3), dtype = float)
		linker_histone_coods3 = np.zeros((n_cores,3), dtype = float)

		for linker_histone_index in xrange(0,n_cores):

			linker_histone_coods1[linker_histone_index,:] = currentcoods[deltaZ+(linker_histone_index*3)+0,:]
			linker_histone_coods2[linker_histone_index,:] = currentcoods[deltaZ+(linker_histone_index*3)+1,:]
			linker_histone_coods3[linker_histone_index,:] = currentcoods[deltaZ+(linker_histone_index*3)+2,:]


	center = np.mean(core_coods, axis = 0)

	links_sum = np.sum(histonedetails[1:args.first_core])
	rem = args.first_core % 2

	colorindex = np.asarray(np.matrix('0,0,1;0,1,0;1,1,0;1,0,0;1,1,1;0.8,0.2,0;0.251,0.8784,0.8157'))

	for core_index in xrange(args.first_core-1,n_cores):

		rem2 = (core_index+1) % 2

		r = core_coods[core_index,:]
		core_x = core_orientations[(core_index*3)+1,:]
		core_y = core_orientations[(core_index*3)+2,:]
		core_z = core_orientations[(core_index*3)+3,:]

		Xc, Yc, Zc = nucleosome_cylinder(r, core_x, core_y, core_z)

		if rem == rem2:
			cc = 1
		elif rem != rem2:
			cc = 5

		Xn, Yn, Zn = wrapper_dna_coods(r, core_x, core_y, core_z)

		


		links_sum = links_sum + histonedetails[core_index+1]





