import argparse
import os
import sys

try:
	import numpy as np 
except:
	print('Error: Unable to access the numpy module\n')
	sys.exit(1)

try:
	import matplotlib.pyplot as plt
	from mpl_toolkits.mplot3d import Axes3D
        from matplotlib.colors import LightSource
        from matplotlib import cm
except:
	print('Error: Unable to access the matplotlib module\n')
	sys.exit(1)


##############################################################################################################
def unitvector(v):
	'''
	Returns the unit-vector of any 3D vector
	'''
	UV = v/np.linalg.norm(v)
	return UV
##############################################################################################################
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

	x = np.concatenate((x, x), axis = 0)
	y = np.concatenate((y, y), axis = 0)

	z = np.zeros((2,n+1), dtype = float)
	z[1,:] = z[1,:] + 0.75

	A = np.zeros((1,41), dtype = float)
	B = np.ones((1,41), dtype = float)


	X = np.concatenate((A, x, A), axis = 0)
	Y = np.concatenate((A, y, A), axis = 0)
	Z = np.concatenate((A, z, B), axis = 0)

	Z = (ro*Z) - (ro/2)

	return X,Y,Z
##############################################################################################################	
##############################################################################################################

def wrapper_dna_coods(r, core_x, core_y, core_z):

	'''
	If provided the co-ordinates of the nucleosome core 'r' and its orientations:
	This function returns co-ordinates of the DNA wrapped around the core

	The function assumes there are 171 base-pairs around each nucleosome with 2 wraps
	Thus, a base pair is placed at 171 equal intervals around a total angle of (4*pi)

	The deltaZ is set at 1.8 units (The height of the cylinder signifying the nucleosome core)
	'''

	ro = 4.8
	d1 = 1.8

	theta_0 = 0.6 * np.pi
	thetas = theta_0 + np.linspace(0, ((4*np.pi) - theta_0), num = 171, retstep = False, endpoint = True)

	Z = d1 - ((thetas - theta_0) * ((2*d1)/((4*np.pi) - theta_0)))
	X = ro * np.cos(thetas - (np.pi/2))
	Y = ro * np.sin(thetas - (np.pi/2))

	Xn = r[0] + (core_x[0]*X) + (core_y[0]*Y) + (core_z[0]*Z)
	Yn = r[1] + (core_x[1]*X) + (core_y[1]*Y) + (core_z[1]*Z)
	Zn = r[2] + (core_x[2]*X) + (core_y[2]*Y) + (core_z[2]*Z)

	return Xn, Yn, Zn
##############################################################################################################
##############################################################################################################
	
def nucleosome_cylinder(r,core_x,core_y,core_z):

	'''
	If provided the co-ordinates of the nucleosome core 'r' and its orientations:
	This function returns the representative co-ordinates of the cylinder describing the nucelosome core

	So, just hit python's representative of surf(Xc,Yc,Zc) to draw the cylinder
	'''

	X, Y, Z = cylinder()

	
	Xc = r[0] + (core_x[0]*X) + (core_y[0]*Y) + (core_z[0]*Z)
	Yc = r[1] + (core_x[1]*X) + (core_y[1]*Y) + (core_z[1]*Z)
	Zc = r[2] + (core_x[2]*X) + (core_y[2]*Y) + (core_z[2]*Z)

	return Xc,Yc,Zc
##############################################################################################################
##############################################################################################################

def draw_DNA(plottedDNA, ax, histonedetails):

	'''
	This is python's version of surface used to draw the DNA
	Provide the 'ax' figure class to it and the co-ordinates of the DNA
	'''

	numDNAdrawn = plottedDNA.shape[0]

	DNAthetas = np.linspace(0, 2*np.pi, 100)
	xc = np.atleast_2d(np.cos(DNAthetas))
	yc = np.atleast_2d(np.sin(DNAthetas))
	zc = np.atleast_2d(np.zeros((xc.size), dtype = float))

	########################################################################################
	#Now plotting the zeroth element.
	#Don't want an if statement slowing down the loop iteration.
	########################################################################################
	BBtheta = np.concatenate((xc,yc,zc), axis = 0)
	FFtheta = BBtheta

	deltaDNA = unitvector(plottedDNA[2,:] - plottedDNA[0,:])
	sph_theta = np.arctan2(deltaDNA[0], deltaDNA[1])
	sph_phi = np.arccos(deltaDNA[2])

	Rotmat = np.zeros((3,3), dtype = float)
	Rotmat[0,0] = np.cos(sph_theta)
	Rotmat[0,1] = np.cos(sph_phi) * np.sin(sph_theta)
	Rotmat[0,2] = np.sin(sph_phi) * np.sin(sph_theta)
	Rotmat[1,0] = -1 * np.sin(sph_theta)
	Rotmat[1,1] = np.cos(sph_phi) * np.cos(sph_theta)
	Rotmat[1,2] = np.sin(sph_phi) * np.cos(sph_theta)
	Rotmat[2,0] = 0.0
	Rotmat[2,1] = -1 * np.sin(sph_phi)
	Rotmat[2,2] = np.cos(sph_phi)


	BB = np.dot(Rotmat, BBtheta)
	FF = np.dot(Rotmat, FFtheta)

	BB = BB - np.transpose(np.atleast_2d(np.mean(BB, axis = 1))) + np.transpose(np.atleast_2d(plottedDNA[0,:]))
	BB = BB - np.transpose(np.atleast_2d(np.mean(BB, axis = 1))) + np.transpose(np.atleast_2d(plottedDNA[0,:]))

	FF = FF - np.transpose(np.atleast_2d(np.mean(FF, axis = 1))) + np.transpose(np.atleast_2d(plottedDNA[1,:]))

	plot_x = np.concatenate((np.atleast_2d(BB[0,:], np.atleast_2d(FF[0,:]))), axis = 0)
	plot_y = np.concatenate((np.atleast_2d(BB[1,:], np.atleast_2d(FF[1,:]))), axis = 0)
	plot_z = np.concatenate((np.atleast_2d(BB[2,:], np.atleast_2d(FF[2,:]))), axis = 0)

	ax.plot_surface(plot_x, plot_y, plot_z, color = 'r', rstride=10, cstride=10, antialiased=True, linewidth = 0, zorder = 5)
	##############################################################################################################
	#Now plotting the rest of the elements until the second last one
	##############################################################################################################
		
	for j in xrange(1,numDNAdrawn-2):
		
		BB = FF
		
		deltaDNA = unitvector(plottedDNA[j+2,:] - plottedDNA[j,:])
		
		sph_theta = np.arctan2(deltaDNA[0], deltaDNA[1])
		sph_phi = np.arccos(deltaDNA[2])

		Rotmat[0,0] = np.cos(sph_theta)
		Rotmat[0,1] = np.cos(sph_phi) * np.sin(sph_theta)
		Rotmat[0,2] = np.sin(sph_phi) * np.sin(sph_theta)
		Rotmat[1,0] = -1 * np.sin(sph_theta)
		Rotmat[1,1] = np.cos(sph_phi) * np.cos(sph_theta)
		Rotmat[1,2] = np.sin(sph_phi) * np.cos(sph_theta)
		#Rotmat[2,0] = 0.0   (Not updating this every iteration.)
		Rotmat[2,1] = -1 * np.sin(sph_phi)
		Rotmat[2,2] = np.cos(sph_phi)

		FFtheta[0,:] = (np.cos(sph_theta) * xc) + (np.sin(sph_theta) * yc)
		FFtheta[1,:] = (-1 * np.sin(sph_theta) * xc) + (np.cos(sph_theta)*yc)
		FF = np.dot(Rotmat, FFtheta)


		FF = FF - np.transpose(np.atleast_2d(np.mean(FF, axis = 1))) + np.transpose(np.atleast_2d(plottedDNA[j+1,:]))

		plot_x = np.concatenate((np.atleast_2d(BB[0,:], np.atleast_2d(FF[0,:]))), axis = 0)
		plot_y = np.concatenate((np.atleast_2d(BB[1,:], np.atleast_2d(FF[1,:]))), axis = 0)
		plot_z = np.concatenate((np.atleast_2d(BB[2,:], np.atleast_2d(FF[2,:]))), axis = 0)

		ax.plot_surface(plot_x, plot_y, plot_z, color = 'r', zorder = 5, rstride=10, cstride=10, antialiased=True, linewidth = 0)

	##############################################################################################################
	#Now plotting the second last DNA base-par. I dont want an if-else statement within the Loop
	##############################################################################################################

	BB = FF
			
	deltaDNA = unitvector(plottedDNA[numDNAdrawn-1,:] - plottedDNA[numDNAdrawn-2,:])
	
	sph_theta = np.arctan2(deltaDNA[0], deltaDNA[1])
	sph_phi = np.arccos(deltaDNA[2])

	Rotmat[0,0] = np.cos(sph_theta)
	Rotmat[0,1] = np.cos(sph_phi) * np.sin(sph_theta)
	Rotmat[0,2] = np.sin(sph_phi) * np.sin(sph_theta)
	Rotmat[1,0] = -1 * np.sin(sph_theta)
	Rotmat[1,1] = np.cos(sph_phi) * np.cos(sph_theta)
	Rotmat[1,2] = np.sin(sph_phi) * np.cos(sph_theta)
	#Rotmat[2,0] = 0.0   (Not updating this every iteration.)
	Rotmat[2,1] = -1 * np.sin(sph_phi)
	Rotmat[2,2] = np.cos(sph_phi)

	FFtheta[0,:] = (np.cos(sph_theta) * xc) + (np.sin(sph_theta) * yc)
	FFtheta[1,:] = (-1 * np.sin(sph_theta) * xc) + (np.cos(sph_theta)*yc)
	FF = np.dot(Rotmat, FFtheta)


	FF = FF - np.transpose(np.atleast_2d(np.mean(FF, axis = 1))) + np.transpose(np.atleast_2d(plottedDNA[j+1,:]))

	plot_x = np.concatenate((np.atleast_2d(BB[0,:], np.atleast_2d(FF[0,:]))), axis = 0)
	plot_y = np.concatenate((np.atleast_2d(BB[1,:], np.atleast_2d(FF[1,:]))), axis = 0)
	plot_z = np.concatenate((np.atleast_2d(BB[2,:], np.atleast_2d(FF[2,:]))), axis = 0)

	ax.plot_surface(plot_x, plot_y, plot_z, color = 'r', zorder = 5, rstride=10, cstride=10, antialiased=True, linewidth = 0)
###################################################################################################################################################
###################################################################################################################################################

def draw_ball(linker_histone_coods1, linker_histone_coods2, linker_histone_coods3, ax, charges):

	sph_rad = 1.0

	phi = np.linspace(0, 2*np.pi, 25, endpoint = True, dtype = float)
	theta = np.linspace(0, np.pi, 25, endpoint = True, dtype = float)

	x = sph_rad * np.outer(np.cos(phi), np.sin(theta))
	y = sph_rad * np.outer(np.sin(phi), np.sin(theta))
	z = sph_rad * np.outer(np.ones(np.size(phi)), np.cos(theta))


	colormat = np.empty((3), dtype = str)
	if charges.size == 0:
		colormat[0] = 'c'
		colormat[1] = 'c'
		colormat[2] = 'c'

	else:
		colormat[np.argmax(charges)] = 'm'
		colormat[np.argmin(charges)] = 'c'
		colormat[3 - np.argmin(charges) - np.argmax(charges)] = 'y'


	for i in xrange(0, linker_histone_coods1.shape[0]):

		plot_x = linker_histone_coods1[i,0] + x
		plot_y = linker_histone_coods1[i,1] + y
		plot_z = linker_histone_coods1[i,2] + z

		ax.plot_surface(plot_x, plot_y, plot_z, color = colormat[0], rstride=1, cstride=1, antialiased=True, linewidth = 0, zorder = 2)

		plot_x = linker_histone_coods2[i,0] + x
		plot_y = linker_histone_coods2[i,1] + y
		plot_z = linker_histone_coods2[i,2] + z

		ax.plot_surface(plot_x, plot_y, plot_z, color = colormat[1], rstride=1, cstride=1, antialiased=True, linewidth = 0, zorder = 2)

		plot_x = linker_histone_coods3[i,0] + x
		plot_y = linker_histone_coods3[i,1] + y
		plot_z = linker_histone_coods3[i,2] + z

		ax.plot_surface(plot_x, plot_y, plot_z, color = colormat[2], rstride=1, cstride=1, antialiased=True, linewidth = 0, zorder = 2)

###################################################################################################################################################
###################################################################################################################################################
		
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
-fc 				: The index of the first core to be plotted. DEFAULT = 1
-lh				: Linker histone present in co-ordinate file (Y/N). DEFAULT = Y
-plh				: Plot linker histone (Y/N). DEFAULT = Y
-fr 				: File with frame numbers to be rendered
-lc 				: File with charges of Linker Histones 1,2,3. Default = Equal Charges

Other options
-------------------------------------------------------------------------
-h, --help			: Show this menu and exit

''')


parser.add_argument('-f', nargs = 1, dest = 'xyzfile', help = argparse.SUPPRESS, required = True)
parser.add_argument('-p', nargs = 1, dest = 'histfile', help = argparse.SUPPRESS, required = True)

parser.add_argument('-fc', nargs = 1, dest = 'first_core', default = [1], type = int, help = argparse.SUPPRESS)
parser.add_argument('-lh', nargs = 1, dest = 'togglelinker', default = ['Y'], help = argparse.SUPPRESS)
parser.add_argument('-plh', nargs = 1, dest = 'plotlinker', default = ['Y'], help = argparse.SUPPRESS)
parser.add_argument('-fr', nargs = 1, dest = 'framefile', default = ['Auto'], help = argparse.SUPPRESS)
parser.add_argument('-lc', nargs = 1, dest = 'linkercharges', default = ['Auto'], help = argparse.SUPPRESS)

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
args.linkercharges = args.linkercharges[0]

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

if args.plotlinker != 'Y' and args.plotlinker != 'N':
	print('Only values of "Y" and "N" are allowed with option -plh\n')
	sys.exit(1)

elif args.togglelinker == 'N' and args.plotlinker == 'Y':
	print('Cannot plot linker Histones not in co-ordinate file. Changing -plh option to "N" from "Y"\n')
	args.plotlinker = 'N'

if args.plotlinker == 'Y' and args.linkercharges != 'Auto':
	if not os.path.isfile(args.linkercharges):
		print('File containing charges of linker histones does not exist. Plotting with equal charges\n')
		args.linkercharges = 'Auto'

	else:
		linker_charge = np.genfromtxt(args.linkercharges, dtype = float, usecols = 0)
		if np.sum(np.int32(np.isnan(linker_charge))) > 0:
			output = 'Check the charges in the file ' + args.linkercharges + '. Not a number entered.\nPlotting with equal charges\n'
			print(output)
			args.linkercharges = 'Auto'

		elif linker_charge.size != 3:
			print('Currently hard-coded for 3 beads per Linker Histone. Plotting with equal charges\n')
			args.linkercharge = 'Auto'

		else:
			linkercharge = np.abs(linker_charge)

###################################################################################################################################################
###################################################################################################################################################

histonedetails = np.genfromtxt(args.histfile, dtype = int, usecols = 0)
n_cores = histonedetails[0]
n_linkers = np.sum(histonedetails) - n_cores

if args.first_core > n_cores:
	Errorour = 'The first core entered with the -fc is greater the number of cores in co-ordinate file - ' + str(n_cores) + '\n'
	print(Errorour)
	sys.exit(1)


if args.togglelinker == 'Y':
	lines_per_frame = (4*n_cores) + (4*n_linkers) + (53*n_cores)
elif args.togglelinker == 'N':
	lines_per_frame = (4*n_cores) + (4*n_linkers) + (50*n_cores)


xyz = np.genfromtxt(args.xyzfile, dtype = float)
n_frames = int(np.floor(xyz.shape[0]/lines_per_frame))

if args.framefile == 'Auto':

	printframe = np.zeros((1), dtype = int)

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
	printframe = np.genfromtxt(args.framefile, usecols = 0, dtype = int)

	if np.sum(np.int32(np.isnan(printframe))) > 0:
		ErrorOut = 'Check the frames in the input file ' + str(args.framefile) + '. Only Integers allowed.\n'
		print(ErrorOut)
		sys.exit(1)

	elif np.amin(printframe) <= 0:
		ErrorOut = 'Check the frames in the input file ' + str(args.framefile) + '. You have entered a negative frame number or a float.\n'
		print(ErrorOut)
		sys.exit(1)

	elif np.amax(printframe) > n_frames:
		ErrorOut = 'Check the frames in the input file ' + str(args.framefile) + '. Entered a frame greater than the total number of frames' + str(n_frames) + '\n'
		print(ErrorOut)
		sys.exit(1)



for i in xrange(0,printframe.size):


	###############################################################################################################################
	#Initializing the plot

	exec('fig' + str(i) + '= plt.figure()')
	exec("ax = fig" + str(i) + ".gca(projection = '3d')")
	ax.set_axis_off()
	###############################################################################################################################

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
	#Because we only want to plot from args.first_core onwards.. 
	#This links_sum is to help skip the linker DNA associated with non-plotted core. 

	rem = args.first_core % 2

	
	plottedDNAX = np.empty((0), dtype = float)
	plottedDNAY = np.empty((0), dtype = float)
	plottedDNAZ = np.empty((0), dtype = float)
	#The linker and nucleosome bound DNA are represented similarly.
	#So, store all the DNA coods for this frame in these three arrays and draw them together. 

	for core_index in xrange(args.first_core-1,n_cores):

		rem2 = (core_index+1) % 2

		r = core_coods[core_index,:]
		core_x = core_orientations[(core_index*3)+0,:]
		core_y = core_orientations[(core_index*3)+1,:]
		core_z = core_orientations[(core_index*3)+2,:]

		Xc, Yc, Zc = nucleosome_cylinder(r, core_x, core_y, core_z)

		if rem == rem2:
			cc = 'b'
		elif rem != rem2:
			cc = '0.75'

		ax.plot_surface(Xc, Yc, Zc, rstride=1, cstride=1, antialiased=True, color = cc, linewidth = 0, zorder = 0)

		Xn, Yn, Zn = wrapper_dna_coods(r, core_x, core_y, core_z)
		

		current_linker_coodsX = dna_linker_coods[links_sum:(links_sum+histonedetails[core_index+1]),0]
		current_linker_coodsY = dna_linker_coods[links_sum:(links_sum+histonedetails[core_index+1]),1]
		current_linker_coodsZ = dna_linker_coods[links_sum:(links_sum+histonedetails[core_index+1]),2]

		plottedDNAX = np.concatenate((plottedDNAX, Xn, current_linker_coodsX), axis = 0)
		plottedDNAY = np.concatenate((plottedDNAY, Yn, current_linker_coodsY), axis = 0)
		plottedDNAZ = np.concatenate((plottedDNAZ, Zn, current_linker_coodsZ), axis = 0)

		links_sum = links_sum + histonedetails[core_index+1]

	plottedDNA = np.concatenate((np.transpose(np.atleast_2d(plottedDNAX)),np.transpose(np.atleast_2d(plottedDNAY)),np.transpose(np.atleast_2d(plottedDNAZ))), axis = 1)
	
	draw_DNA(plottedDNA, ax, histonedetails)

	if args.plotlinker == 'Y':

		plotted_linker_coods1 = linker_histone_coods1[args.first_core-1:,:]
		plotted_linker_coods2 = linker_histone_coods2[args.first_core-1:,:]
		plotted_linker_coods3 = linker_histone_coods3[args.first_core-1:,:]
		##This is to avoid plotting Linker histones attached to nucleosome cores not plotted. 


		if args.linkercharges == 'Auto':
			draw_ball(plotted_linker_coods1, plotted_linker_coods2, plotted_linker_coods3, ax, np.empty((0)))
		else:
			draw_ball(plotted_linker_coods1, plotted_linker_coods2, plotted_linker_coods3, ax, linker_charge)
	
ax.set_aspect('equal')

ls1 = LightSource(225,-45)
ls2 = LightSource(225,45)
ls3 = LightSource(315,45)

ax.view_init(-20,6)

plt.tight_layout()
plt.show()
