try:
	import argparse
except:
	print('Error: Unable to access the argparse module\n')
	sys.exit(1)
	
import os
import sys

try:
	import numpy as np 
except:
	print('Error: Unable to access the numpy module\n')
	sys.exit(1)

try:
	import mayavi.mlab as mlab
except:
	print('Error: Unable to access the mayavi module\nNote that Mayavi also has dependencies on the VTK library\n')
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
    INPUTS:  ro - a vector of radii (Default = 4.0)
             n - number of coordinates to return for each element in r (Default = 40)

    OUTPUTS: x,y,z - coordinates of points
    '''

	ro = 5.00
	n = 40

	points = np.linspace(0,2*np.pi,n+1)
	x = np.atleast_2d((ro-1) * np.cos(points))
	y = np.atleast_2d((ro-1) * np.sin(points))

	x = np.concatenate((x, x), axis = 0)
	y = np.concatenate((y, y), axis = 0)

	z = np.zeros((2,n+1), dtype = float)
	z[1,:] = z[1,:] + 0.875

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

	The function assumes there are 1.7 base-pair turns around each nucleosome 
	Thus, 170 equal intervals are built around a total angle of (4*pi)

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

	So, just hit mayavi's surface plot function on (Xc,Yc,Zc) to draw the cylinder
	'''

	X, Y, Z = cylinder()

	
	Xc = r[0] + (core_x[0]*X) + (core_y[0]*Y) + (core_z[0]*Z)
	Yc = r[1] + (core_x[1]*X) + (core_y[1]*Y) + (core_z[1]*Z)
	Zc = r[2] + (core_x[2]*X) + (core_y[2]*Y) + (core_z[2]*Z)

	#Rotating the cylinder according to the orientiation Matrices core_x, core_y, core_z
	#Translating the cylinder to the coordinate of 'r'

	return Xc,Yc,Zc
##############################################################################################################
##############################################################################################################
def draw_DNA(plottedDNA):

	numDNAbeads = plottedDNA.shape[0]
	cc = (1,0,0)

	mlab.plot3d(plottedDNA[:,0], plottedDNA[:,1], plottedDNA[:,2], color = cc, tube_radius = 0.65, tube_sides = 20)

##############################################################################################################
##############################################################################################################

def draw_LH(linker_histone_coods1, linker_histone_coods2, linker_histone_coods3, charges):

	sph_rad = 1.0

	phi = np.linspace(0, 2*np.pi, 50, endpoint = True, dtype = float)
	theta = np.linspace(0, np.pi, 50, endpoint = True, dtype = float)

	x = sph_rad * np.outer(np.cos(phi), np.sin(theta))
	y = sph_rad * np.outer(np.sin(phi), np.sin(theta))
	z = sph_rad * np.outer(np.ones(np.size(phi)), np.cos(theta))

		
	C = (0.690196, 0.878431, 0.901961) #cyan
	Y = (1, 1, 0) #yellow
	M = (1, 0, 1) #magenta

	cco = np.zeros((3), dtype = str)

	if charges.size == 0:
		cco[0] = 'C'
		cco[1] = 'C'
		cco[2] = 'C'

	else:
		cco[np.argmin(charges)] = 'Y'
		cco[np.argmax(charges)] = 'M'
		cco[3 - (np.argmax(charges) + np.argmin(charges))] = 'C'
		
	for i in xrange(0, linker_histone_coods1.shape[0]):

		
		plot_x = linker_histone_coods1[i,0] + x
		plot_y = linker_histone_coods1[i,1] + y
		plot_z = linker_histone_coods1[i,2] + z

		mlab.mesh(plot_x, plot_y, plot_z, color = eval(cco[0]))

		plot_x = linker_histone_coods2[i,0] + x
		plot_y = linker_histone_coods2[i,1] + y
		plot_z = linker_histone_coods2[i,2] + z

		mlab.mesh(plot_x, plot_y, plot_z, color = eval(cco[1]))

		plot_x = linker_histone_coods3[i,0] + x
		plot_y = linker_histone_coods3[i,1] + y
		plot_z = linker_histone_coods3[i,2] + z

		mlab.mesh(plot_x, plot_y, plot_z, color = eval(cco[2]))

##############################################################################################################
##############################################################################################################
def draw_histtail(histonetail_coods):

	sph_rad = 0.5

	phi = np.linspace(0, 2*np.pi, 25, endpoint = True, dtype = float)
	theta = np.linspace(0, np.pi, 25, endpoint = True, dtype = float)

	x = sph_rad * np.outer(np.cos(phi), np.sin(theta))
	y = sph_rad * np.outer(np.sin(phi), np.sin(theta))
	z = sph_rad * np.outer(np.ones(np.size(phi)), np.cos(theta))

	G = (0.596078, 0.984314, 0.596078)

	for i in xrange(0, histonetail_coods.shape[0]):
		plot_x = x + histonetail_coods[i,0]
		plot_y = y + histonetail_coods[i,1]
		plot_z = z + histonetail_coods[i,2]

		mlab.mesh(plot_x, plot_y, plot_z, color = G)
##############################################################################################################
##############################################################################################################
def calc_pk(core_coods_calc):

	core_array = np.linspace(1, ncalc_cores, ncalc_cores)
	tempmat = np.zeros((6), dtype = float)
	diffmat = np.zeros((5), dtype = float)
	fracmat = np.zeros((5), dtype = float)

	for k in xrange(1,7):

		aa1 = np.polyfit(core_array, core_coods_calc[:,0], k)
		aa2 = np.polyfit(core_array, core_coods_calc[:,1], k)
		aa3 = np.polyfit(core_array, core_coods_calc[:,2], k)

		fitx = np.polyval(aa1, core_array)
		fity = np.polyval(aa2, core_array)
		fitz = np.polyval(aa3, core_array)

		resx = fitx - core_coods_calc[:,0]
		resy = fity - core_coods_calc[:,1]
		resz = fitz - core_coods_calc[:,2]

		tempmat[k-1] = np.sum(np.sqrt(np.square(resx) + np.square(resy) + np.square(resz)))

		if k > 1:
			diffmat[k-2] = tempmat[k-2] - tempmat[k-1]
			fracmat[k-2] = (diffmat[k-2]/tempmat[k-2]) * 100.0

	if np.sum(np.int32(fracmat < 3.34)) == 0:
		pol_deg = 5
	else:
		pol_deg = np.argmax(fracmat < 3.34) + 1

	aa1 = np.polyfit(core_array, core_coods_calc[:,0], pol_deg)
	aa2 = np.polyfit(core_array, core_coods_calc[:,1], pol_deg)
	aa3 = np.polyfit(core_array, core_coods_calc[:,2], pol_deg)

	fitx = np.polyval(aa1, core_array)
	fity = np.polyval(aa2, core_array)
	fitz = np.polyval(aa3, core_array)

	step = 2
	FL_core_array = np.arange(1, n_cores+1, step) - 1
	n_FLcores = np.size(FL_core_array)

	fit_FL = np.zeros((n_FLcores,3), dtype = float)

	fit_FL[:,0] = fitx[FL_core_array]
	fit_FL[:,1] = fity[FL_core_array]
	fit_FL[:,2] = fitz[FL_core_array]

	FL = 0
	for k in xrange(0, n_FLcores-1):
		FL = FL + np.linalg.norm(fit_FL[k+1,:] - fit_FL[k,:])

	pack_ratio = (11*n_cores)/FL
	return pack_ratio

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
-pht 				: Plot core histone tails (Y/N). DEFAULT = N
-plh				: Plot linker histone (Y/N). DEFAULT = Y
-pk 				: Calculate and print the packing ratio. DEFAULT = N
-fr 				: File with frame numbers to be rendered
-lc 				: File with charges of Linker Histones 1,2,3. DEFAULT = Equal Charges

Other options
-------------------------------------------------------------------------
-h, --help			: Show this menu and exit

''')


parser.add_argument('-f', nargs = 1, dest = 'xyzfile', help = argparse.SUPPRESS, required = True)
parser.add_argument('-p', nargs = 1, dest = 'histfile', help = argparse.SUPPRESS, required = True)

parser.add_argument('-fc', nargs = 1, dest = 'first_core', default = [1], type = int, help = argparse.SUPPRESS)
parser.add_argument('-lh', nargs = 1, dest = 'togglelinker', default = ['Y'], help = argparse.SUPPRESS)
parser.add_argument('-plh', nargs = 1, dest = 'plotlinker', default = ['Y'], help = argparse.SUPPRESS)
parser.add_argument('-pht', nargs = 1, dest = 'plothisttail', default = ['N'], help = argparse.SUPPRESS)
parser.add_argument('-pk', nargs = 1, dest = 'calc_pk', default = ['N'], help = argparse.SUPPRESS)
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
args.plothisttail = args.plothisttail[0]
args.calc_pk = args.calc_pk[0]
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

if args.plothisttail != 'Y' and args.plothisttail != 'N':
	print('Only values of "Y" and "N" are allowed with option -pht\n')
	sys.exit(1)

if args.calc_pk != 'Y' and args.calc_pk != 'N':
	print('Only values of "Y" and "N" are allowed with option -pk\n')
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
			linker_charge = np.abs(linker_charge)

###################################################################################################################################################
###################################################################################################################################################

histonedetails = np.genfromtxt(args.histfile, dtype = int, usecols = 0)
n_cores = histonedetails[0]
n_linkers = np.sum(histonedetails) - n_cores

if args.first_core > n_cores:
	Errorout = 'The first core entered with the -fc is greater the number of cores in co-ordinate file - ' + str(n_cores) + '\n'
	print(Errorout)
	sys.exit(1)


if args.togglelinker == 'Y':
	lines_per_frame = (4*n_cores) + (4*n_linkers) + (53*n_cores)
elif args.togglelinker == 'N':
	lines_per_frame = (4*n_cores) + (4*n_linkers) + (50*n_cores)
#Linker Histones have 3 beads.
#Number of linker histones = Number of Cores
#So, number of linker histone beads = Number of Cores * 3


xyz = np.genfromtxt(args.xyzfile, dtype = float)
n_frames = int(np.floor(xyz.shape[0]/lines_per_frame))

if args.framefile == 'Auto':

	printframe = np.zeros((1), dtype = int)
	mlab.options.offscreen = False

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

	#The following three lines are the options for offscreen rendering without GUI.
	#Mayavi uses the vtk library to output. So, deactivating mayavi GUI is insufficient as the vtk GUI stops batch processsing.
	#So, we set the VTK export options here to stop the vtk GUI from asking confirmations. 
	mlab.options.offscreen = True
	try:
		from tvtk.api import tvtk
		exp_options = tvtk.GL2PSExporter(file_format = 'eps', sort = 'bsp', compress = 1)
	except:
		print('ERROR: Unable to access the vtk module.\n')
		sys.exit(1)
	
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
	
	deltaZ = n_cores * 4 #deltaZ....This is to count down the array of co-ordinates. As the Core xyz(s) are at the top.

	for dna_linker_index in xrange(0,n_linkers):
		dna_linker_coods[dna_linker_index,:] = currentcoods[deltaZ+(dna_linker_index*4),:]
		dna_linker_orientations[((dna_linker_index*3)+0),:] = currentcoods[deltaZ+(dna_linker_index*4)+1,:]
		dna_linker_orientations[((dna_linker_index*3)+1),:] = currentcoods[deltaZ+(dna_linker_index*4)+2,:] 
		dna_linker_orientations[((dna_linker_index*3)+2),:] = currentcoods[deltaZ+(dna_linker_index*4)+3,:]

	n_histonetails = n_cores * 50
	#Histones in the core each has a tail. Except the H2A has 2 tails each
	#Total tails = 8 + 2 = 10
	#Each tail is 5 beads. Total beads = 50 per core
	histonetail_coods = np.zeros((n_histonetails,3), dtype = float)

	deltaZ = deltaZ + (n_linkers*4)
	#Updating the deltaZ with the indexes occupied by the linker DNA co-ordinates

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
	#This links_sum is to help skip the indices of the linker DNA associated with non-plotted core. 

	rem = args.first_core % 2

	
	plottedDNAX = np.empty((0), dtype = float)
	plottedDNAY = np.empty((0), dtype = float)
	plottedDNAZ = np.empty((0), dtype = float)
	#The linker and nucleosome bound DNA are represented similarly.
	#Instead of passing them to the plotting function separately
	#Store all the DNA co-ordinates for this frame in these three arrays and pass them together. 
	

	for core_index in xrange(args.first_core-1,n_cores):

		rem2 = (core_index+1) % 2

		r = core_coods[core_index,:]
		core_x = core_orientations[(core_index*3)+0,:]
		core_y = core_orientations[(core_index*3)+1,:]
		core_z = core_orientations[(core_index*3)+2,:]

		Xc, Yc, Zc = nucleosome_cylinder(r, core_x, core_y, core_z)

		if rem == rem2:
			cc = (0,0,1)
			#Colour is in RGB format. So, (0,0,1) is Blue
		elif rem != rem2:
			cc = (0.33333, 0.33333, 0.33333)
			#(0.33333, 0.33333, 0.33333) tuple represents a dark grey colour. 
		
		#This alternative colouring of blue and grey of the nucleosome core allows easier distingushing between the solenoidal and zig-zag chromatin
		#In Solenoidal chromatin, nucleosomes interact with (i-1) and (i+1)
		#In Zig-Zag chromatin, nucleosomes interact with (i-2) and (i+2) 


		mlab.mesh(Xc, Yc, Zc, color = cc)
		
		Xn, Yn, Zn = wrapper_dna_coods(r, core_x, core_y, core_z)
		
		current_linker_coodsX = dna_linker_coods[links_sum:(links_sum+histonedetails[core_index+1]),0]
		current_linker_coodsY = dna_linker_coods[links_sum:(links_sum+histonedetails[core_index+1]),1]
		current_linker_coodsZ = dna_linker_coods[links_sum:(links_sum+histonedetails[core_index+1]),2]

		plottedDNAX = np.concatenate((plottedDNAX, Xn, current_linker_coodsX), axis = 0)
		plottedDNAY = np.concatenate((plottedDNAY, Yn, current_linker_coodsY), axis = 0)
		plottedDNAZ = np.concatenate((plottedDNAZ, Zn, current_linker_coodsZ), axis = 0)

		links_sum = links_sum + histonedetails[core_index+1]

	plottedDNA = np.concatenate((np.transpose(np.atleast_2d(plottedDNAX)),np.transpose(np.atleast_2d(plottedDNAY)),np.transpose(np.atleast_2d(plottedDNAZ))), axis = 1)
	#This plottedDNA variable has wrapper and linker DNA co-ordinates combined.
	
	draw_DNA(plottedDNA)
	
	if args.plotlinker == 'Y':

		plotted_linker_coods1 = linker_histone_coods1[args.first_core-1:,:]
		plotted_linker_coods2 = linker_histone_coods2[args.first_core-1:,:]
		plotted_linker_coods3 = linker_histone_coods3[args.first_core-1:,:]
		##This is to avoid plotting Linker histones attached to the nucleosome cores not plotted. 


		if args.linkercharges == 'Auto':
			draw_LH(plotted_linker_coods1, plotted_linker_coods2, plotted_linker_coods3, np.empty((0)))
		else:
			draw_LH(plotted_linker_coods1, plotted_linker_coods2, plotted_linker_coods3, linker_charge)
	
	if args.plothisttail == 'Y':
		first_plot_hist_tail = (args.first_core-1)*50
		plotted_hist_tail = histonetail_coods[first_plot_hist_tail:,:]
		draw_histtail(plotted_hist_tail)

	if args.calc_pk == 'Y':
		core_coods_calc = core_coods[args.first_core-1:,:]
		ncalc_cores = core_coods_calc.shape[0]

		if ncalc_cores == 1:
			print('Cannot calculate packing ratios for 1 core. Skipping Calculation\n')

		else:
			pack_ratio = calc_pk(core_coods_calc)
			print_text = 'Packing Ratio\n' + str(round(pack_ratio,2))

			mlab.text3d(core_coods_calc[0,0]+10, core_coods_calc[0,1]+10, core_coods_calc[0,2]+10, print_text, scale = 3, color = (0,0,0))


	CF = mlab.gcf()
	###############################################################################################################################
	###############################################################################################################################
	#These are the lighting options. Open to change... 
	camera_light0 = CF.scene.light_manager.lights[0]
	camera_light0.elevation = 40
	camera_light0.azimuth = 60
	camera_light0.intensity = 1.0
	
	camera_light1 = CF.scene.light_manager.lights[1]
	camera_light1.elevation = -40
	camera_light1.azimuth = -60
	camera_light1.intensity = 1.0
	
	camera_light2 = CF.scene.light_manager.lights[2]
	camera_light2.intensity = 0.0
	#Turning off the third light source. 
	#It makes the image too bright for the visualizer to render 'pleasing to the eye'.
	###############################################################################################################################
	###############################################################################################################################
	mlab.view(azimuth = 60, elevation = 45, distance = 200)
	
	#mlab.yaw(25)
	#mlab.pitch(90)
	#Uncomment and modify these two options to rotate the 'View' vertically and horizontally


	mlab.figure(figure = CF, bgcolor=(1.0,1.0,1.0), fgcolor = (1.0,1.0,1.0))
	#Manually setting the BackGround and ForeGround colour to white
	#MayaVi automatically sets it to Grey. 

	fname = 'frame' + str(currentframe) + '.eps'


	if args.framefile == 'Auto':
		mlab.show()

	else:
		#Turn this option off if you want interactive plotting
		#Currently it is turned on for batch processing
		CF.scene.save_gl2ps(fname, exp_options)
