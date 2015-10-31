"""
In order to run this script, you should do the following in advance:

	1- Make sure you have VisIt package on your machine; for LRZ Linux cluster,
	one should load the VisIt module:

			>>> module load visit

	2- Make sure python is working nicely on your system; for LRZ Linux cluster
	, the python module should be loaded:

			>>> module load python

	3- Set the environmental variables: 

		i) Set PYTHONPATH environmental variables to the directory of the VisIt 
		library, since this library is imported to be used in this script. For
		example in order to do that on Linux Cluster in LRZ the following comm-
		and is suggested:

			>>> export PYTHONPATH=/lrz/sys/graphics/visit/2.7.0/current/linux-x86_64/lib/site-packages/

		For more information refer to the VisIt user's guide in the following 
		link: 
			
			http://www.visitusers.org/index.php?title=Using_a_Standard_Python_Interpreter 

		ii) set TBSLAS_RESULT_DIR to the directory that contains the results 
		from simulation using tbslas. (.VTK files directory)

	4- run the script on the machine by invoking python.

			>>> python visualization.py


	IMPORTANT NOTE: make sure you are using the proper system on which Xlib is
	accessible by VisIt; this means you need to run the code on special nodes;
	namely Render Nodes. For instace, for linux cluster in LRZ one should should
	use the following command on the remote visualization nodes:

			>>> rvglrun python visualization.py
	
	For more information, please refer to the LRZ user manual web-page:

		https://www.lrz.de/services/v2c_en/remote_visualisation_en/super_muc_users_en/


###############################################################################
###############################################################################
#### USAGE='USAGE: python PROGRAM <num-of-nodes> <mpi-num-processes>' #########
###############################################################################
###############################################################################
"""
###############################################################################
########################### IMPORT SYSTEM LIBRARIES ###########################
import time
import sys
import os

###############################################################################
########################### SET THE TIME STRING ###############################
TIMESTR = time.strftime("%Y%m%d-%H%M%S")

###############################################################################
####################### CHECK THE ENVIRONMENTAL VARIABLES #####################
try:
    TBSLAS_RESULT_DIR = os.environ['TBSLAS_RESULT_DIR']
except KeyError as e:
    print "Environment variable {0} is not set.".format(e)
    sys.exit()

try:
    PYTHONPATH = os.environ['PYTHONPATH']
except KeyError as e:
    print "Environment variable {0} is not set.".format(e)
    sys.exit()

###############################################################################
########################### INPUT ARGUMENTS ###################################
USAGE='USAGE: python PROGRAM <num-of-nodes> <mpi-num-processes>'

if len(sys.argv) != 3:
    print USAGE
    sys.exit()

TOTAL_NUM_NODES         = int(sys.argv[1])
MPI_TOTAL_NUM_PORCESSES = int(sys.argv[2])


###############################################################################
####################### FIND THE PREFIX IN .VTK FILES #########################
fd=os.listdir(TBSLAS_RESULT_DIR)

if 'advection' in fd[0]:
	prefix='advection'
elif 'advdiff-ss' in fd[0]:
	prefix='advdiff-ss'
elif 'advdiff' in fd[0]:
	prefix='advdiff'
elif 'cubic' in fd[0]:
	prefix='cubic'
elif 'diffusion' in fd[0]:
	prefix='diffusion'
elif 'merge' in fd[0]:
	prefix='merge'
elif 'semilag-tree' in fd[0]:
	prefix='semilag-tree'

###############################################################################
########################## SET THE THE DIRECTORIES ############################
VTK_FILES=TBSLAS_RESULT_DIR+"/"+prefix+"_Vconc_T*_P.pvtu database"
IMAGE_DIR=TBSLAS_RESULT_DIR+"/images-"+TIMESTR
os.makedirs(IMAGE_DIR)


###############################################################################
############################### LAUNCH VISIT ##################################
#sys.path.append("/usr/local/2.9.2/linux-x86_64/lib/site-packages/")
from visit import*
AddArgument("-nn")
AddArgument(str(TOTAL_NUM_NODES))
AddArgument("-np")
AddArgument(str(MPI_TOTAL_NUM_PORCESSES))
AddArgument("-psub")
AddArgument("-pbatch")
AddArgument("-nowin")
Launch()
print "Engine is up and running"
#OpenComputeEngine("localhost", ("-nn", "1", "-np", "1", "-l", "psub", "-p", "pbatch"))

###############################################################################
######################### VISUALIZATION USING VISIT ###########################
#Open the database
visit.OpenDatabase(VTK_FILES, 0)

#Plot the contour
AddPlot("Contour", "cheb_value")
ContourAtts = ContourAttributes()
ContourAtts.contourMethod = ContourAtts.Value
ContourAtts.contourNLevels = 1
ContourAtts.contourValue = 0.5						#set the contour value
ContourAtts.SetMultiColor(0, (255, 0, 0, 255))
ContourAtts.wireframe = 0
ContourAtts.legendFlag = 0
SetPlotOptions(ContourAtts)

AddOperator("Box")
BoxAtts = BoxAttributes()
BoxAtts.minz = 0.5
BoxAtts.inverse = 0
SetOperatorOptions(BoxAtts, 1)

AddPlot("Pseudocolor", "cheb_value")
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.legendFlag = 0
PseudocolorAtts.minFlag = 1
PseudocolorAtts.maxFlag = 1
PseudocolorAtts.min = -0.5
PseudocolorAtts.max = 1.5
SetPlotOptions(PseudocolorAtts)

AddOperator("Box")
BoxAtts = BoxAttributes()
BoxAtts.minz = 0.5
BoxAtts.inverse = 1
SetOperatorOptions(BoxAtts, 1)

#Plot the mesh
AddPlot("Mesh", "mesh", 1, 1)
MeshAtts = MeshAttributes()
MeshAtts.legendFlag = 0
SetPlotOptions(MeshAtts)

#set the view attributes
v=GetView3D()
v.viewNormal = (0.3, 0.45, 0.85)
SetView3D(v)

AnnotationAtts = AnnotationAttributes()
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.databaseInfoFlag = 0
AnnotationAtts.timeInfoFlag = 0
AnnotationAtts.axes3D.visible = 0
AnnotationAtts.axes3D.triadFlag = 0
AnnotationAtts.axes3D.bboxFlag = 0
SetAnnotationAttributes(AnnotationAtts)

#Draw
DrawPlots()

#Set the window attributes
SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 0
SaveWindowAtts.outputDirectory = IMAGE_DIR
SaveWindowAtts.fileName = prefix
SaveWindowAtts.quality = 100
SaveWindowAtts.format = SaveWindowAtts.JPEG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
SaveWindowAtts.screenCapture = 0

#Traverse through states and save images
nts=TimeSliderGetNStates()
for ts in range(0,nts):
	TimeSliderSetState(ts)
	SetSaveWindowAttributes(SaveWindowAtts)
	SaveWindow()

sys.exit()

############################### END OF THE CODE ###############################
