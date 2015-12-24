"""
In order to run this script, you should do the following in advance:

        1- Make sure you have VisIt package on your machine; for LRZ Linux cluster,
        one should load the VisIt module:

                        >>> module load visit

        2- run the script on the machine by invoking visit

                        >>> visit -cli -nowin -s vis.py <vtk-files-dir>


        IMPORTANT NOTE: make sure you are using the proper system on which Xlib is
        accessible by VisIt; this means you need to run the code on special nodes;
        namely Render Nodes. For instace, for linux cluster in LRZ one should should
        use the following command on the remote visualization nodes:

                        >>> rvglrun visit -cli -nowin -s vis.py <vtk-files-dir>

        For more information, please refer to the LRZ user manual web-page:

                https://www.lrz.de/services/v2c_en/remote_visualisation_en/super_muc_users_en/
"""
###############################################################################
# IMPORT SYSTEM LIBRARIES
###############################################################################
import time
import sys
import os

###############################################################################
# IMPORT SYSTEM LIBRARIES
###############################################################################
from visit import *

###############################################################################
# SET THE TIME STRING
###############################################################################
TIMESTR = time.strftime("%Y%m%d-%H%M%S")

###############################################################################
# INPUT ARGUMENTS
###############################################################################
import optparse
parser = optparse.OptionParser()
parser.add_option('-i', '--input',
    action="store", dest="input_dir",
    help="VTK files dir")

options, args = parser.parse_args()

VTK_DIR = options.input_dir

###############################################################################
# FIND THE PREFIX IN .VTK FILES
###############################################################################
fd=os.listdir(VTK_DIR)

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
VTK_FILES=VTK_DIR+"/"+prefix+"_Vconc_T*_P.pvtu database"
IMAGE_DIR=VTK_DIR+"/images-"+TIMESTR
os.makedirs(IMAGE_DIR)


def draw_porous_media():
    #Plot the contour
    AddPlot("Contour", "cheb_value")
    ContourAtts = ContourAttributes()
    ContourAtts.contourMethod = ContourAtts.Value
    ContourAtts.contourNLevels = 1
    ContourAtts.contourValue = 0.5                                          #set the contour value
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

def draw_slice():
    AddPlot("Mesh", "mesh", 1, 1)
    AddPlot("Pseudocolor", "cheb_value", 1, 1)
    AddOperator("Slice", 1)
    SetActivePlots(1)
    SetActivePlots(1)
    SliceAtts = SliceAttributes()
    SliceAtts.originType = SliceAtts.Intercept  # Point, Intercept, Percent, Zone, Node
    SliceAtts.originPoint = (0, 0, 0)
    SliceAtts.originIntercept = 0.5
    SliceAtts.originPercent = 0
    SliceAtts.originZone = 0
    SliceAtts.originNode = 0
    SliceAtts.normal = (0, 0, 1)
    SliceAtts.axisType = SliceAtts.ZAxis  # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
    SliceAtts.upAxis = (0, 1, 0)
    SliceAtts.project2d = 1
    SliceAtts.interactive = 1
    SliceAtts.flip = 0
    SliceAtts.originZoneDomain = 0
    SliceAtts.originNodeDomain = 0
    SliceAtts.meshName = "mesh"
    SliceAtts.theta = 0
    SliceAtts.phi = 90
    SetOperatorOptions(SliceAtts, 1)

    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.timeInfoFlag = 0
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.triadFlag = 0
    AnnotationAtts.axes3D.bboxFlag = 0
    SetAnnotationAttributes(AnnotationAtts)

    DrawPlots()

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':
    OpenDatabase(VTK_FILES, 0)

    draw_slice()
    # draw_porous_media()

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
