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
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input_dir', action='store')
args, unknown = parser.parse_known_args()
VTK_DIR = args.input_dir

###############################################################################
# FIND THE PREFIX IN .VTK FILES
###############################################################################
fd=os.listdir(VTK_DIR)

for item in fd:

    if 'advection' in item:
        prefix='advection'
    elif 'advdiff-ss' in item:
        prefix='advdiff-ss'
    elif 'advdiff' in item:
        prefix='advdiff'
    elif 'cubic' in item:
        prefix='cubic'
    elif 'diffusion' in item:
        prefix='diffusion'
    elif 'merge' in item:
        prefix='merge'
    elif 'semilag-tree' in item:
        prefix='semilag-tree'

###############################################################################
########################## SET THE THE DIRECTORIES ############################
CONC_VTK_FILES=VTK_DIR+"/"+prefix+"_Vconc_T*_P.pvtu database"
RHO_VTK_FILES =VTK_DIR+"/"+"stokes_rho_0_.pvtu"
VEL_VTK_FILES =VTK_DIR+"/"+"stokes_vel_0_.pvtu"

IMAGE_DIR=VTK_DIR+"/images-"+TIMESTR
os.makedirs(IMAGE_DIR)


def draw_porous_media():
    #Plot the contour
    AddPlot("Contour", "cheb_value")
    ContourAtts = ContourAttributes()
    ContourAtts.contourMethod = ContourAtts.Value
    ContourAtts.contourNLevels = 1
    ContourAtts.contourValue = 1e+9                                          #set the contour value
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

    #Draw
    DrawPlots()

def draw_porous_media_IV():

    AddPlot("Pseudocolor", "cheb_value", 0, 0)
    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.colorTableName = "OrRd"
    PseudocolorAtts.opacityType = PseudocolorAtts.Constant
    PseudocolorAtts.opacity = 1
    PseudocolorAtts.legendFlag = 0
    SetPlotOptions(PseudocolorAtts)

    AddOperator("Resample", 0)
    SetActivePlots(0)
    ResampleAtts = ResampleAttributes()
    ResampleAtts.useExtents = 1
    ResampleAtts.is3D = 1
    ResampleAtts.samplesX = 100
    ResampleAtts.samplesY = 100
    ResampleAtts.samplesZ = 100
    ResampleAtts.tieResolver = ResampleAtts.random
    ResampleAtts.distributedResample = 0
    SetOperatorOptions(ResampleAtts, 1)

    AddOperator("Isovolume", 0)
    IsovolumeAtts = IsovolumeAttributes()
    IsovolumeAtts.lbound = 1e+9
    IsovolumeAtts.ubound = 1e+37
    IsovolumeAtts.variable = "default"
    SetOperatorOptions(IsovolumeAtts, 1)

    DrawPlots()


def cut_media(i=0):

    AddOperator("Clip", 0)
    SetActivePlots(0)
    ClipAtts = ClipAttributes()
    ClipAtts.plane1Status = 1
    ClipAtts.plane2Status = 0
    ClipAtts.plane3Status = 0

    if (i==0):
	ClipAtts.plane1Origin = (0.5, 0.5, 0.5)
	ClipAtts.plane1Normal = (0, 1, 1)    
    else:
	ClipAtts.plane1Origin = (0.7, 0.5, 0.5)
	ClipAtts.plane1Normal = (-1, 0, 0)

    ClipAtts.planeInverse = 0
    ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
    SetOperatorOptions(ClipAtts, 0)

    DrawPlots()


def draw_slice():
    AddPlot("Mesh", "mesh", 1, 1)
    
    AddOperator("Slice", 1)
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
    SliceAtts.interactive = 0
    SliceAtts.flip = 0
    SliceAtts.originZoneDomain = 0
    SliceAtts.originNodeDomain = 0
    SliceAtts.meshName = "mesh"
    SliceAtts.theta = 0
    SliceAtts.phi = 90
    SetOperatorOptions(SliceAtts, 1)

    AddPlot("Pseudocolor", "cheb_value", 1, 1)

    DrawPlots()


def draw_velocity(i=0):

    OpenDatabase(VTK_FILES, 1)
    ActivateDatabase(VTK_FILES)

    AddPlot("Streamline", "cheb_value", 0, 0)
    StreamlineAtts = StreamlineAttributes()

    SetActivePlots(1)
    StreamlineAtts.sourceType = StreamlineAtts.SpecifiedBox
    StreamlineAtts.boxExtents = (0.9, 1, 0, 1, 0, 1)
    StreamlineAtts.useWholeBox = 0
    StreamlineAtts.coloringMethod = StreamlineAtts.ColorBySpeed
    StreamlineAtts.colorTableName = "Default"
    StreamlineAtts.legendFlag = 0
    StreamlineAtts.lightingFlag = 0
    StreamlineAtts.legendMinFlag = 0
    StreamlineAtts.legendMaxFlag = 1
    StreamlineAtts.legendMin = 0
    StreamlineAtts.legendMax = 0.05
    StreamlineAtts.showSeeds = 0
    StreamlineAtts.fillInterior = 1
    StreamlineAtts.randomSamples = 1
    StreamlineAtts.randomSeed = 0
    StreamlineAtts.numberOfRandomSamples = 2000
    SetPlotOptions(StreamlineAtts)

    if (i==1):
	SetActivePlots(1)
	AddOperator("Clip", 0)
	ClipAtts = ClipAttributes()
	ClipAtts.plane1Status = 1
	ClipAtts.plane2Status = 0
	ClipAtts.plane3Status = 0
	ClipAtts.plane1Origin = (0.5, 0.5, 0.5)
	ClipAtts.plane1Normal = (0, 1, 1)
	ClipAtts.planeInverse = 0
	ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
	SetOperatorOptions(ClipAtts, 0)

    DrawPlots()


def set_view():

    #set the view attributes
    v=GetView3D()
    #v.viewNormal = (0.3, 0.45, 0.85)
    v.viewNormal = (-7, 5, 4)
    v.viewUp = (5, 8 , -2)
    SetView3D(v)


def save_images():

    #set annotation attributes
    AnnotationAtts = AnnotationAttributes()
    AnnotationAtts.userInfoFlag = 0
    AnnotationAtts.databaseInfoFlag = 0
    AnnotationAtts.timeInfoFlag = 0
    AnnotationAtts.axes3D.visible = 0
    AnnotationAtts.axes3D.triadFlag = 0
    AnnotationAtts.axes3D.bboxFlag = 0
    SetAnnotationAttributes(AnnotationAtts)

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

def draw_concentration_field():

    AddPlot("Volume", "cheb_value", 0, 0)
    VolumeAtts = VolumeAttributes()
    VolumeAtts.legendFlag = 0
    VolumeAtts.lightingFlag = 0

    #adding control points for color setting
    p = ColorControlPoint()
    VolumeAtts.colorControlPoints.AddControlPoints(p)
    p.colors = (208, 209, 230, 255)
    p.position = 0
    VolumeAtts.colorControlPoints.AddControlPoints(p)
    p.colors = (166, 189, 219, 255)
    p.position = 0.2
    VolumeAtts.colorControlPoints.AddControlPoints(p)
    p.colors = (116, 169, 207, 255)
    p.position = 0.4
    VolumeAtts.colorControlPoints.AddControlPoints(p)
    p.colors = (54, 144, 192, 255)
    p.position = 0.6
    VolumeAtts.colorControlPoints.AddControlPoints(p)
    p.colors = (5, 112, 176, 255)
    p.position = 0.75
    VolumeAtts.colorControlPoints.AddControlPoints(p)
    p.colors = (4, 90, 141, 255)
    p.position = 0.875
    VolumeAtts.colorControlPoints.AddControlPoints(p)
    p.colors = (2, 56, 88, 255)
    p.position = 1
    VolumeAtts.colorControlPoints.AddControlPoints(p)
    VolumeAtts.colorControlPoints.smoothing = VolumeAtts.colorControlPoints.Linear
    VolumeAtts.opacityAttenuation = 0.5
    VolumeAtts.opacityMode = VolumeAtts.FreeformMode

    VolumeAtts.resampleFlag = 1
    VolumeAtts.useColorVarMin = 1
    VolumeAtts.colorVarMin = 0.1
    VolumeAtts.useColorVarMax = 1
    VolumeAtts.colorVarMax = 0.2
    VolumeAtts.samplesPerRay = 1000
    VolumeAtts.rendererType = VolumeAtts.RayCasting
    VolumeAtts.gradientType = VolumeAtts.SobelOperator  # CenteredDifferences, SobelOperator
    VolumeAtts.num3DSlices = 200
    VolumeAtts.limitsMode = VolumeAtts.OriginalData
    VolumeAtts.rendererSamples = 3
    VolumeAtts.transferFunctionDim = 1
    VolumeAtts.lowGradientLightingReduction = VolumeAtts.Lowest  # Off, Lowest, Lower, Low, Medium, High, Higher, Highest
    SetPlotOptions(VolumeAtts) 
    
    DrawPlots()

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':

#     OpenDatabase(RHO_VTK_FILES, 0)
#     OpenDatabase(VEL_VTK_FILES, 1)
    OpenDatabase(CONC_VTK_FILES, 2)

#     draw_porous_media()
#     draw_porous_media_IV()
#     cut_media()

#     ActivateDatabase(CONC_VTK_FILES)
#     draw_velocity()

#     ActivateDatabase(CONC_VTK_FILES)
#     draw_concentration_field()
    draw_slice()

#     set_view()
    save_images()

    sys.exit()
