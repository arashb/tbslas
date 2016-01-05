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
########################## SET THE THE DIRECTORIES ############################
CONC_VTK_FILES=VTK_DIR+"/"+"conc_T*_P.pvtu database"
RHO_VTK_FILES =VTK_DIR+"/"+"stokes_rho_0_.pvtu"
VEL_VTK_FILES =VTK_DIR+"/"+"stokes_vel_0_.pvtu"

IMAGE_DIR=VTK_DIR+"/images-"+TIMESTR
os.makedirs(IMAGE_DIR)

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
    SaveWindowAtts.fileName = "image_"
    SaveWindowAtts.quality = 100
    SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
    SaveWindowAtts.width = 2048
    SaveWindowAtts.height = 2048
    SaveWindowAtts.screenCapture = 0

    #Traverse through states and save images
    nts=TimeSliderGetNStates()
    for ts in range(0,nts):
        TimeSliderSetState(ts)
        SetSaveWindowAttributes(SaveWindowAtts)
        SaveWindow()


