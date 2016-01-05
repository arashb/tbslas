#*************************************************************************
# Copyright (C) 2015 by Arash Bakhtiari
# You may not use this file except in compliance with the License.
# You obtain a copy of the License in the LICENSE file.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#*************************************************************************

###############################################################################
# IMPORT LOCAL LIBRARIES
###############################################################################
from visit import *

def set_view():

    #set the view attributes
    v=GetView3D()
    #v.viewNormal = (0.3, 0.45, 0.85)
    v.viewNormal = (-7, 5, 4)
    v.viewUp = (5, 8 , -2)
    SetView3D(v)

def save_images(image_dir):

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
    SaveWindowAtts.outputDirectory = image_dir
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
