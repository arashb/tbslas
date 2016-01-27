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
from math  import *
from visit import *

def set_view(theta=7*pi/12):
    
    phi = pi/4
    #set the view attributes
    v=GetView3D()
    v.imageZoom = 1.0
    v.viewNormal = (cos(theta), cos(phi)*sin(theta), sin(phi)*sin(theta))
    v.viewUp = (0, 0, 1)
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


def change_view_and_save(image_dir, theta_i=pi/4, theta_f=7*pi/12):

    phi = pi/4
    d_theta=0.01
    r=int(floor((theta_f-theta_i)/d_theta))+1
    v=GetView3D()

    theta=theta_i
    for i in range(r):
    #set the view attributes
	v.imageZoom = 0.8+i*0.2/r
        v.viewNormal = (cos(theta), cos(phi)*sin(theta), sin(phi)*sin(theta))
        v.viewUp = (0, 0, 1)
        SetView3D(v)
	theta=theta+d_theta
	save_images(image_dir)

