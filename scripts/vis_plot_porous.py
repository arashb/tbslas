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

from vis_plot_utils import *

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
    ResampleAtts = ResampleAttributes()
    ResampleAtts.useExtents = 1
    ResampleAtts.is3D = 1
    ResampleAtts.samplesX = 200
    ResampleAtts.samplesY = 200
    ResampleAtts.samplesZ = 200
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


def cut_porous_media(i=0):

    AddOperator("Clip", 0)
    ClipAtts = ClipAttributes()
    ClipAtts.plane1Status = 1
    ClipAtts.plane2Status = 0
    ClipAtts.plane3Status = 0

    if (i==0):
        ClipAtts.plane1Origin = (0.5, 0.5, 0.5)
        ClipAtts.plane1Normal = (0, 1, 1)
        ClipAtts.planeInverse = 0
    elif (i==1):
        ClipAtts.plane1Origin = (0.5, 0.5, 0.5)
        ClipAtts.plane1Normal = (0, 1, 1)
        ClipAtts.planeInverse = 1

    else:
        ClipAtts.plane1Origin = (0.7, 0.5, 0.5)
        ClipAtts.plane1Normal = (-1, 0, 0)
        ClipAtts.planeInverse = 0


    ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
    SetOperatorOptions(ClipAtts, 0)

    DrawPlots()

def translate_porous():

    AddOperator("Transform", 0)

    TransformAtts = TransformAttributes()
    TransformAtts.doTranslate = 1
    TransformAtts.translateX = 0
    TransformAtts.translateY = 0
    TransformAtts.translateZ = 0

    SetOperatorOptions(TransformAtts, 0)

    DrawPlots()

def translate_and_save(image_dir, material_plot_number=0, velocity_plot_number=1):
    
    d_x=0.01
    d_y=0.01
    d_z=0.0
    r=int(floor((2.0)/d_x))

    x=0
    y=0
    z=0
    for i in range(r):

	x=x+d_x
	y=y+d_y
	z=z+d_z

	SetActivePlots(material_plot_number)
	TransformAtts = TransformAttributes()
	TransformAtts.doTranslate = 1
	TransformAtts.translateX = x
	TransformAtts.translateY = y
	TransformAtts.translateZ = z
	SetOperatorOptions(TransformAtts, 0)

	SetActivePlots(velocity_plot_number)
	TransformAtts = TransformAttributes()
	TransformAtts.doTranslate = 1
	TransformAtts.translateX = x
	TransformAtts.translateY = y
	TransformAtts.translateZ = z
	SetOperatorOptions(TransformAtts, 0)

	StreamlineAtts = StreamlineAttributes()
	StreamlineAtts.sourceType = StreamlineAtts.SpecifiedBox
	StreamlineAtts.boxExtents = (0.9+x, 1+x, 0+y, 1+y, 0+z, 1+z)
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

        ClipAtts = ClipAttributes()
	ClipAtts.plane1Status = 1
        ClipAtts.plane2Status = 0
        ClipAtts.plane3Status = 0
        ClipAtts.plane1Origin = (0.5+x, 0.5+y, 0.5+z)
        ClipAtts.plane1Normal = (0, 1, 1)
        ClipAtts.planeInverse = 1
        ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
        SetOperatorOptions(ClipAtts, 0)

	save_images(image_dir)


def draw_porous_velocity(plot_number=2,i=0):

    AddPlot("Streamline", "cheb_value", 0, 0)
    StreamlineAtts = StreamlineAttributes()

    SetActivePlots(plot_number)
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

    if (i==0):
        SetActivePlots(plot_number)
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

    elif (i==1):
        SetActivePlots(plot_number)
        AddOperator("Clip", 0)
        ClipAtts = ClipAttributes()
        ClipAtts.plane1Status = 1
        ClipAtts.plane2Status = 0
        ClipAtts.plane3Status = 0
        ClipAtts.plane1Origin = (0.5, 0.5, 0.5)
        ClipAtts.plane1Normal = (0, 1, 1)
        ClipAtts.planeInverse = 1
        ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1
        SetOperatorOptions(ClipAtts, 0)

    DrawPlots()


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

def draw_three_concentration_fields(plot_number=2,color='b'):
  
     AddPlot("Volume", "cheb_value", 0, 0)

     SetActivePlots(plot_number)
     VolumeAtts = VolumeAttributes()
     VolumeAtts.legendFlag = 0
     VolumeAtts.lightingFlag = 0

     VolumeAtts.colorControlPoints.ClearControlPoints()
     p = ColorControlPoint()
     for j in range (0, 9):
             VolumeAtts.colorControlPoints.AddControlPoints(p)

     if color=='b':
          VolumeAtts.colorControlPoints.GetControlPoints(0).colors = (8, 64, 129, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(0).position = 0
          VolumeAtts.colorControlPoints.GetControlPoints(1).colors = (8, 104, 172, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(1).position = 0.125
          VolumeAtts.colorControlPoints.GetControlPoints(2).colors = (43, 140, 190, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(2).position = 0.25
          VolumeAtts.colorControlPoints.GetControlPoints(3).colors = (78, 179, 211, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(3).position = 0.325
          VolumeAtts.colorControlPoints.GetControlPoints(4).colors = (123, 204, 196, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(4).position = 0.5
          VolumeAtts.colorControlPoints.GetControlPoints(5).colors = (168, 221, 181, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(5).position = 0.675
          VolumeAtts.colorControlPoints.GetControlPoints(6).colors = (204, 235, 197, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(6).position = 0.75
          VolumeAtts.colorControlPoints.GetControlPoints(7).colors = (224, 243, 219, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(7).position = 0.875
          VolumeAtts.colorControlPoints.GetControlPoints(8).colors = (247, 252, 240, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(8).position = 1

     elif color=='g':
          VolumeAtts.colorControlPoints.GetControlPoints(0).colors = (0, 68, 27, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(0).position = 0
          VolumeAtts.colorControlPoints.GetControlPoints(1).colors = (0, 109, 44, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(1).position = 0.125
          VolumeAtts.colorControlPoints.GetControlPoints(2).colors = (35, 139, 69, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(2).position = 0.25
          VolumeAtts.colorControlPoints.GetControlPoints(3).colors = (65, 171, 93, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(3).position = 0.325
          VolumeAtts.colorControlPoints.GetControlPoints(4).colors = (116, 196, 118, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(4).position = 0.5
          VolumeAtts.colorControlPoints.GetControlPoints(5).colors = (161, 217, 155, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(5).position = 0.675
          VolumeAtts.colorControlPoints.GetControlPoints(6).colors = (199, 233, 192, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(6).position = 0.75
          VolumeAtts.colorControlPoints.GetControlPoints(7).colors = (229, 245, 224, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(7).position = 0.875
          VolumeAtts.colorControlPoints.GetControlPoints(8).colors = (247, 252, 245, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(8).position = 1

     elif color=='y':
          VolumeAtts.colorControlPoints.GetControlPoints(0).colors = (128, 128, 0, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(0).position = 0
          VolumeAtts.colorControlPoints.GetControlPoints(1).colors = (160, 160, 32, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(1).position = 0.125
          VolumeAtts.colorControlPoints.GetControlPoints(2).colors = (192, 192, 64, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(2).position = 0.25
          VolumeAtts.colorControlPoints.GetControlPoints(3).colors = (223, 223, 96, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(3).position = 0.325
          VolumeAtts.colorControlPoints.GetControlPoints(4).colors = (255, 255, 128, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(4).position = 0.5
          VolumeAtts.colorControlPoints.GetControlPoints(5).colors = (255, 255, 160, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(5).position = 0.675
          VolumeAtts.colorControlPoints.GetControlPoints(6).colors = (255, 255, 192, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(6).position = 0.75
          VolumeAtts.colorControlPoints.GetControlPoints(7).colors = (255, 255, 223, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(7).position = 0.875
          VolumeAtts.colorControlPoints.GetControlPoints(8).colors = (255, 255, 255, 255)
          VolumeAtts.colorControlPoints.GetControlPoints(8).position = 1


     VolumeAtts.colorControlPoints.smoothing = VolumeAtts.colorControlPoints.Linear  # None, Linear, CubicSpline
     VolumeAtts.colorControlPoints.equalSpacingFlag = 1
     VolumeAtts.colorControlPoints.discreteFlag = 0
     VolumeAtts.opacityAttenuation = 1
     VolumeAtts.opacityMode = VolumeAtts.FreeformMode  # FreeformMode, GaussianMode, ColorTableMode
     VolumeAtts.resampleFlag = 1
     VolumeAtts.resampleTarget = 5000000
     VolumeAtts.useColorVarMin = 1
     VolumeAtts.colorVarMin = 0.3
     VolumeAtts.useColorVarMax = 1
     VolumeAtts.colorVarMax = 1.5
     VolumeAtts.rendererType = VolumeAtts.Splatting  # Splatting, Texture3D, RayCasting, RayCastingIntegration, SLIVR, RayCastingSLIVR, Tuvok
     VolumeAtts.gradientType = VolumeAtts.SobelOperator  # CenteredDifferences, SobelOperator
     VolumeAtts.scaling = VolumeAtts.Skew  # Linear, Log, Skew
     VolumeAtts.skewFactor = 0.001
     VolumeAtts.lowGradientLightingReduction = VolumeAtts.Lowest  # Off, Lowest, Lower, Low, Medium, High, Higher, Highest
     SetPlotOptions(VolumeAtts)

     DrawPlots()

