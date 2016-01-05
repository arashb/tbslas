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


def cut_porous_media(i=0):

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


def draw_porous_velocity(i=0):

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

    if (i==0):
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
