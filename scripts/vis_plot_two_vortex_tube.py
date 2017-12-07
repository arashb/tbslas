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

# def draw_taylor_green_concentration_field(plot_number):
#     AddPlot("Contour", "cheb_value", 0, 0)
#     SetActivePlots(plot_number)
#     ContourAtts = ContourAttributes()
#     ContourAtts.defaultPalette.equalSpacingFlag = 1
#     ContourAtts.defaultPalette.discreteFlag = 1
#     ContourAtts.colorType = ContourAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
#     ContourAtts.singleColor = (51, 204, 204, 255)
#     ContourAtts.contourValue = (0.5)
#     ContourAtts.contourMethod = ContourAtts.Value  # Level, Value, Percent
#     ContourAtts.legendFlag = 0
#     SetPlotOptions(ContourAtts)
#     DrawPlots()

def draw_two_vortex_vorticity(plot_number,i):
    AddPlot("Contour", "cheb_value_magnitude", 0, 0)
    SetActivePlots(plot_number)
    ContourAtts = ContourAttributes()
    ContourAtts.colorType = ContourAtts.ColorByColorTable  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
    ContourAtts.colorTableName = "OrRd"
    ContourAtts.legendFlag = 1
    ContourAtts.contourNLevels = 5
    ContourAtts.contourPercent = (30)
    ContourAtts.contourMethod = ContourAtts.Percent  # Level, Value, Percent
    # ContourAtts.contourMethod = ContourAtts.Level  # Level, Value, Percent
    SetPlotOptions(ContourAtts)

    View3DAtts = View3DAttributes()
    View3DAtts.viewNormal = (-0.26083, 0.406859, 0.875462)
    View3DAtts.focus = (0.5, 0.5, 0.5)
    View3DAtts.viewUp = (0.283676, 0.899119, -0.333337)
    View3DAtts.viewAngle = 30
    View3DAtts.parallelScale = 0.866025
    View3DAtts.nearPlane = -1.73205
    View3DAtts.farPlane = 1.73205
    View3DAtts.imagePan = (0, 0)
    View3DAtts.imageZoom = 1
    View3DAtts.perspective = 1
    View3DAtts.eyeAngle = 2
    View3DAtts.centerOfRotationSet = 0
    View3DAtts.centerOfRotation = (0.5, 0.5, 0.5)
    View3DAtts.axis3DScaleFlag = 0
    View3DAtts.axis3DScales = (1, 1, 1)
    View3DAtts.shear = (0, 0, 1)
    View3DAtts.windowValid = 1
    SetView3D(View3DAtts)

    DrawPlots()
