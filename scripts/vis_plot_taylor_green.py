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

def draw_taylor_green_concentration_field(plot_number):
    AddPlot("Contour", "cheb_value", 0, 0)
    SetActivePlots(plot_number)
    ContourAtts = ContourAttributes()
    ContourAtts.defaultPalette.equalSpacingFlag = 1
    ContourAtts.defaultPalette.discreteFlag = 1
    ContourAtts.colorType = ContourAtts.ColorBySingleColor  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
    ContourAtts.singleColor = (51, 204, 204, 255)
    ContourAtts.contourValue = (0.5)
    ContourAtts.contourMethod = ContourAtts.Value  # Level, Value, Percent
    ContourAtts.legendFlag = 0
    SetPlotOptions(ContourAtts)
    DrawPlots()

def draw_taylor_green_velocity(plot_number,i):
    AddPlot("Contour", "cheb_value_magnitude", 0, 0)
    SetActivePlots(plot_number)
    ContourAtts = ContourAttributes()
    ContourAtts.colorType = ContourAtts.ColorByColorTable  # ColorBySingleColor, ColorByMultipleColors, ColorByColorTable
    ContourAtts.colorTableName = "OrRd"
    ContourAtts.legendFlag = 0
    ContourAtts.contourNLevels = 10
    ContourAtts.contourMethod = ContourAtts.Level  # Level, Value, Percent
    SetPlotOptions(ContourAtts)

    if (i==1):
	
	AddOperator("Clip", 1)
	ClipAtts = ClipAttributes()
	ClipAtts.funcType = ClipAtts.Plane  # Plane, Sphere
	ClipAtts.plane1Status = 1
	ClipAtts.plane2Status = 0
	ClipAtts.plane3Status = 0
	ClipAtts.plane1Origin = (0.5, 0.5, 0.5)
	ClipAtts.plane2Origin = (0, 0, 0)
	ClipAtts.plane3Origin = (0, 0, 0)
	ClipAtts.plane1Normal = (0, 1, 1)
	ClipAtts.plane2Normal = (0, 1, 0)
	ClipAtts.plane3Normal = (0, 0, 1)
	ClipAtts.planeInverse = 0
	ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1  # None, Plane1, Plane2, Plane3
	SetOperatorOptions(ClipAtts, 1)
	DrawPlots()
