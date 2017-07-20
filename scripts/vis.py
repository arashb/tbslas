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

#*************************************************************************
# In order to run this script, you should do the following in advance:
#
#         1- Make sure you have VisIt package on your machine; for LRZ Linux cluster,
#         one should load the VisIt module:
#
#                         >>> module load visit
#
#         2- run the script on the machine by invoking visit
#
#                         >>> visit -cli -nowin -s vis.py -i<vtk-files-dir>
#
#
#         IMPORTANT NOTE: make sure you are using the proper system on which Xlib is
#         accessible by VisIt; this means you need to run the code on special nodes;
#         namely Render Nodes. For instace, for linux cluster in LRZ one should should
#         use the following command on the remote visualization nodes:
#
#                         >>> rvglrun visit -cli -nowin -s vis.py -i<vtk-files-dir>
#
#         For more information, please refer to the LRZ user manual web-page:
#
#                 https://www.lrz.de/services/v2c_en/remote_visualisation_en/super_muc_users_en/
#*************************************************************************

############################################################################
# IMPORT SYSTEM LIBRARIES
############################################################################
import time
import sys
import os

############################################################################
# IMPORT LOCAL LIBRARIES
############################################################################
from visit import *

from vis_plot_utils        import *
from vis_plot_slice        import *
from vis_plot_porous       import *
from vis_plot_taylor_green import *
from vis_plot_two_vortex_tube import *

############################################################################
# INPUT ARGUMENTS
############################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input_dir', action='store')
args, unknown = parser.parse_known_args()

############################################################################
# SET THE TIME STRING
############################################################################
TIMESTR = time.strftime("%Y%m%d-%H%M%S")

############################################################################
# DATABASES
############################################################################
VTK_DIR   = args.input_dir
IMAGE_DIR = VTK_DIR+"/images-"+TIMESTR
os.makedirs(IMAGE_DIR)

CON_VTK_FILES    = VTK_DIR+"/"+"conc_T*_P.pvtu database"

CON_VTK_FILES1_0 = VTK_DIR+"/"+"conc01_T0000_P.pvtu"
CON_VTK_FILES2_0 = VTK_DIR+"/"+"conc02_T0000_P.pvtu"
CON_VTK_FILES3_0 = VTK_DIR+"/"+"conc03_T0000_P.pvtu"

CON_VTK_FILES1   = VTK_DIR+"/"+"conc01_T*_P.pvtu database"
CON_VTK_FILES2   = VTK_DIR+"/"+"conc02_T*_P.pvtu database"
CON_VTK_FILES3   = VTK_DIR+"/"+"conc03_T*_P.pvtu database"

RHO_VTK_FILES    = VTK_DIR+"/"+"stokes_rho_0_.pvtu"
VEL_VTK_FILES    = VTK_DIR+"/"+"stokes_vel_0_.pvtu"
VEL_VTK_FILES    = VTK_DIR+"/"+"vel_T*_P.pvtu database"
VOR_VTK_FILES    = VTK_DIR+"/"+"vort_T*_P.pvtu database"

## uncomment for taylor-green
#CON_VTK_FILES    = VTK_DIR+"/"+"conc_T*_P.pvtu database"
#VEL_VTK_FILES    = VTK_DIR+"/"+"velocity_T0000_P.pvtu"

############################################################################
# VISUALIZATION SCENARIOS
############################################################################
def vis_slice(vtk_files, output_dir):
    OpenDatabase(vtk_files)
    draw_slice()
    save_images(output_dir)

def vis_porous(rho_vtk_files, vel_vtk_files, conc_vtk_files, output_dir):
    OpenDatabase(rho_vtk_files, 0)
    draw_porous_media_IV()
    cut_porous_media()

    OpenDatabase(vel_vtk_files, 1)
    ActivateDatabase(vel_vtk_files)
    draw_porous_velocity()

    OpenDatabase(conc_vtk_files, 2)
    ActivateDatabase(conc_vtk_files)
    draw_concentration_field()

    set_view()
    save_images(output_dir)


def vis_porous_three_spheres(rho_vtk_files, vel_vtk_files, conc_vtk_files1, conc_vtk_files2, conc_vtk_files3, output_dir):
    OpenDatabase(rho_vtk_files, 0)
    draw_porous_media_IV()
    cut_porous_media()

    OpenDatabase(vel_vtk_files, 1)
    ActivateDatabase(vel_vtk_files)
    draw_porous_velocity()

    OpenDatabase(conc_vtk_files1, 2)
    ActivateDatabase(conc_vtk_files1)
    draw_three_concentration_fields(2, 'b')

    OpenDatabase(conc_vtk_files2, 3)
    ActivateDatabase(conc_vtk_files2)
    draw_three_concentration_fields(3, 'g')

    OpenDatabase(conc_vtk_files3, 4)
    ActivateDatabase(conc_vtk_files3)
    draw_three_concentration_fields(4, 'y')

    set_view()
    save_images(output_dir)

def vis_porous_three_spheres_initial_camera_rotation(rho_vtk_files, vel_vtk_files, conc_vtk_files1, conc_vtk_files2, conc_vtk_files3, output_dir):
    OpenDatabase(rho_vtk_files, 0)
    SetActivePlots(0)
    draw_porous_media_IV()
    cut_porous_media()

    SetActivePlots(1)
    draw_porous_media_IV()
    cut_porous_media(1)
    translate_porous()

    OpenDatabase(vel_vtk_files, 0)
    ActivateDatabase(vel_vtk_files)
    draw_porous_velocity(2,0)

    OpenDatabase(vel_vtk_files, 0)
    ActivateDatabase(vel_vtk_files)
    draw_porous_velocity(3,1)
    SetActivePlots(3)
    translate_porous()

    OpenDatabase(conc_vtk_files2, 0)
    ActivateDatabase(conc_vtk_files1)
    draw_three_concentration_fields(4, 'b')

    OpenDatabase(conc_vtk_files2, 0)
    ActivateDatabase(conc_vtk_files2)
    draw_three_concentration_fields(5, 'g')

    OpenDatabase(conc_vtk_files3, 0)
    ActivateDatabase(conc_vtk_files3)
    draw_three_concentration_fields(6, 'y')

    change_view_and_save(output_dir)

    ToggleLockViewMode()
    ToggleMaintainViewMode()
    translate_and_save(output_dir, 1, 3)

def vis_taylor_green(vel_vtk_files, conc_vtk_files, output_dir):
    OpenDatabase(vel_vtk_files, 0)
    draw_taylor_green_velocity(1,0)
    # OpenDatabase(conc_vtk_files, 0)
    # draw_taylor_green_concentration_field(1)
    set_view(8*pi/12)
    save_images(output_dir)

def vis_two_vortex_tubes(vor_vtk_files, output_dir):
    OpenDatabase(vor_vtk_files, 0)
    draw_two_vortex_vorticity(1,0);
    # draw_taylor_green_velocity(1,0)
    # set_view(8*pi/12)
    save_images(output_dir)

############################################################################
# MAIN
############################################################################
if __name__ == '__main__':
    ########################################################################
    # PLOTS
    ########################################################################
    # vis_slice(CON_VTK_FILES, IMAGE_DIR)

    #vis_porous(RHO_VTK_FILES, VEL_VTK_FILES, CON_VTK_FILES, IMAGE_DIR)
    #vis_porous_three_spheres(RHO_VTK_FILES, VEL_VTK_FILES, CON_VTK_FILES1, CON_VTK_FILES2, CON_VTK_FILES3, IMAGE_DIR)
    #vis_porous_three_spheres_initial_camera_rotation(RHO_VTK_FILES, VEL_VTK_FILES, CON_VTK_FILES1_0, CON_VTK_FILES2_0, CON_VTK_FILES3_0, IMAGE_DIR)

    # vis_taylor_green(VEL_VTK_FILES ,CON_VTK_FILES, IMAGE_DIR)
    vis_taylor_green(VOR_VTK_FILES ,CON_VTK_FILES, IMAGE_DIR)

    # vis_two_vortex_tubes(VOR_VTK_FILES, IMAGE_DIR)

    sys.exit()
