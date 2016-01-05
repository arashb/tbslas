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
from vis_utils import *
from vis_slice_plots import *
from vis_porous_plots import *


def vis_slice(vtk_files):
    OpenDatabase(vtk_files)
    draw_slice()
    save_images()


def vis_porous(rho_vtk_files, vel_vtk_files, conc_vtk_files):

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
    save_images() 

################################################################################
# MAIN
################################################################################
if __name__ == '__main__':

    vis_slice(RHO_VTK_FILES)
    #vis_porous(RHO_VTK_FILES, VEL_VTK_FILES, CONC_VTK_FILES,)
    sys.exit()
