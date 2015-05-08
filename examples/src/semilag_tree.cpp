#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <string>

#include <pvfmm_common.hpp>
#include <mpi_tree.hpp>
#include <cheb_node.hpp>
#include <utils.hpp>
#include <vector.hpp>
#include <cheb_utils.hpp>
#include <profile.hpp>

#include <utils/common.h>
#include <tree/semilag_tree.h>
#include <tree/utils_tree.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;

double tcurr = 0;

const char* OUTPUT_FILE_FORMAT = "%s/%s-VAR_%s-TS_%04d-RNK";
const char* OUTPUT_FILE_PREFIX = "semilag-tree";

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  parse_command_line_options(argc, argv);

  {
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    pvfmm::Profile::Enable(true);

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();
    sim_config->vtk_filename_format     = OUTPUT_FILE_FORMAT;
    sim_config->vtk_filename_prefix     = OUTPUT_FILE_PREFIX;
    sim_config->vtk_filename_variable   = "conc";
    // =========================================================================
    // PRINT METADATA
    // =========================================================================
    if (!myrank) {
      MetaData_t::Print();
    }
    // =========================================================================
    // INIT THE TREES
    // =========================================================================
    Tree_t tvel_curr(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  tbslas::get_vorticity_field<double,3>,
                                  3,
                                  tvel_curr);
    tvel_curr.ConstructLET(pvfmm::FreeSpace);

    char out_name_buffer[300];
    if (sim_config->vtk_save) {
      snprintf(out_name_buffer,
               sizeof(out_name_buffer),
               sim_config->vtk_filename_format.c_str(),
               tbslas::get_result_dir().c_str(),
               sim_config->vtk_filename_prefix.c_str(),
               "vel",
               0);
      tvel_curr.Write2File(out_name_buffer, sim_config->vtk_order);
    }
    Tree_t tconc_curr(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  get_gaussian_field_cylinder_atT<double,3>,
                                  1,
                                  tconc_curr);
    if (sim_config->vtk_save) {
      snprintf(out_name_buffer,
               sizeof(out_name_buffer),
               sim_config->vtk_filename_format.c_str(),
               tbslas::get_result_dir().c_str(),
               sim_config->vtk_filename_prefix.c_str(),
               sim_config->vtk_filename_variable.c_str(),
               0);
      tconc_curr.Write2File(out_name_buffer, sim_config->vtk_order);
    }
    if (sim_config->dt == 0) {
      double vel_max_value;
      int vel_max_depth;
      tbslas::GetMaxTreeValues<Tree_t>
          (tvel_curr, vel_max_value, vel_max_depth);

      if (!myrank)
        printf("%d: VEL MAX VALUE: %f VEL MAX DEPTH:%d\n",
               myrank,
               vel_max_value,
               vel_max_depth);

      double conc_max_value;
      int conc_max_depth;
      tbslas::GetMaxTreeValues<Tree_t>
          (tconc_curr, conc_max_value, conc_max_depth);

      if (!myrank)
        printf("%d:CON MAX VALUE: %f CON MAX DEPTH:%d\n",
               myrank,
               conc_max_value,
               conc_max_depth);

      double cfl      = 1;
      double dx_min   = pow(0.5, conc_max_depth);
      sim_config->dt = (cfl * dx_min)/vel_max_value;
    }
    // =========================================================================
    // RUN
    // =========================================================================
    tbslas::RunSemilagSimulationInSitu(&tvel_curr,
                                       &tconc_curr,
                                       sim_config->total_num_timestep,
                                       sim_config->dt,
                                       sim_config->num_rk_step,
                                       true);

    //Output Profiling results.
    pvfmm::Profile::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
