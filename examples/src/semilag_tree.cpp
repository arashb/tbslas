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

#include <utils/common.h>
// Enable profiling
#define __TBSLAS_PROFILE__ 5
#include <utils/profile.h>
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

  // =========================================================================
  // COMMAND LINE OPTIONS
  // =========================================================================
  commandline_option_start(argc, argv);

  omp_set_num_threads(
      atoi(commandline_option(argc, argv, "-omp", "1", false,
                              "-omp  <int> = (1)    : Number of OpenMP threads.")));

  size_t N =
      (size_t)strtod(
          commandline_option(argc, argv, "-N", "1", true,
                             "-N    <int>          : Number of point sources."),
          NULL);

  size_t M =
      (size_t)strtod(
          commandline_option(argc, argv, "-M", "1", false,
                             "-M    <int>          : Number of points per octant."),
          NULL);

  int q =
      strtoul(
          commandline_option(argc, argv, "-q", "14", false,
                             "-q    <int> = (14)   : Chebyshev order (+ve integer)."),
          NULL,10);

  int d =
      strtoul(
          commandline_option(argc, argv, "-d", "15", false,
                             "-d    <int> = (15)   : Maximum tree depth."),
          NULL,10);

  double tol =
      strtod(
          commandline_option(argc, argv, "-tol", "1e-5", false,
                             "-tol <real> = (1e-5) : Tolerance for adaptive refinement.")
          ,NULL);
  bool adap =
      (commandline_option(argc, argv, "-adap", NULL, false,
                          "-adap                : Adaptive tree refinement." )!=NULL);

  int tn =
      strtoul(
          commandline_option(argc, argv, "-tn", "1", false,
                             "-tn   <int> = (1)    : Number of time steps."),
          NULL,10);

  bool cubic =
      (commandline_option(argc, argv, "-cubic", NULL, false,
                          "-cubic               : Cubic Interpolation  used to evaluate tree values.")!=NULL);

  double dt =
      strtod(
          commandline_option(argc, argv,  "-dt",  "0", false,
                             "-tol <real> = (1e-5) : Temporal resolution." ), NULL);

  commandline_option_end(argc, argv);

  {
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    tbslas::Profile<double>::Enable(true, &comm);

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    tbslas::SimConfig* sim_config     = tbslas::SimConfigSingleton::Instance();
    sim_config->cubic                 = cubic;
    sim_config->total_num_timestep    = tn;
    sim_config->dt                    = dt; //3.14/12;
    sim_config->num_rk_step           = 1;
    sim_config->vtk_filename_format   = OUTPUT_FILE_FORMAT;
    sim_config->vtk_filename_prefix   = OUTPUT_FILE_PREFIX;
    sim_config->vtk_filename_variable = "conc";
    sim_config->vtk_order             = q;

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
    tbslas::ConstructTree<Tree_t>(N, M, q, d, adap, tol, comm,
                                  tbslas::get_vorticity_field<double,3>,
                                  3,
                                  tvel_curr);
    tvel_curr.ConstructLET(pvfmm::FreeSpace);

    char out_name_buffer[300];
    snprintf(out_name_buffer,
             sizeof(out_name_buffer),
             sim_config->vtk_filename_format.c_str(),
             tbslas::get_result_dir().c_str(),
             sim_config->vtk_filename_prefix.c_str(),
             "vel",
             0);
    tvel_curr.Write2File(out_name_buffer, sim_config->vtk_order);

    Tree_t tconc_curr(comm);
    tbslas::ConstructTree<Tree_t>(N, M, q, d, adap, tol, comm,
                                  get_gaussian_field_cylinder_atT<double,3>,
                                  1,
                                  tconc_curr);
    snprintf(out_name_buffer,
             sizeof(out_name_buffer),
             sim_config->vtk_filename_format.c_str(),
             tbslas::get_result_dir().c_str(),
             sim_config->vtk_filename_prefix.c_str(),
             sim_config->vtk_filename_variable.c_str(),
             0);
    tconc_curr.Write2File(out_name_buffer, q);

    if (dt == 0) {
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
    // Tree_t* result;
    tbslas::RunSemilagSimulationInSitu(&tvel_curr,
                                       &tconc_curr,
                                       true,
                                       true);

    //Output Profiling results.
    tbslas::Profile<double>::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
