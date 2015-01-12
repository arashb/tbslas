// *************************************************************************
// Copyright (C) 2014 by Arash Bakhtiari
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

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
#include <utils/metadata.h>
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
const char* OUTPUT_FILE_PREFIX = "cubic";

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;

  // Read command line options.
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

  double dt =
      strtod(
          commandline_option(argc, argv,  "-dt",  "0.1e-2", false,
                             "-tol <real> = (1e-5) : Temporal resolution." ), NULL);

  bool cubic =
      (commandline_option(argc, argv, "-cubic", NULL, false,
                          "-cubic               : Cubic Interpolation  used to evaluate tree values.")!=NULL);

  int cuf =
      strtoul(
          commandline_option(argc, argv, "-cuf", "4", false,
                             "-cuf   <int> = (4)    : Upsampling factor used for cubic interpolation."),
          NULL,10);

  bool cubic_analytical =
      (commandline_option(argc, argv, "-ca", NULL, false,
                          "-ca                  : Analytical values used in cubic interpolation upsampling.")!=NULL);

  commandline_option_end(argc, argv);

  {
    int myrank;
    MPI_Comm_rank(comm, &myrank);
    // tbslas::Profile<double>::Enable(true, &comm);

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();
    sim_config->total_num_timestep      = tn;
    sim_config->dt                      = dt;
    sim_config->num_rk_step             = 1;
    sim_config->vtk_filename_format     = OUTPUT_FILE_FORMAT;
    sim_config->vtk_filename_prefix     = OUTPUT_FILE_PREFIX;
    sim_config->vtk_filename_variable   = "conc";
    sim_config->vtk_order               = q;
    sim_config->use_cubic               = cubic;
    sim_config->cubic_upsampling_factor = cuf;
    sim_config->cubic_use_analytical    = cubic_analytical;
    // =========================================================================
    // PRINT METADATA
    // =========================================================================
    // if (!myrank) {
    //   MetaData_t::Print();
    // }
    // =========================================================================
    // CONSTRUCT TREE
    // =========================================================================
    tcurr = 0;
    Tree_t tree(comm);
    tbslas::ConstructTree<Tree_t>(N, M, q, d, adap, tol, comm,
                                  get_gaussian_field_cylinder_atT<double,3>,
                                  1,
                                  tree);

    // =========================================================================
    // compute error
    // =========================================================================
    std::vector<double> grid_points;
    // tbslas::CollectChebTreeGridPoints(tree, grid_points);

    std::vector<double> node_pos = pvfmm::cheb_nodes<double>(q, 3);

    // int grid_res = q*sim_config->cubic_upsampling_factor;

    // std::vector<double> node_pos(grid_res*grid_res*grid_res*3);
    // tbslas::get_reg_grid_points<double, 3>(grid_res,
    //                                        node_pos.data());

    tbslas::CollectTreeGridPoints(tree,
                                  node_pos,
                                  grid_points);

    double l_inf_error, l_two_error;
    tbslas::ComputeTreeError(tree,
                             get_gaussian_field_cylinder_atT<double,3>,
                             grid_points,
                             l_inf_error,
                             l_two_error);
    int tnln = tbslas::CountNumLeafNodes(tree);

    if(!myrank)
      printf("TOL: %2.5e LTWO: %2.5e LINF = %2.5e TNLN: %d\n",
             tol, l_two_error, l_inf_error, tnln);

    //Output Profiling results.
    // tbslas::Profile<double>::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
