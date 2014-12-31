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

double tcurr = 0;

const char* OUTPUT_FILE_FORMAT = "%s/%s-VAR_%s-TS_%04d-RNK";
const char* OUTPUT_FILE_PREFIX = "zalesak";

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;

  // Read command line options.
  commandline_option_start(argc, argv);
  omp_set_num_threads( atoi(commandline_option(argc, argv,  "-omp",     "1", false, "-omp  <int> = (1)    : Number of OpenMP threads."          )));
  size_t   N = (size_t)strtod(commandline_option(argc, argv,    "-N",     "1",  true, "-N    <int>          : Number of point sources."           ),NULL);
  size_t   M = (size_t)strtod(commandline_option(argc, argv,    "-M",     "1", false, "-M    <int>          : Number of points per octant."       ),NULL);
  int      q =         strtoul(commandline_option(argc, argv,    "-q",    "14", false, "-q    <int> = (14)   : Chebyshev order (+ve integer)."     ),NULL,10);
  int      d =         strtoul(commandline_option(argc, argv,    "-d",    "15", false, "-d    <int> = (15)   : Maximum tree depth."                ),NULL,10);
  double tol =         strtod(commandline_option(argc, argv,  "-tol",  "1e-5", false, "-tol <real> = (1e-5) : Tolerance for adaptive refinement." ),NULL);
  bool  adap =               (commandline_option(argc, argv, "-adap",    NULL, false, "-adap                : Adaptive tree refinement."          )!=NULL);
  int     tn =         strtoul(commandline_option(argc, argv,    "-tn",    "1", false, "-tn    <int> = (1)   : Number of time steps."     ),NULL,10);
  double dt  =         strtod(commandline_option(argc, argv,  "-dt",  "0.1e-2", false, "-tol <real> = (1e-5) : Temporal resolution." ), NULL);

  commandline_option_end(argc, argv);

  {
    int myrank;
    MPI_Comm_rank(comm, &myrank);

    // tbslas::Profile<double>::Enable(true, &comm);

    // =========================================================================
    // INIT THE VELOCITY TREE
    // =========================================================================    
    Tree_t tvel_curr(comm);
    tbslas::ConstructTree<Tree_t>(N, M, q, d, adap, tol, comm,
                                  tbslas::get_vorticity_field<double,3>,
                                  3,
                                  tvel_curr);

    // char out_name_buffer[300];
    // snprintf(out_name_buffer, sizeof(out_name_buffer),
    //          "%s/sltree_vel_%d_", tbslas::get_result_dir().c_str(), 0);
    // tvel_curr.Write2File(out_name_buffer, q);

    // =========================================================================
    // INIT THE CONCENTRATION TREE
    // =========================================================================
    tcurr = 0;
    Tree_t tconc_curr(comm);
    tbslas::ConstructTree<Tree_t>(N, M, q, d, adap, tol, comm,
                                  get_gaussian_field_cylinder_atT<double,3>,
                                  1,
                                  tconc_curr);

    // clone tree
    Tree_t tconc_next(comm);
    tbslas::CloneTree<Tree_t>(tconc_curr, tconc_next, 1);

    // =========================================================================
    // PREPARE SIMULATION PARAMETERS
    // =========================================================================
    struct tbslas::SimParam<double> sim_param;
    sim_param.dt                 = dt;
    sim_param.total_num_timestep = tn;
    sim_param.num_rk_step        = 1;
    sim_param.vtk_filename_format = OUTPUT_FILE_FORMAT;
    sim_param.vtk_filename_prefix = OUTPUT_FILE_PREFIX;

    Tree_t* tresult;
    tbslas::RunSemilagSimulation(&tvel_curr,
                                 &tconc_curr,
                                 &tconc_next,
                                 &sim_param,
                                 &tresult,
                                 true,
                                 true);

    // =========================================================================
    // compute error
    // =========================================================================
    double l_inf_error, l_two_error;
    tcurr = tn*dt;
    tbslas::ComputeTreeError(*tresult,
                             get_gaussian_field_cylinder_atT<double,3>,
                             l_inf_error,
                             l_two_error);
    int tnln = tbslas::CountNumLeafNodes(*tresult);
    if(!myrank)
      printf("TOL: %2.10f LTWO: %2.10f LINF = %2.10f TNLN: %d\n", tol, l_two_error, l_inf_error, tnln);

    //Output Profiling results.
    // tbslas::Profile<double>::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
