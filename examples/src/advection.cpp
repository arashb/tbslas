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
#include <profile.hpp>

#include <utils/common.h>
#include <utils/metadata.h>
#include <utils/reporter.h>
#include <tree/semilag_tree.h>
#include <tree/utils_tree.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;
double tcurr = 0;

void (*fn_vel)(const double* , int , double*)=NULL;
void (*fn_con)(const double* , int , double*)=NULL;


template <class Real_t>
void gaussian_kernel(const Real_t* coord,
             int n,
             Real_t* out) {
  const Real_t xc = 0.7;
  const Real_t yc = 0.7;
  const Real_t zc = 0.7;
  const int a  = -160;
  const Real_t amp = 1.0;
  tbslas::gaussian_kernel(coord,
              n,
              out,
              xc,
              yc,
              zc);
}

template <class Real_t>
void diffusion_kernel_atT(const Real_t* coord,
              int n,
              Real_t* out) {
  const Real_t TBSLAS_DIFF_COEFF = 0.0001;
  const Real_t amp = 1e-2;
  const Real_t xc = 0.5;
  const Real_t yc = 0.5;
  const Real_t zc = 0.5;
  tbslas::diffusion_kernel(coord,
                           n,
                           out,
                           TBSLAS_DIFF_COEFF,
                           tcurr+25,
                           amp,
                           xc,
                           yc,
                           zc);
}

template <class Real_t>
void get_hopf_field(const Real_t* coord,
            int n,
            Real_t* out) {
  const Real_t xc = 0.5;
  const Real_t yc = 0.5;
  const Real_t zc = 0.5;
  tbslas::get_hopf_field(coord,
             n,
             out,
             xc,
             yc,
             zc);
}

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;
  int myrank;
  MPI_Comm_rank(comm, &myrank);

  parse_command_line_options(argc, argv);

  int   test = strtoul(commandline_option(argc, argv, "-test",     "1", false,
                                          "-test <int> = (1)    : 1) Gaussian profile 2) Zalesak disk"),NULL,10);
  int   merge = strtoul(commandline_option(argc, argv, "-merge",     "1", false,
                                          "-merge <int> = (1)    : 1) no merge 2) complete merge 3) Semi-Merge"),NULL,10);

  {
    tbslas::SimConfig* sim_config       = tbslas::SimConfigSingleton::Instance();
    pvfmm::Profile::Enable(sim_config->profile);
    // =========================================================================
    // PRINT METADATA
    // =========================================================================
    if (!myrank) {
      MetaData_t::Print();
    }

    // =========================================================================
    // TEST CASE
    // =========================================================================
    int max_depth_vel=0;
    int max_depth_con=0;
    pvfmm::BoundaryType bc;
    switch(test) {
      case 1:
        fn_vel = tbslas::get_vorticity_field<double,3>;
        fn_con = get_gaussian_field_cylinder_atT<double,3>;
        bc = pvfmm::FreeSpace;
        break;
      case 2:
        fn_vel = tbslas::get_vorticity_field<double,3>;
        fn_con = get_slotted_cylinder_atT<double,3>;
        bc = pvfmm::FreeSpace;
        break;
      case 3:
        fn_vel = tbslas::get_vel_field_hom_y<double,3>;
        fn_con = get_gaussian_field_cylinder_atT<double,3>;
        bc = pvfmm::Periodic;
        break;
      case 4:                     // regular V, regular C
        max_depth_vel=6;
        max_depth_con=6;
        fn_vel = tbslas::get_vel_field_hom_x<double,3>;
        fn_con = tbslas::get_linear_field_y<double,3>;
        bc = pvfmm::Periodic;
        break;
      case 5:         // regular V, irregular C
        max_depth_vel=6;
        fn_vel = tbslas::get_vorticity_field<double,3>;
        fn_con = gaussian_kernel<double>;
        bc = pvfmm::FreeSpace;
        break;
      case 6:         // irregular V, irregular C
        fn_vel = get_hopf_field<double>;
        fn_con = gaussian_kernel<double>;
        bc = pvfmm::FreeSpace;
        break;
      case 7:  // scaling test case -> uniform fields
        fn_vel = tbslas::get_vel_field_hom_x<double,3>;
        fn_con = tbslas::get_linear_field_y<double,3>;
        bc = pvfmm::Periodic;
        break;
    }

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    sim_config->vtk_filename_prefix     = "advection";
    sim_config->vtk_filename_variable   = "conc";
    sim_config->bc = bc;

    // =========================================================================
    // INIT THE VELOCITY TREE
    // =========================================================================
    Tree_t tvel_curr(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  max_depth_vel?max_depth_vel:sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  fn_vel,
                                  3,
                                  tvel_curr);

    // =========================================================================
    // INIT THE CONCENTRATION TREE
    // =========================================================================
    tcurr = 0;
    Tree_t tconc_curr(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  max_depth_con?max_depth_con:sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  fn_con,
                                  1,
                                  tconc_curr);
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

      snprintf(out_name_buffer,
               sizeof(out_name_buffer),
               sim_config->vtk_filename_format.c_str(),
               tbslas::get_result_dir().c_str(),
               sim_config->vtk_filename_prefix.c_str(),
               sim_config->vtk_filename_variable.c_str(),
               0);
      tconc_curr.Write2File(out_name_buffer, sim_config->vtk_order);
    }

    // =========================================================================
    // RUN
    // =========================================================================
    // set the input_fn to NULL -> needed for adaptive refinement
    std::vector<Node_t*>  ncurr_list = tconc_curr.GetNodeList();
    for(int i = 0; i < ncurr_list.size(); i++) {
      ncurr_list[i]->input_fn = (void (*)(const double* , int , double*))NULL;
    }

    switch(merge) {
      case 2:
        pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
        tbslas::MergeTree(tvel_curr, tconc_curr);
        pvfmm::Profile::Toc();
        break;
      case 3:
        pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
        tbslas::SemiMergeTree(tvel_curr, tconc_curr);
        pvfmm::Profile::Toc();
        break;
    }

    double in_al2,in_rl2,in_ali,in_rli;
    CheckChebOutput<Tree_t>(&tconc_curr,
                            fn_con,
                            1,
                            in_al2,in_rl2,in_ali,in_rli,
                            std::string("Input"));
    typedef tbslas::Reporter<double> Rep;
    if(!myrank) {
      Rep::AddData("TOL", sim_config->tree_tolerance);
      Rep::AddData("ChbOrder", sim_config->tree_chebyshev_order);
      Rep::AddData("MaxDEPTHC", (max_depth_con)?max_depth_con:sim_config->tree_max_depth);
      Rep::AddData("MaxDEPTHV", (max_depth_vel)?max_depth_vel:sim_config->tree_max_depth);

      Rep::AddData("DT", sim_config->dt);
      Rep::AddData("TN", sim_config->total_num_timestep);

      Rep::AddData("InAL2", in_al2);
      Rep::AddData("InRL2", in_rl2);
      Rep::AddData("InALINF", in_ali);
      Rep::AddData("InRLINF", in_rli);
    }


    int num_leaves = tbslas::CountNumLeafNodes(tconc_curr);

    int timestep = 1;
    for (; timestep < sim_config->total_num_timestep+1; timestep++) {
      pvfmm::Profile::Tic("SL", &sim_config->comm, false, 5);
      tbslas::SolveSemilagInSitu(tvel_curr,
                                 tconc_curr,
                                 timestep,
                                 sim_config->dt,
                                 sim_config->num_rk_step);
      pvfmm::Profile::Toc();

      // refine the tree according to the computed values
      pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
      tconc_curr.RefineTree();
      pvfmm::Profile::Toc();

      //Write2File
      if (sim_config->vtk_save) {
        tconc_curr.Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(), sim_config->vtk_order);
      }

      // print error every 100 time steps
      if (timestep % 100 == 0) {
        tcurr = timestep*sim_config->dt;
        double al2,rl2,ali,rli;
        CheckChebOutput<Tree_t>(&tconc_curr,
                                fn_con,
                                1,
                                al2,rl2,ali,rli,
                                std::string("Output_TN" + std::to_string(static_cast<long long>(timestep))));

      }

      // (SEMI) MERGE TO FIX IMBALANCE
      switch(merge) {
        case 2:
          pvfmm::Profile::Tic("CMerge", &sim_config->comm, false, 5);
          tbslas::MergeTree(tvel_curr, tconc_curr);
          pvfmm::Profile::Toc();
          break;
        case 3:
          pvfmm::Profile::Tic("SMerge", &sim_config->comm, false, 5);
          tbslas::SemiMergeTree(tvel_curr, tconc_curr);
          pvfmm::Profile::Toc();
          break;
      }
    }  // end for

    // =========================================================================
    // COMPUTE ERROR
    // =========================================================================
    tcurr = sim_config->total_num_timestep*sim_config->dt;
    double al2,rl2,ali,rli;
    CheckChebOutput<Tree_t>(&tconc_curr,
                            fn_con,
                            1,
                            al2,rl2,ali,rli,
                            std::string("Output"));

    // =========================================================================
    // REPORT RESULTS
    // =========================================================================
    if(!myrank) {
      Rep::AddData("OutAL2", al2);
      Rep::AddData("OutRL2", rl2);
      Rep::AddData("OutALINF", ali);
      Rep::AddData("OutRLINF", rli);
      Rep::Report();
    }
    //Output Profiling results.
    pvfmm::Profile::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
