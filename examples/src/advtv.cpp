// *************************************************************************
// Copyright (C) 2016 by Arash Bakhtiari
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
#include <vector.hpp>
#include <cheb_utils.hpp>
#include <profile.hpp>

#include <utils.hpp>
#include <field_wrappers.h>

#include <utils/common.h>
#include <utils/metadata.h>
#include <utils/reporter.h>
#include <utils/fields.h>

#include <tree/tree_semilag.h>
#include <tree/tree_utils.h>
#include <tree/tree_set_functor.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

typedef tbslas::MetaData<std::string,
                         std::string,
                         std::string> MetaData_t;

double tcurr_init = 0.0;
double tcurr = 0.0;

void (*fn_vel)(const double* , int , double*)=NULL;
void (*fn_con)(const double* , int , double*)=NULL;

int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;
  int np;
  MPI_Comm_size(comm, &np);
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
    sim_config->num_rk_step = 1;
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
        fn_vel = get_vorticity_field_wrapper<double>;
        fn_con = get_gaussian_field_atT<double,3>;
        // fn_con = get_gaussian_field_cylinder_atT<double,3>;
        // bc = pvfmm::FreeSpace;
        bc = pvfmm::Periodic;
        break;
      case 2:
        fn_vel = tbslas::get_vorticity_field<double,3>;
        fn_con = get_slotted_cylinder_atT<double,3>;
        bc = pvfmm::FreeSpace;
        // bc = pvfmm::Periodic;
        break;
      case 3:
        fn_vel = tbslas::get_vel_field_hom_y<double,3>;
        fn_con = get_gaussian_field_atT<double,3>;
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
        fn_con = get_guassian_kernel_wraper<double>;
        // bc = pvfmm::FreeSpace;
        bc = pvfmm::Periodic;
        break;
      case 6:         // irregular V, irregular C
        fn_vel = get_hopf_field_wrapper<double>;
        fn_con = get_guassian_kernel_wraper<double>;
        // bc = pvfmm::FreeSpace;
        bc = pvfmm::Periodic;
        break;
      case 7:  // scaling test case -> uniform fields
        fn_vel = tbslas::get_vel_field_hom_x<double,3>;
        fn_con = tbslas::get_linear_field_y<double,3>;
        bc = pvfmm::Periodic;
        break;
      case 8:
        fn_vel = get_taylor_green_field_wrapper<double>;
        fn_con = get_guassian_kernel_wraper<double>;
        bc = pvfmm::Periodic;
        break;
      case 9:
        fn_vel = get_taylor_green_field_wrapper<double>;
        fn_con = get_multiple_guassian_kernel_wraper<double>;
        bc = pvfmm::Periodic;
        break;
      case 10:
        fn_vel = get_vorticity_field_tv_wrapper<double>;
        fn_con = get_gaussian_field_tv_wrapper<double>;
        // fn_con = get_gaussian_field_cylinder_atT<double,3>;
        // bc = pvfmm::FreeSpace;
        bc = pvfmm::Periodic;
        break;

    }

    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    sim_config->vtk_filename_variable   = "conc";
    sim_config->bc = bc;

    tbslas::new_nodes<Tree_t::Real_t>(sim_config->tree_chebyshev_order, 3);

    // =========================================================================
    // INIT THE VELOCITY TREES
    // =========================================================================
    std::vector<double>  vel_field_times;
    std::vector<Tree_t*> vel_field_elems;

    tcurr = -sim_config->dt;
    Tree_t* tree;
    for ( int i  = 0; i < 4; i++) {
      tree = new Tree_t(comm);
      tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                    sim_config->tree_num_points_per_octanct,
                                    sim_config->tree_chebyshev_order,
                                    sim_config->tree_max_depth,
                                    sim_config->tree_adap,
                                    sim_config->tree_tolerance,
                                    comm,
                                    fn_vel,
                                    3,
                                    *tree);
      // if (sim_config->vtk_save_rate) {
      //   tree1.Write2File(tbslas::GetVTKFileName(1, "vel").c_str(),
      //                    sim_config->vtk_order);
      // }
      vel_field_times.push_back(tcurr);
      vel_field_elems.push_back(tree);
      tcurr += sim_config->dt;
    }
    tbslas::FieldSetFunctor<double, Tree_t> vel_field_functor(vel_field_elems, vel_field_times);

    if (sim_config->vtk_save_rate) {
      for (int i = 1; i < 4; i++) {
        vel_field_elems[i]->Write2File(tbslas::GetVTKFileName(i-1, "vel").c_str(),
                                       sim_config->vtk_order);
      }
    }

    // =========================================================================
    // INIT THE CONCENTRATION TREE
    // =========================================================================
    tcurr = 0.0;
    Tree_t tcon(comm);
    tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                  sim_config->tree_num_points_per_octanct,
                                  sim_config->tree_chebyshev_order,
                                  max_depth_con?max_depth_con:sim_config->tree_max_depth,
                                  sim_config->tree_adap,
                                  sim_config->tree_tolerance,
                                  comm,
                                  fn_con,
                                  1,
                                  tcon);

    if (sim_config->vtk_save_rate) {
      tcon.Write2File(tbslas::GetVTKFileName(0, sim_config->vtk_filename_variable).c_str(),
                      sim_config->vtk_order);
    }

    // =========================================================================
    // RUN
    // =========================================================================
    // set the input_fn to NULL -> needed for adaptive refinement
    std::vector<Node_t*>  ncurr_list = tcon.GetNodeList();
    for(int i = 0; i < ncurr_list.size(); i++) {
      ncurr_list[i]->input_fn = (void (*)(const double* , int , double*))NULL;
    }

    double in_al2,in_rl2,in_ali,in_rli;
    CheckChebOutput<Tree_t>(&tcon,
                            fn_con,
                            1,
                            in_al2,in_rl2,in_ali,in_rli,
                            std::string("Input"));

    int con_noct_sum = 0;
    int con_noct_max = 0;
    int con_noct_min = 0;

    int vel_noct_sum = 0;
    int vel_noct_max = 0;
    int vel_noct_min = 0;

    // =====================================================================
    // ESTIMATE THE PROBLEM SIZE -> NUMBER OF TREES' OCTANTS
    // =====================================================================
    if (sim_config->profile) {
      int con_noct = tbslas::CountNumLeafNodes(tcon);
      con_noct_sum += con_noct;
      con_noct_max = con_noct;
      con_noct_min = con_noct;

      // int vel_noct = tbslas::CountNumLeafNodes(tvel);
      // vel_noct_sum += vel_noct;
      // vel_noct_max = vel_noct;
      // vel_noct_min = vel_noct;
    }

    Tree_t* tvel_new;

    for (int timestep = 1; timestep < sim_config->total_num_timestep+1; timestep++) {

      // =====================================================================
      // (SEMI) MERGE TO FIX IMBALANCE
      // =====================================================================
      // TODO:
      // pvfmm::Profile::Tic("Merge", &sim_config->comm, false, 5);
      // switch(merge) {
      //   case 2:
      //     tbslas::MergeTree(tvel, tcon);
      //     break;
      //   case 3:
      //     tbslas::SemiMergeTree(tvel, tcon);
      //     break;
      // }
      // pvfmm::Profile::Toc();

      // =====================================================================
      // ESTIMATE THE PROBLEM SIZE -> NUMBER OF TREES' OCTANTS
      // =====================================================================
      pvfmm::Profile::Tic("CountLVS", &sim_config->comm, false, 5);
      if (sim_config->profile) {
        int con_noct = tbslas::CountNumLeafNodes(tcon);
        con_noct_sum += con_noct;
        if (con_noct > con_noct_max) con_noct_max = con_noct;
        if (con_noct < con_noct_min) con_noct_min = con_noct;

        // int vel_noct = tbslas::CountNumLeafNodes(tvel);
        // vel_noct_sum += vel_noct;
        // if (vel_noct > vel_noct_max) vel_noct_max = vel_noct;
        // if (vel_noct < vel_noct_min) vel_noct_min = vel_noct;
      }
      pvfmm::Profile::Toc();

      pvfmm::Profile::Tic(std::string("Solve_TN" + tbslas::ToString(static_cast<long long>(timestep))).c_str(), &comm, true);
      {
        // ===================================================================
        // SOLVE SEMILAG
        // ===================================================================
        pvfmm::Profile::Tic("SLM", &sim_config->comm, false, 5);
        tbslas::SolveSemilagInSitu(vel_field_functor,
                                   tcon,
                                   timestep,
                                   sim_config->dt,
                                   sim_config->num_rk_step);
        pvfmm::Profile::Toc();
      }
      pvfmm::Profile::Toc();        // solve

      // ===================================================================
      // REFINE TREE
      // ===================================================================
      pvfmm::Profile::Tic("RefineTree", &sim_config->comm, false, 5);
      tcon.RefineTree();
      pvfmm::Profile::Toc();

      pvfmm::Profile::Tic("Balance21", &sim_config->comm, false, 5);
      tcon.Balance21(sim_config->bc);
      pvfmm::Profile::Toc();

      // =====================================================================
      // Update velocity functor
      // =====================================================================
      pvfmm::Profile::Tic("UpdateVel", &sim_config->comm, false, 5);
      tcurr = (timestep+2)*sim_config->dt;
      tvel_new  = new Tree_t(comm);
      tbslas::ConstructTree<Tree_t>(sim_config->tree_num_point_sources,
                                    sim_config->tree_num_points_per_octanct,
                                    sim_config->tree_chebyshev_order,
                                    sim_config->tree_max_depth,
                                    sim_config->tree_adap,
                                    sim_config->tree_tolerance,
                                    comm,
                                    fn_vel,
                                    3,
                                    *tvel_new);
      if (sim_config->vtk_save_rate && ((timestep+2) % sim_config->vtk_save_rate == 0)) {
        tvel_new->Write2File(tbslas::GetVTKFileName((timestep+2), "vel").c_str(),
                         sim_config->vtk_order);
      }
      vel_field_functor.update(tvel_new, tcurr);
      pvfmm::Profile::Toc();

      // =====================================================================
      // Write2File and print error every N time steps
      // =====================================================================
      if (sim_config->vtk_save_rate) {
        if ( timestep % sim_config->vtk_save_rate == 0) {
          tcon.Write2File(tbslas::GetVTKFileName(timestep, sim_config->vtk_filename_variable).c_str(),
                          sim_config->vtk_order);
          tcurr = timestep*sim_config->dt;
          double al2,rl2,ali,rli;
          CheckChebOutput<Tree_t>(&tcon,
                                  fn_con,
                                  1,
                                  al2,rl2,ali,rli,
                                  std::string("Output_TN" + tbslas::ToString(static_cast<long long>(timestep))));
        }
      }

    }  // end for

    // =========================================================================
    // COMPUTE ERROR
    // =========================================================================
    tcurr = sim_config->total_num_timestep*sim_config->dt;
    double al2,rl2,ali,rli;
    CheckChebOutput<Tree_t>(&tcon,
                            fn_con,
                            1,
                            al2,rl2,ali,rli,
                            std::string("Output"));

    // =========================================================================
    // REPORT RESULTS
    // =========================================================================
    int tcon_max_depth = 0;
    int tvel_max_depth = 0;
    tbslas::GetTreeMaxDepth(tcon, tcon_max_depth);
    // tbslas::GetTreeMaxDepth(tvel, tvel_max_depth);

    typedef tbslas::Reporter<double> Rep;
    if(!myrank) {
      Rep::AddData("NP", np, tbslas::REP_INT);
      Rep::AddData("OMP", sim_config->num_omp_threads, tbslas::REP_INT);

      Rep::AddData("TOL", sim_config->tree_tolerance);
      Rep::AddData("Q", sim_config->tree_chebyshev_order, tbslas::REP_INT);

      Rep::AddData("MaxD", sim_config->tree_max_depth, tbslas::REP_INT);
      Rep::AddData("CMaxD", tcon_max_depth, tbslas::REP_INT);
      Rep::AddData("VMaxD", tvel_max_depth, tbslas::REP_INT);

      // Rep::AddData("CBC", sim_config->use_cubic?1:0, tbslas::REP_INT);
      // Rep::AddData("CUF", sim_config->cubic_upsampling_factor, tbslas::REP_INT);

      Rep::AddData("DT", sim_config->dt);
      Rep::AddData("TN", sim_config->total_num_timestep, tbslas::REP_INT);
      Rep::AddData("TEST", test, tbslas::REP_INT);
      Rep::AddData("MERGE", merge, tbslas::REP_INT);

      Rep::AddData("InAL2", in_al2);
      Rep::AddData("OutAL2", al2);

      Rep::AddData("InRL2", in_rl2);
      Rep::AddData("OutRL2", rl2);

      Rep::AddData("InALINF", in_ali);
      Rep::AddData("OutALINF", ali);

      Rep::AddData("InRLINF", in_rli);
      Rep::AddData("OutRLINF", rli);

      Rep::AddData("CMinNOCT", con_noct_min, tbslas::REP_INT);
      Rep::AddData("CAvgNOCT", con_noct_sum/(sim_config->total_num_timestep+1),
                   tbslas::REP_INT); // NUMBER OF TIMESTEPS + INITIAL TREE
      Rep::AddData("CMaxNOCT", con_noct_max, tbslas::REP_INT);

      Rep::AddData("VMinNOCT", vel_noct_min, tbslas::REP_INT);
      Rep::AddData("VAvgNOCT", vel_noct_sum/(sim_config->total_num_timestep+1),
                   tbslas::REP_INT); // NUMBER OF TIMESTEPS + INITIAL TREE
      Rep::AddData("VMaxNOCT", vel_noct_max, tbslas::REP_INT);

      Rep::Report();
    }

    // DEALLOCATE MEMORY
    // TODO: deallocate the

    //Output Profiling results.
    pvfmm::Profile::print(&comm);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
