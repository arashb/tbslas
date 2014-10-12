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

#include <semilag/semilag.h>
#include <semilag/utils.h>

typedef pvfmm::Cheb_Node<double> Node_t;
typedef pvfmm::MPI_Tree<Node_t> Tree_t;

const int DATA_DOF=1;

template<typename real_t>
std::vector<int>
isOutside(Node_t* n,
          const std::vector<real_t> x,
          const std::vector<real_t> y,
          const std::vector<real_t> z) {
  assert((x.size() == y.size()) && (y.size() == z.size()));

  real_t* node_coord = n->Coord();
  int depth          = n->Depth();
  real_t length      = static_cast<real_t>(std::pow(0.5, depth));

  real_t xmin = node_coord[0];
  real_t xmax = xmin + length;
  real_t ymin = node_coord[1];
  real_t ymax = ymin + length;
  real_t zmin = node_coord[2];
  real_t zmax = zmin + length;

  std::vector<int> out_index_list;
  for (int i = 0; i < x.size(); i++) {
    if ( x[i] < xmin || x[i] > xmax) {
      out_index_list.push_back(i);
      continue;
    }
    if ( y[i] < ymin || y[i] > ymax) {
      out_index_list.push_back(i);
      continue;
    }
    if ( z[i] < zmin || z[i] > zmax) {
      out_index_list.push_back(i);
      continue;
    }
  }
  return out_index_list;
}

namespace tbslas {

template<typename real_t>
class FieldFunctor {
 public:
  explicit FieldFunctor(Node_t* node): node_(node) {
  }

  virtual ~FieldFunctor() {
  }

  void operator () (const real_t* points_pos,
                    int num_points,
                    real_t* out) {
    //slas::get_gaussian_field<double,3>(points_pos, num_points, out);
    // TODO: remove this part
    // std::vector<real_t> x_c,y_c,z_c;
    // for (int i = 0; i < num_points; i++) {
    //   x_c.push_back(points_pos[i*COORD_DIM+0]);
    //   y_c.push_back(points_pos[i*COORD_DIM+1]);
    //   z_c.push_back(points_pos[i*COORD_DIM+2]);
    // }
    // std::vector<int> outsiders = isOutside(node_, x_c, y_c, z_c);
    // assert(outsiders.empty());

    // TODO: optimize! you do not need to get the values of each point separately
    for (int i = 0; i < num_points; i++) {
      std::vector<real_t> x,y,z;
      x.push_back(points_pos[i*COORD_DIM+0]);
      y.push_back(points_pos[i*COORD_DIM+1]);
      z.push_back(points_pos[i*COORD_DIM+2]);
      node_->ReadVal(x, y, z, &out[i*DATA_DOF]);
    }
  }

 private:
  Node_t* node_;
};

template <class real_t>
void semilag_construct_tree(const size_t N,
                            const size_t M,
                            const int cheb_deg,
                            const int depth,
                            const bool adap,
                            const real_t tol,
                            const MPI_Comm& comm,
                            Tree_t& tree) {

  //Various parameters.
  typename Node_t::NodeData tree_data;
  tree_data.dim=COORD_DIM;
  tree_data.max_depth=depth;
  tree_data.cheb_deg=cheb_deg;

  //Set input function pointer
  //tree_data.input_fn=fn<real_t>;
  tree_data.input_fn=slas::get_gaussian_field<real_t,3>;
  tree_data.data_dof=DATA_DOF;
  tree_data.tol=tol;

  //Set source coordinates.
  std::vector<real_t> pt_coord;
  pt_coord=point_distrib<real_t>(UnifGrid,N,comm);
  tree_data.max_pts=M; // Points per octant.
  tree_data.pt_coord=pt_coord;

  //initialize with input data.
  tree.Initialize(&tree_data);
  tree.RefineTree();
  tree.Balance21(pvfmm::FreeSpace);
  //tree.RedistNodes();
}

template <typename real_t>
void semilag_advect_tree(Tree_t& tree, const int cheb_deg, Tree_t& tree_next) {
  int sdim        = 3;
  int timestep    = 1;
  real_t dt       = 0.5;
  int num_rk_step = 1;

  Node_t* n = tree.PostorderFirst();
  Node_t* n_next = tree_next.PostorderFirst();
  while (n != NULL) {
    if (n->IsLeaf() && !n->IsGhost()) {
      // compute chebychev points positions on the fly
      std::vector<real_t> points_pos = pvfmm::cheb_nodes<real_t>(cheb_deg, sdim);

      int depth          = n->Depth();
      real_t length      = static_cast<real_t>(std::pow(0.5, depth));
      real_t* node_coord = n->Coord();

      printf("NODE: [%f, %f, %f]\n",
             node_coord[0],
             node_coord[1],
             node_coord[2]);
      // TODO: figure aut a way to optimize this part
      // scale the cheb points
      int num_points = points_pos.size()/COORD_DIM;
      for (int i = 0; i < num_points; i++) {
        points_pos[i*COORD_DIM+0] = node_coord[0] + length * points_pos[i*COORD_DIM+0];
        points_pos[i*COORD_DIM+1] = node_coord[1] + length * points_pos[i*COORD_DIM+1];
        points_pos[i*COORD_DIM+2] = node_coord[2] + length * points_pos[i*COORD_DIM+2];
      }

      std::vector<real_t> points_val(num_points);
      FieldFunctor<real_t> con_field(tree.RootNode());

      //field(points_pos.data(), num_points, points_val.data());
      slas::semilag_rk2(slas::get_vorticity_field<double,3>,
                        con_field,
                        points_pos,
                        sdim,
                        timestep,
                        dt,
                        num_rk_step,
                        points_val
                        );

      pvfmm::cheb_approx<real_t, real_t>(points_val.data(),
                                         cheb_deg,
                                         DATA_DOF,
                                         &(n_next->ChebData()[0])
                                         );
    }
    n = tree.PostorderNxt(n);
    n_next = tree_next.PostorderNxt(n_next);
  }
}

}  // namespace tbslas


int main (int argc, char **argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm=MPI_COMM_WORLD;

  // Read command line options.
  commandline_option_start(argc, argv);
  omp_set_num_threads( atoi(commandline_option(argc, argv,  "-omp",     "1", false, "-omp  <int> = (1)    : Number of OpenMP threads."          )));
  size_t   N=(size_t)strtod(commandline_option(argc, argv,    "-N",     "1",  true, "-N    <int>          : Number of point sources."           ),NULL);
  size_t   M=(size_t)strtod(commandline_option(argc, argv,    "-M",     "1", false, "-M    <int>          : Number of points per octant."       ),NULL);
  int      q=       strtoul(commandline_option(argc, argv,    "-q",    "14", false, "-q    <int> = (14)   : Chebyshev order (+ve integer)."     ),NULL,10);
  int      d=       strtoul(commandline_option(argc, argv,    "-d",    "15", false, "-d    <int> = (15)   : Maximum tree depth."                ),NULL,10);
  double tol=        strtod(commandline_option(argc, argv,  "-tol",  "1e-5", false, "-tol <real> = (1e-5) : Tolerance for adaptive refinement." ),NULL);
  bool  adap=              (commandline_option(argc, argv, "-adap",    NULL, false, "-adap                : Adaptive tree refinement."          )!=NULL);
  commandline_option_end(argc, argv);

  {
    Tree_t tree_current(comm);
    tbslas::semilag_construct_tree<double>(N, M, q, d, adap, tol, comm, tree_current);
    Tree_t tree_next(comm);
    tbslas::semilag_construct_tree<double>(N, M, q, d, adap, tol, comm, tree_next);
    //Find error in Chebyshev approximation.
    // CheckChebOutput<Tree_t>(&tree,
    //                         (typename TestFn<double>::Fn_t) &slas::get_gaussian_field<double,3>,
    //                         DATA_DOF,
    //                         "Input");
    tree_current.Write2File("result/output_00_",4);
    tbslas::semilag_advect_tree<double>(tree_current, q, tree_next);
    tree_next.Write2File("result/output_01_",4);
  }

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
