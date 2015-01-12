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

#ifndef SRC_SEMILAG_CUBIC_INTERP_POLICY_H_
#define SRC_SEMILAG_CUBIC_INTERP_POLICY_H_

#include <vector>
#include <utils/profile.h>

namespace tbslas {

template<typename real_t>
class CubicInterpPolicy {
 public:
  CubicInterpPolicy();
  virtual ~CubicInterpPolicy();

 public:
  void
  interp(const std::vector<real_t>& grid_points,
         const int sdim,
         const std::vector<real_t>& grid_values,
         const int vdim,
         const std::vector<real_t>& query_points,
         std::vector<real_t>& query_values) const;

  void
  interp(const std::vector<real_t>& grid_points,
         const int sdim,
         const std::vector< std::vector<real_t> >& grid_values,
         const int vdim,
         const std::vector<real_t>& time,
         const std::vector<real_t>& query_points,
         const real_t query_time,
         std::vector<real_t>& query_values) const;

  static void
  interp(const std::vector<real_t>& reg_grid_vals,
         const int data_dof,
         const real_t spacing,
         const std::vector<real_t>& query_points,
         std::vector<real_t>& query_values) {
    // number of query points
    int n_pts = query_points.size()/3;
    // regular grid parameters
    int reg_grid_resolution = 1.0/spacing + 1;
    int reg_grid_num_points = std::pow(reg_grid_resolution, 3);
    query_values.resize(n_pts*data_dof);

    assert(reg_grid_num_points*data_dof == reg_grid_vals.size());
    const real_t* coord_ptr = query_points.data();

    for(int data_dim_cnt = 0; data_dim_cnt < data_dof; data_dim_cnt++) {
      const real_t* reg_grid_vals_dim =
          &reg_grid_vals[data_dim_cnt*reg_grid_num_points];
      for ( int pi = 0; pi < n_pts; pi++) {

        //in case that target coordinate is outside the unit domain->set to 0
        if ( coord_ptr[pi*COORD_DIM+0] < 0 || coord_ptr[pi*COORD_DIM+0] > 1.0 ||
             coord_ptr[pi*COORD_DIM+1] < 0 || coord_ptr[pi*COORD_DIM+1] > 1.0 ||
             coord_ptr[pi*COORD_DIM+2] < 0 || coord_ptr[pi*COORD_DIM+2] > 1.0
             ) {
          query_values[pi*data_dof + data_dim_cnt] = 0;
          continue;
        }

        real_t xq = coord_ptr[pi*COORD_DIM+0];
        real_t yq = coord_ptr[pi*COORD_DIM+1];
        real_t zq = coord_ptr[pi*COORD_DIM+2];

        // scale to [0,1] in local node
        // Real_t xq =(coord_ptr[pi*COORD_DIM+0]-c[0])*s;
        // Real_t yq =(coord_ptr[pi*COORD_DIM+1]-c[1])*s;
        // Real_t zq =(coord_ptr[pi*COORD_DIM+2]-c[2])*s;

        real_t dx(xq/spacing);
        real_t dy(yq/spacing);
        real_t dz(zq/spacing);

        // Calculate the corresponding lower-bound grid indices.
        int xi = static_cast<int>(std::floor(dx));
        int yi = static_cast<int>(std::floor(dy));
        int zi = static_cast<int>(std::floor(dz));
        assert(xi >= 0 && yi >= 0  && zi >= 0);
        assert(xi < reg_grid_resolution &&
               yi < reg_grid_resolution &&
               zi < reg_grid_resolution);

        int xshift, yshift, zshift;
        if (xi == 0) xshift = 0;
        else if (xi == reg_grid_resolution -2)
          xshift = xi -2;
        else xshift = xi - 1;

        if (yi == 0) yshift = 0;
        else if (yi == reg_grid_resolution -2)
          yshift = yi -2;
        else yshift = yi - 1;

        if (zi == 0) zshift = 0;
        else if (zi == reg_grid_resolution -2)
          zshift = zi -2;
        else zshift = zi - 1;

        real_t xx[4];
        real_t yy[4];
        real_t zz[4];
        for (int i = 0; i < 4; i++) {
          // xx[i] = reg_grid_coord_1d[xshift + i];
          // yy[i] = reg_grid_coord_1d[yshift + i];
          // zz[i] = reg_grid_coord_1d[zshift + i];
          xx[i] = (xshift + i)*spacing;
          yy[i] = (yshift + i)*spacing;
          zz[i] = (zshift + i)*spacing;

        }

        int xli, yli, zli;
        int index;
        real_t pp[4][4][4];
        for(int i = 0; i < 4; i++) {
          xli = xshift + i;
          for ( int j = 0; j < 4; j++) {
            yli = yshift + j;
            for ( int k = 0; k < 4; k++) {
              zli = zshift + k;
              index = xli + yli*reg_grid_resolution + zli*reg_grid_resolution*reg_grid_resolution;
              pp[i][j][k] = reg_grid_vals_dim[index];
            }
          }
        }
        query_values[pi*data_dof + data_dim_cnt] =
            tbslas::CubicInterpPolicy<real_t>::InterpCubic3D(xq, yq, zq,
                                                             xx, yy, zz,
                                                             pp);
      }  // for target point
    }  // for data_dof
  }
 private:

  inline static real_t
  cubic_interp_1d(real_t p[4], real_t x);

  inline static real_t
  cubic_interp_2d(real_t p[4][4], real_t x, real_t y);

  inline static real_t
  h00(real_t t) {
    return 2*t*t*t - 3*t*t + 1;
  }

  inline static real_t
  h10(real_t t) {
    return t*t*t - 2*t*t + t;
  }

  inline static real_t
  h01(real_t t) {
    return -2*t*t*t + 3*t*t;
  }

  inline static real_t
  h11(real_t t) {
    return t*t*t - t*t;
  }

  inline static real_t
  tangent(real_t tk_1,
          real_t tk,
          real_t tk1,
          real_t pk_1,
          real_t pk,
          real_t pk1) {
    return (pk1 - pk)*0.5/(tk1 - tk) + (pk - pk_1)*0.5/(tk - tk_1) ;
  }

 public:

  static real_t
  InterpCubic1D(real_t x,
                real_t xx[4],
                real_t pp[4]
                ) {
    real_t mk = CubicInterpPolicy::tangent(xx[0], xx[1], xx[2],
                                           pp[0], pp[1], pp[2]);

    real_t mk1 = CubicInterpPolicy::tangent(xx[1], xx[2], xx[3],
                                            pp[1], pp[2], pp[3]);
    real_t t = (x - xx[1])/(xx[2] - xx[1]);
    // printf("T: %f, x: %f, xx[1]: %f xx[2]: %f\n", t, x, xx[1], xx[2] );

    real_t val =
        CubicInterpPolicy::h00(t)*pp[1] +
        CubicInterpPolicy::h10(t)*(xx[2] - xx[1])*mk +
        CubicInterpPolicy::h01(t)*pp[2] +
        CubicInterpPolicy::h11(t)*(xx[2] - xx[1])*mk1;
    // if (val < -700)
      // printf("HT: %f, mk: %f, mk1: %f\n", t, mk, mk1 );
    return val;
  }

  static real_t
  InterpCubic2D(real_t x,
                real_t y,
                real_t xx[4],
                real_t yy[4],
                real_t pp[4][4]
                ) {
    real_t val[4];
    val[0] = CubicInterpPolicy::InterpCubic1D(y, yy, pp[0]);
    val[1] = CubicInterpPolicy::InterpCubic1D(y, yy, pp[1]);
    val[2] = CubicInterpPolicy::InterpCubic1D(y, yy, pp[2]);
    val[3] = CubicInterpPolicy::InterpCubic1D(y, yy, pp[3]);
    return CubicInterpPolicy::InterpCubic1D(x, xx, val);
  }

  static real_t
  InterpCubic3D(real_t x,
                real_t y,
                real_t z,
                real_t xx[4],
                real_t yy[4],
                real_t zz[4],
                real_t pp[4][4][4]
                ) {
    real_t val[4];
    val[0] = CubicInterpPolicy::InterpCubic2D(y, z, yy, zz, pp[0]);
    val[1] = CubicInterpPolicy::InterpCubic2D(y, z, yy, zz, pp[1]);
    val[2] = CubicInterpPolicy::InterpCubic2D(y, z, yy, zz, pp[2]);
    val[3] = CubicInterpPolicy::InterpCubic2D(y, z, yy, zz, pp[3]);
    return CubicInterpPolicy::InterpCubic1D(x, xx, val);
  }

};  // class CubicInterpPolicy

}  // namespace tbslas

#include "cubic_interp_policy.inc"

#endif  // SRC_SEMILAG_CUBIC_INTERP_POLICY_H_
