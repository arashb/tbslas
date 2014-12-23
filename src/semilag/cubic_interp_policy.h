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
