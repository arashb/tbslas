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
#ifndef __FIELD_WRAPPERS_H__
#define __FIELD_WRAPPERS_H__

#include <utils/fields.h>


extern double TBSLAS_DIFF_COEFF;
extern double TBSLAS_ALPHA;
extern double EXP_ALPHA;
extern double tcurr_init;
extern double tcurr;

template<typename real_t, int sdim>
void
get_gaussian_field_atT(const real_t* points_pos,
                       int num_points,
                       real_t* out) {
  real_t xc      = 0.6;
  real_t yc      = 0.5;
  real_t zc      = 0.5;

  real_t theta0 = atan((yc-0.5)/(xc-0.5));
  real_t r = sqrt((xc-0.5)*(xc-0.5) + (yc-0.5)*(yc-0.5));
  xc = 0.5+r*cos(theta0 + tcurr);
  yc = 0.5+r*sin(theta0 + tcurr);
  const real_t theta   = 0.0;
  const real_t sigma_x = 0.06;
  const real_t sigma_y = 0.06;
  const real_t sigma_z = 0.06;
  const real_t A       = 1.0;

  tbslas::get_gaussian_field_3d<real_t, sdim>(points_pos,
                                              num_points,
                                              out,
                                              xc,
                                              yc,
                                              zc,
                                              A,
                                              sigma_x,
                                              sigma_y,
                                              sigma_z );


}

template<typename real_t, int sdim>
void
get_gaussian_field_cylinder_atT(const real_t* points_pos,
                       int num_points,
                       real_t* out) {
  real_t xc      = 0.6;
  real_t yc      = 0.5;
  real_t zc      = 0.5;

  real_t theta0 = atan((yc-0.5)/(xc-0.5));
  real_t r = sqrt((xc-0.5)*(xc-0.5) + (yc-0.5)*(yc-0.5));
  xc = 0.5+r*cos(theta0 + tcurr);
  yc = 0.5+r*sin(theta0 + tcurr);
  const real_t theta   = 0.0;
  const real_t sigma_x = 0.06;
  const real_t sigma_y = 0.06;
  const real_t sigma_z = 0.06;
  const real_t A       = 1.0;

  tbslas::get_gaussian_field_cylinder<real_t, sdim>(points_pos,
                                                    num_points,
                                                    out,
                                                    xc,
                                                    yc,
                                                    theta,
                                                    sigma_x,
                                                    sigma_y,
                                                    A);
}


template<typename real_t, int sdim>
void
get_slotted_cylinder_atT(const real_t* points_pos,
                         int num_points,
                         real_t* out) {
  real_t xc = 0.5;
  real_t yc = 0.5;
  real_t zc = 0.5;
  real_t R  = 0.3;
  real_t w  = 0.1;
  real_t a  = tcurr;
  tbslas::get_slotted_cylinder<real_t, sdim>(points_pos,
                                             num_points,
                                             out,
                                             xc,
                                             yc,
                                             zc,
                                             R,
                                             w,
                                             a);
}

template <class Real_t>
void get_guassian_kernel_wraper(const Real_t* coord,
                     int n,
                     Real_t* out) {
  const Real_t xc  = 0.7;
  const Real_t yc  = 0.7;
  const Real_t zc  = 0.7;
  tbslas::gaussian_kernel(coord, n, out, xc, yc, zc);
}

template <class Real_t>
void get_multiple_guassian_kernel_wraper(const Real_t* coord,
                     int n,
                     Real_t* out) {
  // FIRST GAUSSIAN
  const Real_t xc1  = 0.7;
  const Real_t yc1  = 0.7;
  const Real_t zc1  = 0.7;
  std::vector<Real_t> out1(n);
  tbslas::gaussian_kernel(coord, n, out1.data(), xc1, yc1, zc1);
  // FIRST GAUSSIAN
  const Real_t xc2  = 0.3;
  const Real_t yc2  = 0.3;
  const Real_t zc2  = 0.3;
  std::vector<Real_t> out2(n);
  tbslas::gaussian_kernel(coord, n, out2.data(), xc2, yc2, zc2);
  // FIRST GAUSSIAN
  const Real_t xc3  = 0.3;
  const Real_t yc3  = 0.3;
  const Real_t zc3  = 0.7;
  std::vector<Real_t> out3(n);
  tbslas::gaussian_kernel(coord, n, out3.data(), xc3, yc3, zc3);
  for (int i = 0; i < n; i++) {
    out[i] = out1[i] + out2[i] + out3[i];
  }
}

template <class Real_t>
void get_taylor_green_field_wrapper(const Real_t* coord,
                            int n,
                            Real_t* out) {
  tbslas::get_taylor_green_field(coord, n, out);
}

template <class Real_t>
void get_exp_alpha_field_wrapper(const Real_t* coord,
                                  int n,
                                  Real_t* out) {
  const Real_t xc = 0.5;
  const Real_t yc = 0.5;
  const Real_t zc = 0.55;
  const Real_t R = 0.3;
  const Real_t alpha = EXP_ALPHA;
  tbslas::get_exp_alpha_field(coord, n, out, xc, yc, zc, R, alpha);
}

template <class Real_t>
void get_hopf_field_wrapper(const Real_t* coord,
                            int n,
                            Real_t* out) {
  const Real_t xc = 0.5;
  const Real_t yc = 0.5;
  const Real_t zc = 0.5;
  tbslas::get_hopf_field(coord, n, out, xc, yc, zc);
}

template <class Real_t>
void get_diffusion_kernel_atT(const Real_t* coord,
                              int n,
                              Real_t* out) {
  const Real_t amp = 1e-2;
  const Real_t xc = 0.5;
  const Real_t yc = 0.5;
  const Real_t zc = 0.5;
  tbslas::diffusion_kernel(coord,
                           n,
                           out,
                           TBSLAS_DIFF_COEFF,
                           tcurr,
                           amp,
                           xc,
                           yc,
                           zc);
}

template <class Real_t>
void get_diffusion_kernel_hopf(const Real_t* coord,
                               int n,
                               Real_t* out) {
  const Real_t amp = 1e-2;
  const Real_t xc = 0.6;
  const Real_t yc = 0.6;
  const Real_t zc = 0.6;
  double time_curr = 2.5;
  tbslas::diffusion_kernel(coord,
                           n,
                           out,
                           TBSLAS_DIFF_COEFF,
                           time_curr,
                           amp,
                           xc,
                           yc,
                           zc);
}

template <class Real_t>
void get_gaussian_kernel_wrapper(const Real_t* coord,
                     int n,
                     Real_t* out) {
  const Real_t xc  = 0.5;
  const Real_t yc  = 0.5;
  const Real_t zc  = 0.55;
  tbslas::gaussian_kernel(coord, n, out, xc, yc, zc);
}

template <class Real_t>
void get_diffusion_kernel_atT_hom(const Real_t* coord,
                                  int n,
                                  Real_t* out) {
  const Real_t amp = 1e-2;
  const Real_t xc = 0.5+(tcurr-tcurr_init)*-0.5;
  const Real_t yc = 0.5;
  const Real_t zc = 0.5;
  tbslas::diffusion_kernel(coord,
                           n,
                           out,
                           TBSLAS_DIFF_COEFF,
                           tcurr,
                           amp,
                           xc,
                           yc,
                           zc);
}

template <class Real_t>
void fn_input_t2(const Real_t* coord,
                 int n,
                 Real_t* out) {
  tbslas::gaussian_kernel_diffusion_input(coord,
                                          n,
                                          out,
                                          TBSLAS_ALPHA);
}

template <class Real_t>
void fn_poten_t2(const Real_t* coord,
                 int n,
                 Real_t* out) {
  tbslas::gaussian_kernel(coord,
                          n,
                          out);
}

template <class Real_t>
void
get_vorticity_field_tv_wrapper(const Real_t* points_pos,
                               int num_points,
                               Real_t* points_values) {
  tbslas::get_vorticity_field_tv<Real_t,3>(points_pos,
                                           num_points,
                                           tcurr,
                                           points_values);

}

#endif
