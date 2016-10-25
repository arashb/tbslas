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

template <class real_t>
void
get_vorticity_field_wrapper(const real_t* points_pos,
                               int num_points,
                               real_t* points_values) {
  real_t omega = 1;
  real_t freq = 0;
  tbslas::get_vorticity_field_tv<real_t,3>(points_pos,
                                           num_points,
                                           tcurr,
                                           points_values,
                                           omega,
                                           freq);

}
template<typename real_t, int sdim>
void
get_gaussian_field_atT(const real_t* points_pos,
                       int num_points,
                       real_t* out) {

  real_t omega = 1.0;

  real_t xc_i = 0.6;
  real_t yc_i = 0.5;
  real_t zc_i = 0.5;

  // compute the time integral
  real_t theta0 = atan((yc_i-0.5)/(xc_i-0.5));
  real_t r      = sqrt((xc_i-0.5)*(xc_i-0.5) + (yc_i-0.5)*(yc_i-0.5));
  real_t theta  = theta0 + omega*tcurr;

  real_t xc = 0.5+r*cos(theta);
  real_t yc = 0.5+r*sin(theta);
  real_t zc = zc_i;

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

template <class real_t>
void
get_vorticity_field_tv_wrapper(const real_t* points_pos,
                               int num_points,
                               real_t* points_values) {
  real_t omega = 1;
  real_t freq = 2*PI;
  tbslas::get_vorticity_field_tv<real_t,3>(points_pos,
                                           num_points,
                                           tcurr,
                                           points_values,
                                           omega,
                                           freq);
}

template <class real_t>
void get_taylor_green_field_tv_wrapper(const real_t* coord,
                            int n,
                            real_t* out) {
  real_t freq = 2*PI;
  tbslas::get_taylor_green_field_tv(coord,
				    n,
				    tcurr,
				    out,
				    freq);

}


template<typename real_t>
void
get_gaussian_field_tv_wrapper(const real_t* points_pos,
                              const int num_points,
                              real_t* out) {
  real_t omega = 1;
  real_t freq = 2*PI;

  // initial positions
  real_t xci      = 0.6;
  real_t yci      = 0.5;
  real_t zci      = 0.5;

  // compute the time integral
  /* real_t theta0 = atan((yci-0.5)/(xci-0.5)); */
  /* real_t r      = sqrt((xci-0.5)*(xci-0.5) + (yci-0.5)*(yci-0.5)); */
  /* real_t theta  = theta0 + omega*(tcurr - tcurr_init - (cos(freq*tcurr) - cos(freq*tcurr_init))/freq); */

  /* real_t xc = 0.5+r*cos(theta); */
  /* real_t yc = 0.5+r*sin(theta); */
  /* real_t zc = zci; */

  real_t xc = xci;
  real_t yc = yci;
  real_t zc = zci;

  const real_t sigma_x = 0.06;
  const real_t sigma_y = 0.06;
  const real_t sigma_z = 0.06;
  const real_t A       = 1.0;

  tbslas::get_gaussian_field_3d<real_t,3>(points_pos,
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

template <class real_t>
void get_diffusion_kernel_tv_wrapper(const real_t* coord,
                                     int n,
                                     real_t* out) {
  const real_t amp = 1e-2;

  real_t omega = 1;
  real_t freq = 2*PI;

  // initial positions
  real_t xci      = 0.6;
  real_t yci      = 0.5;
  real_t zci      = 0.5;

  // compute the time integral
  real_t theta0 = atan((yci-0.5)/(xci-0.5));
  real_t r      = sqrt((xci-0.5)*(xci-0.5) + (yci-0.5)*(yci-0.5));
  real_t theta  = theta0 + omega*(tcurr - tcurr_init - (cos(freq*tcurr) - cos(freq*tcurr_init))/freq);

  real_t xc = 0.5+r*cos(theta);
  real_t yc = 0.5+r*sin(theta);
  real_t zc = zci;

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

template <class real_t>
void get_guassian_kernel_wraper(const real_t* coord,
                     int n,
                     real_t* out) {
  const real_t xc  = 0.7;
  const real_t yc  = 0.7;
  const real_t zc  = 0.7;
  tbslas::gaussian_kernel(coord, n, out, xc, yc, zc);
}

template <class real_t>
void get_multiple_guassian_kernel_wraper(const real_t* coord,
                     int n,
                     real_t* out) {
  // FIRST GAUSSIAN
  const real_t xc1  = 0.7;
  const real_t yc1  = 0.7;
  const real_t zc1  = 0.7;
  std::vector<real_t> out1(n);
  tbslas::gaussian_kernel(coord, n, out1.data(), xc1, yc1, zc1);
  // FIRST GAUSSIAN
  const real_t xc2  = 0.3;
  const real_t yc2  = 0.3;
  const real_t zc2  = 0.3;
  std::vector<real_t> out2(n);
  tbslas::gaussian_kernel(coord, n, out2.data(), xc2, yc2, zc2);
  // FIRST GAUSSIAN
  const real_t xc3  = 0.3;
  const real_t yc3  = 0.3;
  const real_t zc3  = 0.7;
  std::vector<real_t> out3(n);
  tbslas::gaussian_kernel(coord, n, out3.data(), xc3, yc3, zc3);
  for (int i = 0; i < n; i++) {
    out[i] = out1[i] + out2[i] + out3[i];
  }
}

template <class real_t>
void get_taylor_green_field_wrapper(const real_t* coord,
                            int n,
                            real_t* out) {
  tbslas::get_taylor_green_field(coord, n, out);
}


template <class real_t>
void get_exp_alpha_field_wrapper(const real_t* coord,
                                  int n,
                                  real_t* out) {
  const real_t xc = 0.5;
  const real_t yc = 0.5;
  const real_t zc = 0.55;
  const real_t R = 0.3;
  const real_t alpha = EXP_ALPHA;
  tbslas::get_exp_alpha_field(coord, n, out, xc, yc, zc, R, alpha);
}

template <class real_t>
void get_hopf_field_wrapper(const real_t* coord,
                            int n,
                            real_t* out) {
  const real_t xc = 0.5;
  const real_t yc = 0.5;
  const real_t zc = 0.5;
  tbslas::get_hopf_field(coord, n, out, xc, yc, zc);
}

template <class real_t>
void get_diffusion_kernel_atT(const real_t* coord,
                              int n,
                              real_t* out) {
  const real_t amp = 1e-2;
  const real_t xc = 0.5;
  const real_t yc = 0.5;
  const real_t zc = 0.5;
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

template <class real_t>
void get_diffusion_kernel_hopf(const real_t* coord,
                               int n,
                               real_t* out) {
  const real_t amp = 1e-2;
  const real_t xc = 0.6;
  const real_t yc = 0.6;
  const real_t zc = 0.6;
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

template <class real_t>
void get_gaussian_kernel_wrapper(const real_t* coord,
                     int n,
                     real_t* out) {
  const real_t xc  = 0.5;
  const real_t yc  = 0.5;
  const real_t zc  = 0.55;
  tbslas::gaussian_kernel(coord, n, out, xc, yc, zc);
}

template <class real_t>
void get_diffusion_kernel_atT_hom(const real_t* coord,
                                  int n,
                                  real_t* out) {
  const real_t amp = 1e-2;
  const real_t xc = 0.5+(tcurr-tcurr_init)*-0.5;
  const real_t yc = 0.5;
  const real_t zc = 0.5;
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

template <class real_t>
void fn_input_t2(const real_t* coord,
                 int n,
                 real_t* out) {
  tbslas::gaussian_kernel_diffusion_input(coord,
                                          n,
                                          out,
                                          TBSLAS_ALPHA);
}

template <class real_t>
void fn_poten_t2(const real_t* coord,
                 int n,
                 real_t* out) {
  tbslas::gaussian_kernel(coord,
                          n,
                          out);

}

template <class real_t>
void get_taylor_green_field_tv_ns_wrapper(const real_t* points_pos,
				  int num_points,
				  real_t* points_values) {

  real_t freq = 1;
  real_t A = 1;
  real_t B = 1;
  real_t C = 1;
  real_t a = 2*PI;
  real_t b = 2*PI;
  real_t c = 2*PI;
  real_t t = tcurr;
  const real_t* p;
  real_t x,y,z;
  real_t t_factor = cos(freq*t);
  /* const real_t shiftx = 0.3; */
  /* const real_t shifty = 0.3; */
  /* const real_t shiftz = 0.3; */

  for (int i = 0; i < num_points; i++) {
    p = &points_pos[i*COORD_DIM];
    x = -0.5*PI + p[0]*PI;
    y = -0.5*PI + p[1]*PI;
    z = -0.5*PI + p[2]*PI;
    /* x = p[0]-shiftx; */
    /* y = p[1]-shifty; */
    /* z = p[2]-shiftz; */
    points_values[i*3+0] = A*cos(a*x)*sin(b*y)*sin(c*z)*t_factor;
    points_values[i*3+1] = B*sin(a*x)*cos(b*y)*sin(c*z)*t_factor;
    points_values[i*3+2] = C*sin(a*x)*sin(b*y)*cos(c*z)*t_factor;
  }
}


template <class real_t>
void get_two_vortex_tubes_vorticity_field_ns(const real_t* points_pos,
					     int num_points,
					     real_t* points_values) {
  real_t r2 = 0;
  real_t r = 0;
  real_t R = 0;
  real_t r_cutoff = 0.666;

  real_t K = 0.5*exp(2)*log(2);
  real_t A = 0.2;
  
  real_t omega0 = 26.093;
  real_t omega = 0;

  real_t xc, yc, alpha;
  const real_t* p;
  real_t x,y,z;

  real_t scale_lentgh = 2*PI;
  real_t scale_start = scale_lentgh/2.0;
  for (int i = 0; i < num_points; i++) {
    p = &points_pos[i*COORD_DIM];
    
    x = -scale_start + scale_lentgh*p[0];
    y = -scale_start + scale_lentgh*p[1];
    z = -scale_start + scale_lentgh*p[2];

    // x = p[0];
    // y = p[1];
    // z = p[2];

    // FIRST TUBE
    {
      xc = 0.866;
      yc = 0;
      alpha = 2*PI/3.0;
      r2 = (x - (xc + A*cos(alpha)*(1+cos(z))))*(x - (xc + A*cos(alpha)*(1+cos(z))))+
	(y - (yc + A*sin(alpha)*(1+cos(z))))*(y - (yc + A*sin(alpha)*(1+cos(z))));
      R = sqrt(r2)/r_cutoff;
      // omega = 0;
      if (R<1) {
	omega = omega0*(1 - exp(-K*exp(1/(R-1))/R));
	points_values[i*3+0] = -omega*A*cos(alpha)*sin(z);
	points_values[i*3+1] = -omega*A*sin(alpha)*sin(z);
	points_values[i*3+2] =  omega;
      }
    }

    // SECOND TUBE
    {
      xc = -0.866;
      yc = 0;
      alpha = PI/3.0;
      r2 = (x - (xc + A*cos(alpha)*(1+cos(z))))*(x - (xc + A*cos(alpha)*(1+cos(z))))+
    	(y - (yc + A*sin(alpha)*(1+cos(z))))*(y - (yc + A*sin(alpha)*(1+cos(z))));
      R = sqrt(r2)/r_cutoff;
      // omega = 0;
      if (R<1) {
    	omega = omega0*(1 - exp(-K*exp(1/(R-1))/R));
    	points_values[i*3+0] += -omega*A*cos(alpha)*sin(z);
    	points_values[i*3+1] += -omega*A*sin(alpha)*sin(z);
    	points_values[i*3+2] +=  omega;
      }
    }
  }
}
template <class real_t>
void get_two_vortex_tubes_curl_vorticity_field_ns(const real_t* points_pos,
						  int num_points,
						  real_t* points_values) {
  real_t r2 = 0;
  real_t r = 0;
  real_t R = 0;
  real_t r_cutoff = 0.666;

  real_t K = 0.5*exp(2)*log(2);
  real_t A = 0.2;
  
  real_t omega0 = 26.093;
  real_t omega = 0;

  real_t xc, yc, alpha;
  const real_t* p;
  real_t x,y,z;

  real_t scale_lentgh = 2*PI;
  real_t scale_start = scale_lentgh/2.0;
  for (int i = 0; i < num_points; i++) {
    p = &points_pos[i*COORD_DIM];
    
    x = -scale_start + scale_lentgh*p[0];
    y = -scale_start + scale_lentgh*p[1];
    z = -scale_start + scale_lentgh*p[2];

    // x = p[0];
    // y = p[1];
    // z = p[2];

    // FIRST TUBE
    {
      xc = 0.866;
      yc = 0;
      alpha = 2*PI/3.0;

      real_t xcenter = xc + A*cos(alpha)*(1+cos(z));
      real_t ycenter = yc + A*sin(alpha)*(1+cos(z));
      real_t xshifted = x - xcenter;
      real_t yshifted = y - ycenter;
      r2 = xshifted*xshifted + yshifted*yshifted;
      r = sqrt(r2);
      R = r/r_cutoff;

      real_t dRdx = xshifted/(r*r_cutoff);
      real_t dRdy = yshifted/(r*r_cutoff);
      real_t dRdz = (A*cos(alpha)*sin(z)*xshifted + A*sin(alpha)*sin(z)*yshifted)/(r*r_cutoff);

      // real_t dRdx = xshifted/r;
      // real_t dRdy = yshifted/r;
      // real_t dRdz = (A*cos(alpha)*sin(z)*xshifted + A*sin(alpha)*sin(z)*yshifted)/(r);

      real_t dwdR = -omega0*exp(-K*exp(1/(R-1))/R)*exp(1/(R-1))*( K/(R*R) + K/(R*(R-1)*(R-1)) );
    
      real_t dwdx = dwdR*dRdx;
      real_t dwdy = dwdR*dRdy;
      real_t dwdz = dwdR*dRdz;

      omega = 0;
      if (R<1) {
    	omega = omega0*(1 - exp(-K*exp(1/(R-1))/R));
	real_t tmp =  cos(z)*omega + sin(z)*dwdz ;
	points_values[i*3+0] = scale_lentgh*(dwdy + A*sin(alpha)*tmp) ;
	points_values[i*3+1] = scale_lentgh*(-A*cos(alpha)*tmp - dwdx);
	points_values[i*3+2] = scale_lentgh*(-A*sin(alpha)*sin(z)*dwdx + A*cos(alpha)*sin(z)*dwdy);
      }
    }

    // SECOND TUBE
    {
      xc = -0.866;
      yc = 0;
      alpha = PI/3.0;

      real_t xcenter = xc + A*cos(alpha)*(1+cos(z));
      real_t ycenter = yc + A*sin(alpha)*(1+cos(z));
      real_t xshifted = x - xcenter;
      real_t yshifted = y - ycenter;
      r2 = xshifted*xshifted + yshifted*yshifted;
      r = sqrt(r2);
      R = r/r_cutoff;


      real_t dRdx = xshifted/(r*r_cutoff);
      real_t dRdy = yshifted/(r*r_cutoff);
      real_t dRdz = (A*cos(alpha)*sin(z)*xshifted + A*sin(alpha)*sin(z)*yshifted)/(r*r_cutoff);

      // real_t dRdx = xshifted/r;
      // real_t dRdy = yshifted/r;
      // real_t dRdz = (A*cos(alpha)*sin(z)*xshifted + A*sin(alpha)*sin(z)*yshifted)/(r);

      real_t dwdR = -omega0*exp(-K*exp(1/(R-1))/R)*exp(1/(R-1))*( K/(R*R) + K/(R*(R-1)*(R-1)) );
    
      real_t dwdx = dwdR*dRdx;
      real_t dwdy = dwdR*dRdy;
      real_t dwdz = dwdR*dRdz;

      omega = 0;
      if (R<1) {
    	omega = omega0*(1 - exp(-K*exp(1/(R-1))/R));
    	real_t tmp =  cos(z)*omega + sin(z)*dwdz ;
    	points_values[i*3+0] += scale_lentgh*(dwdy + A*sin(alpha)*tmp) ;
    	points_values[i*3+1] += scale_lentgh*(-A*cos(alpha)*tmp - dwdx);
    	points_values[i*3+2] += scale_lentgh*(-A*sin(alpha)*sin(z)*dwdx + A*cos(alpha)*sin(z)*dwdy);
      }
    }
  }
}

template <class real_t>
void get_taylor_green_field_ns(const real_t* points_pos,
			       int num_points,
			       real_t* points_values) {
  real_t A = 2.0/sqrt(3.0);
  real_t B = 2.0/sqrt(3.0);
  real_t C = 2.0/sqrt(3.0);
  real_t a = 2*PI;
  real_t b = 2*PI;
  real_t c = 2*PI;
  real_t theta = 0.0;
  real_t m = +2*PI/3.0;
  real_t n = -2*PI/3.0;
  real_t q = 0;
  const real_t* p;
  real_t x,y,z;
  for (int i = 0; i < num_points; i++) {
    p = &points_pos[i*COORD_DIM];
    x = p[0];
    y = p[1];
    z = p[2];
    points_values[i*3+0] = A*sin(theta+m)*sin(a*x)*cos(b*y)*cos(c*z);
    points_values[i*3+1] = B*sin(theta+n)*cos(a*x)*sin(b*y)*cos(c*z);
    points_values[i*3+2] = C*sin(theta+q)*cos(a*x)*cos(b*y)*sin(c*z);
  }
}

#endif
