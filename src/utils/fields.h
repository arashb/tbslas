// *************************************************************************
// Copyright (C) 2015 by Arash Bakhtiari
// You may not use this file except in compliance with the License.
// You obtain a copy of the License in the LICENSE file.

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// *************************************************************************

#ifndef SRC_UTILS_FIELDS_H_
#define SRC_UTILS_FIELDS_H_

#define PI 3.14159265

namespace tbslas {

//////////////////////////////////////////////////////////////////////
// VARIOUS FIELD GENERATOR FUNCTIONS
//////////////////////////////////////////////////////////////////////
///////////////////////////////////
// TAYLOR-GREEN VORTEX
///////////////////////////////////
// u = A*cos(a*x)sin(b*y)sin(c*z)
// v = B*sin(a*x)cos(b*y)sin(c*z)
// w = C*sin(a*x)sin(b*y)cos(c*z)
///////////////////////////////////
template<typename real_t>
void
get_taylor_green_field(const real_t* points_pos,
                       int num_points,
                       real_t* points_values,
                       real_t A = 1,
                       real_t B = 1,
                       real_t C = 1,
                       real_t a = 2*PI,
                       real_t b = 2*PI,
                       real_t c = 2*PI) {
  const real_t* p;
  for (int i = 0; i < num_points; i++) {
    p = &points_pos[i*COORD_DIM];
    points_values[i*3+0] = A*cos(a*p[0])*sin(b*p[1])*sin(c*p[2]);
    points_values[i*3+1] = B*sin(a*p[0])*cos(b*p[1])*sin(c*p[2]);
    points_values[i*3+2] = C*sin(a*p[0])*sin(b*p[1])*cos(c*p[2]);
  }
}

///////////////////////////////////
// VORTICITY FIELD
///////////////////////////////////
///////////////////////////////////
template<typename real_t, int sdim>
void
get_vorticity_field(const real_t* points_pos,
                    int num_points,
                    real_t* points_values) {
  real_t omega       = 1;
  real_t time_factor = 1;//+sin(2*3.14159*time);
  for (int i = 0; i < num_points; i++) {
    points_values[i*3+0] =
        omega*(0.5 - points_pos[i*sdim+1])*time_factor;
    points_values[i*3+1] =
        omega*(points_pos[i*sdim+0] - 0.5)*time_factor;
    points_values[i*3+2] = 0;
  }
}

template<typename real_t, int sdim>
void
get_linear_field_y(const real_t* points_pos,
                 int num_points,
                 real_t* points_values) {
  real_t omega = 1;
  real_t time_factor = 1;//+sin(2*3.14159*time);
  for (int i = 0; i < num_points; i++) {
    points_values[i] =
        omega*(0.5 - points_pos[i*sdim+1])*time_factor;
    // points_values[i*3+1] =
    //     omega*(points_pos[i*sdim+0] - 0.5)*time_factor;
    // points_values[i*3+2] = 0;
  }
}

template<typename real_t>
void
get_hopf_field(const real_t* points_pos,
           int num_points,
           real_t* points_values,
           const real_t xc = 0.5,
           const real_t yc = 0.5,
           const real_t zc = 0.5,
           const real_t A = 0.05,
           const real_t r = 0.1) {
  real_t x,y,z;
  for (int i = 0; i < num_points; i++) {
    const real_t* c = &points_pos[i*COORD_DIM];
    x = c[0] - xc;
    y = c[1] - yc;
    z = c[2] - zc;
    real_t r_2 = x*x+y*y+z*z;
    real_t coeff = A/((r*r + r_2)*(r*r + r_2));
    points_values[i*3+0] = coeff*(2*(-r*y+x*z));
    points_values[i*3+1] = coeff*(2*(r*x+y*z));
    points_values[i*3+2] = coeff*(r*r - x*x - y*y + z*z);
  }
}

template<typename real_t, int sdim>
void
get_vel_field_hom_zero(const real_t* points_pos,
            int num_points,
            real_t* points_values) {
  for (int i = 0; i < num_points; i++) {
    points_values[i*3+0] = 0;
    points_values[i*3+1] = 0;
    points_values[i*3+2] = 0;
  }
}

template<typename real_t, int sdim>
void
get_vel_field_hom_x(const real_t* points_pos,
            int num_points,
            real_t* points_values) {
  real_t omega = -0.5;
  real_t time_factor = 1;//+sin(2*3.14159*time);
  for (int i = 0; i < num_points; i++) {
    points_values[i*3+0] = omega*time_factor;
    points_values[i*3+1] = 0;
    points_values[i*3+2] = 0;
  }
}

template<typename real_t, int sdim>
void
get_vel_field_hom_y(const real_t* points_pos,
            int num_points,
            real_t* points_values) {
  real_t omega = -0.5;
  real_t time_factor = 1;//+sin(2*3.14159*time);
  for (int i = 0; i < num_points; i++) {
    points_values[i*3+0] = 0;
    points_values[i*3+1] = omega*time_factor;
    points_values[i*3+2] = 0;
  }
}

template<typename real_t, int sdim>
void
get_slotted_cylinder(const real_t* points_pos,
                     int num_points,
                     real_t* out,
                     const real_t xc = 0.5,
                     const real_t yc = 0.5,
                     const real_t zc = 0.5,
                     const real_t R  = 0.3,
                     const real_t w  = 0.1,
                     const real_t a  = 0.0) {
  for (int i = 0; i < num_points; i++) {
    out[i] = 0;
    real_t r_2 = (points_pos[i*sdim+0]-xc) * (points_pos[i*sdim+0]-xc) +
        (points_pos[i*sdim+1]-yc) * (points_pos[i*sdim+1]-yc) +
        (points_pos[i*sdim+2]-zc) * (points_pos[i*sdim+2]-zc);
    if ( r_2 < R*R) {
      out[i] = 1;
      real_t d = cos(a)*(points_pos[i*sdim+1]-yc) - sin(a)*(points_pos[i*sdim+0]-xc);
      if ( fabs(d) <  w) {
        real_t x = sin(a)*(points_pos[i*sdim+1]-yc) + cos(a)*(points_pos[i*sdim+0]-xc);
        if( x + w > 0)
          out[i]  = 0;
      }
    }
  }
}

template<typename real_t, int sdim>
void
get_gaussian_field_cylinder(const real_t* points_pos,
                            int num_points,
                            real_t* out,
                            const real_t xc      = 0.6 ,
                            const real_t yc      = 0.5 ,
                            const real_t theta   = 0.0 ,
                            const real_t sigma_x = 0.06,
                            const real_t sigma_y = 0.06,
                            const real_t A       = 1.0) {
  real_t cos_theta   = cos(theta);
  real_t cos_2theta  = cos(2*theta);
  real_t cos_theta_2 = cos_theta * cos_theta;
  real_t sin_theta   = sin(theta);
  real_t sin_2theta  = sin(2*theta);
  real_t sin_theta_2 = sin_theta * sin_theta;
  real_t sigma_x_2   = sigma_x*sigma_x;
  real_t sigma_y_2   = sigma_y*sigma_y;

  real_t a = cos_theta_2*0.5/sigma_x_2 + sin_theta_2*0.5/sigma_y_2;
  real_t b = -sin_2theta*0.5/sigma_x_2 + sin_2theta*0.5/sigma_y_2;
  real_t c = sin_theta_2*0.5/sigma_x_2 + cos_theta_2*0.5/sigma_y_2;

  for (int i = 0; i < num_points; i++) {
    out[i]  =
        A*exp(-(a * (points_pos[i*sdim+0]-xc) * (points_pos[i*sdim+0]-xc) +
                b * (points_pos[i*sdim+0]-xc) * (points_pos[i*sdim+1]-yc) +
                c * (points_pos[i*sdim+1]-yc) * (points_pos[i*sdim+1]-yc)
                )
              );
  }
}


template<typename real_t, int sdim>
void
get_gaussian_field_3d(const real_t* points_pos,
                      const int num_points,
                      real_t* out,
                      const real_t xc      = 0.8 ,
                      const real_t yc      = 0.5 ,
                      const real_t zc      = 0.5 ,
                      const real_t A       = 1.0 ,
                      const real_t sigma_x = 0.06,
                      const real_t sigma_y = 0.35,
                      const real_t sigma_z = 0.35) {
  for (int i = 0; i < num_points; i++) {
    out[i]  = A*exp(-(((points_pos[i*sdim+0]-xc) * (points_pos[i*sdim+0]-xc))/(2*sigma_x*sigma_x) +
                      ((points_pos[i*sdim+1]-yc) * (points_pos[i*sdim+1]-yc))/(2*sigma_y*sigma_y) +
                      ((points_pos[i*sdim+2]-zc) * (points_pos[i*sdim+2]-zc))/(2*sigma_z*sigma_z)
                      )
                    );
  }
}

template <class Real_t>
void gaussian_kernel_diffusion_input(const Real_t* coord,
                                     int n,
                                     Real_t* out,
                                     const Real_t alpha,
                                     const Real_t xc = 0.5,
                                     const Real_t yc = 0.5,
                                     const Real_t zc = 0.5,
                                     const int a  = -160,
                                     const Real_t amp = 1.0) { //Input function
  int dof=1;
  for(int i=0;i<n;i++) {
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2=(c[0]-xc)*(c[0]-xc)+(c[1]-yc)*(c[1]-yc)+(c[2]-yc)*(c[2]-yc);
      out[i*dof+0]=-amp*(2*a*r_2+3)*2*a*exp(a*r_2)/alpha+exp(a*r_2);
    }
  }
}

// out = amp * exp(a*R^2)
template <class Real_t>
void gaussian_kernel(const Real_t* coord,
                     int n,
                     Real_t* out,
                     const Real_t xc = 0.5,
                     const Real_t yc = 0.5,
                     const Real_t zc = 0.5,
                     const int a  = -160,
                     const Real_t amp = 1.0) { //Output potential
  int dof=1;
  for(int i=0;i<n;i++) {
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2=(c[0]-xc)*(c[0]-xc)+(c[1]-yc)*(c[1]-yc)+(c[2]-zc)*(c[2]-zc);
      out[i*dof+0]=amp*exp(a*r_2);
    }
  }
}

// 1/(4*pi*t)^(3/2)*exp(-r^2/4*d*t)
template <class Real_t>
void diffusion_kernel(const Real_t* coord,
                      int n,
                      Real_t* out,
                      const Real_t diff_coeff,
                      const Real_t t,
                      const Real_t amp = 1e-2,
                      const Real_t xc = 0.5,
                      const Real_t yc = 0.5,
                      const Real_t zc = 0.5) { //Input function
  assert(t!=0);
  // const Real_t PI = 3.141592653589;
  int dof        = 1;
  Real_t dt4     = 4*diff_coeff*t;
  Real_t pidt4_3 = std::pow(PI*dt4, 3);
  Real_t nn      = 1.0/sqrt(pidt4_3);
  for(int i=0;i<n;i++) {
    const Real_t* c=&coord[i*COORD_DIM];
    {
      Real_t r_2   = (c[0]-xc)*(c[0]-xc)+(c[1]-yc)*(c[1]-yc)+(c[2]-zc)*(c[2]-zc);
      out[i*dof+0] = amp*nn*exp(-r_2/dt4);
    }
  }
}

}  // namespace tbslas

#endif // SRC_UTILS_FIELDS_H_
