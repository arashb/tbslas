// Copyright (C) 2014 by Arash Bakhtiari

#include <vector>

#include "semilag/utils.h"
#include "semilag/cubic_interp_policy.h"
#include "semilag/vec_field.h"
#include "semilag/semilag.h"

typedef tbslas::VecField<double,3,3> VFieldD33 ;
typedef tbslas::VecField<double,3,1> VFieldD31 ;
typedef std::vector<double> VecD               ;

int main(int argc, char *argv[]) {
  VFieldD33 vel_field;
  VFieldD31 con_field;
  size_t domain_res = 16;
  size_t num_points_dim = domain_res+1;
  double dx         = 1.0/domain_res;
  double dt         = dx;
  int tn            = 10;
  double tf         = tn*dt;
  int num_rk_step   = 1;

  const size_t num_points            = num_points_dim*num_points_dim*num_points_dim;
  const VecD points_pos              = tbslas::generate_reg_grid_points<double,3>(num_points_dim);
  const VecD points_values_vorticity = tbslas::generate_vorticity_field(points_pos);
  const VecD points_values_gaussian  = tbslas::generate_gaussian_field(points_pos);
  tbslas::CubicInterpPolicy<double> cubic_interp_policy;

  vel_field.init(points_pos, points_values_vorticity);
  con_field.init(points_pos, points_values_gaussian);

  for(int timestep = 1; timestep <= tn; timestep++) {
    vel_field.push_back_values(tbslas::generate_vorticity_field(points_pos, timestep*dt),
                               timestep*dt);
  }

  for(int timestep = 1; timestep <= tn; timestep++) {
    std::cout << "****************************************" << std::endl;
    std::cout << "-> timestep: " << timestep
              << " time: " << dt*timestep << std::endl;

    tbslas::semilag_rk2(vel_field,
                        cubic_interp_policy,
                        timestep,
                        dt,
                        num_rk_step,
                        con_field
                        );
  }

  vel_field.save("semilag_vel_");
  con_field.save("semilag_con_");

  return 0;
}
