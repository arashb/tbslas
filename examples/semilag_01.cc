// Copyright (C) 2014 by Arash Bakhtiari

#include <vector>
#include <cmath>

#include "semilag/utils.h"
#include "semilag/cubic_interp_policy.h"
#include "semilag/vec_field.h"
#include "semilag/semilag.h"

typedef slas::VecField<double,3,3> VFieldD33;
typedef slas::VecField<double,3,1> VFieldD31;
typedef std::vector<double> VecD;

int main(int argc, char *argv[]) {
  VFieldD33 vel_field;
  VFieldD31 con_field;

  size_t domain_res     = 16;
  size_t num_points_dim = domain_res+1;
  double dx             = 1.0/domain_res;
  double dt             = dx;
  int tn                = 10;
  double tf             = tn*dt;
  int num_rk_step       = 1;

  const size_t num_points            = (size_t)std::pow(num_points_dim, 3);
  const VecD points_pos              = slas::generate_reg_grid_points<double,3>(num_points_dim);
  const VecD points_values_vorticity = slas::generate_vorticity_field(points_pos);
  const VecD points_values_gaussian  = slas::generate_gaussian_field(points_pos);

  slas::CubicInterpPolicy<double> cubic_interp_policy;

  vel_field.init(points_pos, points_values_vorticity);
  con_field.init(points_pos, points_values_gaussian);

  for(int timestep = 1; timestep <= tn; timestep++) {
    // vel_field.push_back_values(slas::generate_vorticity_field(points_pos, timestep*dt),
    //                            timestep*dt);
    vel_field.push_back_values(slas::generate_vorticity_field(points_pos),
                               timestep*dt);
  }

  for(int timestep = 1; timestep <= tn; timestep++) {
    std::cout << "****************************************" << std::endl;
    std::cout << "-> timestep: " << timestep
              << " time: " << dt*timestep << std::endl;

    slas::semilag_rk2(vel_field,
                      cubic_interp_policy,
                      timestep,
                      dt,
                      num_rk_step,
                      con_field);
  }

  int xyz_mem_layout = 0;
  vel_field.save("semilag_01_vel_", xyz_mem_layout);
  con_field.save("semilag_01_con_", xyz_mem_layout);

  return 0;
}
