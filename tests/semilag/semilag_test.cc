// Copyright (C) 2014 by Arash Bakhtiari
#include <vector>
#include <cmath>

#include <gtest/gtest.h>

#include "semilag/utils.h"
#include "semilag/cubic_interp_policy.h"
#include "semilag/vec_field.h"
#include "semilag/semilag.h"

typedef slas::VecField<double,3,3> VFieldD33;
typedef slas::VecField<double,3,1> VFieldD31;
typedef std::vector<double> VecD              ;

#define PI 3.14159

class SemilagTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
  }

  const size_t domain_res      = 16;
  const size_t num_points_dim  = domain_res+1;
  const int cfl                = 1;
  const double dx              = 1.0/domain_res;
  const double dt              = cfl*dx;

  // const double tf              = 2*PI;
  // const int tn                 = tf/dt;

  const int tn                 = 1;
  const double tf              = tn*dt;

  const int num_rk_step        = 1;

  const int sdim               = 3;
  const size_t num_points      = (size_t)std::pow(num_points_dim, sdim);
  const size_t num_points_size = num_points*sdim;

};

TEST_F(SemilagTest, SteadyVorticityWithFunctor) {
  VecD points_val(num_points);  // one dimensional concentration values

  // init position of grid points
  VecD points_pos(num_points_size);
  slas::get_reg_grid_points<double,3>(num_points_dim, points_pos.data());

  // init velocity values
  VecD points_val_vorticity(num_points_size);
  slas::get_vorticity_field<double,3>(points_pos.data(), num_points,
                                      points_val_vorticity.data());

  VecD points_val_gaussian(num_points);
  slas::get_gaussian_field<double,3>(points_pos.data(), num_points,
                                     points_val_gaussian.data());

  VFieldD33 vel_field;
  VFieldD31 con_field;

  vel_field.init(points_pos, points_val_vorticity);
  con_field.init(points_pos, points_val_gaussian);

  for(int timestep = 1; timestep <= tn; timestep++) {
    slas::get_vorticity_field<double,3>(points_pos.data(),
                                        num_points,
                                        points_val_vorticity.data()
                                        );
    vel_field.push_back_values(points_val_vorticity,
                               timestep*dt);
  }

  for(int timestep = 1; timestep <= tn; timestep++) {
    std::cout << "****************************************" << std::endl;
    std::cout << "-> timestep: " << timestep
              << " time: " << dt*timestep << std::endl;

    slas::semilag_rk2(slas::get_vorticity_field<double,3>,
                      slas::get_gaussian_field<double,3>,
                      points_pos,
                      sdim,
                      timestep,
                      dt,
                      num_rk_step,
                      points_val
                      );

    // vel_field.push_back_values(points_val_vorticity, dt*timestep);
    con_field.push_back_values(points_val, dt*timestep);
  }

  vel_field.save("semilag_functor_vel_");
  con_field.save("semilag_functor_con_");

  // VFieldD33 vel_field_expected;
  // VFieldD31 con_field_expected;

  // for (int i = 0; i < num_points; i++) {
  //   std::cout << "expected: " << expected_val[i] << " computed: " << points_val[1] << std::endl;
  //   ASSERT_NEAR(expected_val[i], points_val[i], 0.001);
  // }
}

TEST_F(SemilagTest, SteadyVorticityWithVecField) {

}
