/*******************************************************************************
 * BSD 3-Clause License
 *
 * Copyright (c) 2021, Qiayuan Liao
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * * Neither the name of the copyright holder nor the names of its
 *   contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *******************************************************************************/

//
// Created by qiayuan on 4/3/20.
//

#include <ros/ros.h>
#include <geometry_msgs/Twist.h>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/integrate_const.hpp>
#include <std_msgs/Bool.h>
#include <std_msgs/Float32.h>
#include <realtime_tools/realtime_buffer.h>
#include <geometry_msgs/TransformStamped.h>
#include <geometry_msgs/PointStamped.h>
#include <dynamic_reconfigure/server.h>
#include <rm_common/linear_interpolation.h>
#include <rm_common/ros_utilities.h>
#include <rm_common/ori_tool.h>
#include <rm_msgs/TrackData.h>
#include <std_msgs/Float64.h>

struct BallisticConfig
{
  double mass, radius, gun_len, drag_coff, Cd, air_density, g;
  double initial_vel_near, initial_vel_far, max_simulation_time, max_integration_step;
  double newton_convergence_tol, finite_difference_eps, max_newton_step;
  int max_newton_iterations;
};

class BallisticSolver
{
public:
  explicit BallisticSolver(ros::NodeHandle& nh);
  ~BallisticSolver() = default;
  bool used_fallback_;
  double simulate(double pitch_angle_, double initial_vel_, double target_dis, double target_hgt);
  bool solver(double& yaw, double& pitch);

private:
  BallisticConfig config_{};
  realtime_tools::RealtimeBuffer<BallisticConfig> config_rt_buffer_;
  rm_common::LinearInterp output_pitch_match_lut_;
  typedef boost::numeric::odeint::runge_kutta4<std::array<double, 6>> stepper_;
};
BallisticSolver::BallisticSolver(ros::NodeHandle& controller_nh) : used_fallback_(false)
{
  config_ = { .mass = getParam(controller_nh, "mass", 0.0445),
              .radius = getParam(controller_nh, "radius", 0.02125),
              .gun_len = getParam(controller_nh, "gun_len", 0.08),
              .drag_coff= getParam(controller_nh, "drag_coff", 1.2),
              .Cd = getParam(controller_nh, "Cd", 0.4),
              .air_density = getParam(controller_nh, "air_density", 1.2),
              .g = getParam(controller_nh, "g", 9.81),
              .initial_vel_near = getParam(controller_nh, "initial_vel_near", 14.0),
              .initial_vel_far = getParam(controller_nh, "initial_vel_far", 15.7),
              .max_simulation_time = getParam(controller_nh, "max_simulation_time", 2.5),
              .max_integration_step = getParam(controller_nh, "max_integration_step", 0.01),
              .newton_convergence_tol = getParam(controller_nh, "newton_convergence_tol", 2e-5),
              .finite_difference_eps = getParam(controller_nh, "finite_difference_eps", 1e-4),
              .max_newton_step = getParam(controller_nh, "max_newton_step", 0.04),
              .max_newton_iterations = getParam(controller_nh, "max_newton_iterations", 2)};
  config_rt_buffer_.initRT(config_);
  XmlRpc::XmlRpcValue lut_config;
  if (controller_nh.getParam("/controllers/gimbal_controller/ballistic_solver/output_pitch_match", lut_config))
    output_pitch_match_lut_.init(lut_config);
}
bool BallisticSolver::solver(double& yaw, double& pitch)
{
  geometry_msgs::Vector3 launch2target;
  launch2target.x = 19;
  launch2target.y = 6;
  launch2target.z = 1.2 - 0.513468 +0.05;
  double target_dis = std::sqrt(launch2target.x * launch2target.x + launch2target.y * launch2target.y) - config_.gun_len;
  double target_hgt = launch2target.z;
  double initial_vel = target_dis <= 16.5 ? config_.initial_vel_near : config_.initial_vel_far;
  yaw = std::atan2(launch2target.y, launch2target.x);
  // Use LUT to get initial guess.Originally defined: negative pitch indicates upward angle (head up)
  double initial_pitch = -output_pitch_match_lut_.output(target_dis);
  std::cout << "initial_pitch: " << initial_pitch << std::endl;
  // Error function for root finding
  auto error_function = [&](double pitch_angle) -> double {
    return simulate(pitch_angle, initial_vel, target_dis, target_hgt);
  };
  used_fallback_ = true;
  double current_pitch = initial_pitch;
  double error = error_function(current_pitch);
  // Check if initial guess is already good enough
  if (std::abs(error) < config_.newton_convergence_tol)
  {
    std::cout<<"initial pitch is well"<<std::endl;
    pitch = -initial_pitch;
    used_fallback_ = false;
    return true;
  }
  // Newton iteration
  for (int iter = 0; iter < config_.max_newton_iterations; ++iter)
  {
    double h = config_.finite_difference_eps;
    double f_plus  = error_function(current_pitch + h);
    // double f_minus = error_function(current_pitch - h);
    // double jacobian = (f_plus - f_minus) / (2 * h);
    double jacobian= (f_plus - error) / h;
    double delta = error / jacobian;
    delta = (delta < -config_.max_newton_step) ? -config_.max_newton_step :
            (delta > config_.max_newton_step)  ? config_.max_newton_step :
                                                 delta;
    double update_pitch = current_pitch - delta;
    double update_error = error_function(update_pitch);
    // Check convergence after update
    if (std::abs(update_error) < config_.newton_convergence_tol)
    {
      pitch = -update_pitch;
      used_fallback_ = false;
      return true;
    }
    // Update for next iteration
    current_pitch = update_pitch;
    error = update_error;
    std::cout<<"newton"<<std::endl;
  }
  pitch = -initial_pitch;
  return true;
}

double BallisticSolver::simulate(double pitch_angle, double initial_vel, double target_dis, double target_hgt)
{
  if (initial_vel <= 0.0 || target_dis <= 0.0)
    return 1e6;
  BallisticConfig config = *config_rt_buffer_.readFromRT();
  std::array<double, 6> state = {{0.0, 0.0, 0.0,initial_vel * std::cos(pitch_angle),0.0,initial_vel * std::sin(pitch_angle)}};
  double t = 0.0;
  double t_max = config.max_simulation_time;
  double dt = config.max_integration_step;
  stepper_ stepper;
  double x_prev = 0.0, z_prev = 0.0;
  double x_curr = 0.0, z_curr = 0.0;
  double z_at_target = 0.0;
  auto system_func = [config](const std::array<double, 6>& state, std::array<double, 6>& dsdt, double t) {
    double vx = state[3], vy = state[4], vz = state[5];
    double speed_sq = vx*vx + vy*vy + vz*vz;
    dsdt[0] = vx;
    dsdt[1] = vy;
    dsdt[2] = vz;
    double ax = 0.0, ay = 0.0, az = -config.g;
    if (speed_sq > 1e-10) {
      double speed = std::sqrt(speed_sq);
      double F_drag =
          config.drag_coff * 0.5 * config.air_density * config.Cd * M_PI * config.radius * config.radius * speed_sq;
      double inv_speed = 1.0 / speed;
      ax += (-F_drag * vx * inv_speed) / config.mass;
      ay += (-F_drag * vy * inv_speed) / config.mass;
      az += (-F_drag * vz * inv_speed) / config.mass;
    }
    dsdt[3] = ax;
    dsdt[4] = ay;
    dsdt[5] = az;
  };
  while (t < t_max && state[0] < target_dis)
  {
    x_prev = x_curr;
    z_prev = z_curr;
    x_curr = state[0];
    z_curr = state[2];
    stepper.do_step(system_func, state, t, dt);
    t += dt;
  }
  if (state[0] >= target_dis) {
    z_at_target = z_prev + (z_curr - z_prev) * (target_dis - x_prev) / (x_curr - x_prev);
    std::cout<<t<<std::endl;
  }
  std::cout << "z_at_target = " << z_at_target << std::endl;
  return z_at_target - target_hgt;
}
int main(int argc, char** argv)
{
  ros::init(argc, argv, "ballistic_solver");
  ros::NodeHandle nh;
  BallisticSolver solver(nh);
  ros::Rate rate(10);
  while (ros::ok())
  {
    // ros::Publisher pub = nh.advertise<std_msgs::Float64>("/ballistic_solver_result", 10);
    // std_msgs::Float64 msg;
    // pub.publish(msg);
    double yaw, pitch;
    bool success = solver.solver(yaw, pitch);
    if (success)
      std::cout << "yaw = " << yaw << ", pitch = " << pitch << std::endl;
    rate.sleep();
  }
  return 0;
}