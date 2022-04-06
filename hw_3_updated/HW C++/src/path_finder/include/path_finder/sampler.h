/*
Copyright (C) 2022 Hongkai Ye (kyle_yeh@163.com)
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
#ifndef _BIAS_SAMPLER_
#define _BIAS_SAMPLER_

#include <ros/ros.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <random>
#include <cmath>

class BiasSampler
{
public:
  BiasSampler()
  {
    std::random_device rd;
    gen_ = std::mt19937_64(rd());
    uniform_rand_ = std::uniform_real_distribution<double>(0.0, 1.0);
    normal_rand_ = std::normal_distribution<double>(0.0, 1.0);
    range_.setZero();
    origin_.setZero();
  };

  void setSamplingRange(const Eigen::Vector3d origin, const Eigen::Vector3d range)
  {
    origin_ = origin;
    range_ = range;
  }

  void samplingOnce(Eigen::Vector3d &sample)
  {
    sample[0] = uniform_rand_(gen_);
    sample[1] = uniform_rand_(gen_);
    sample[2] = uniform_rand_(gen_);
    // range : the range of the whole map
    sample.array() *= range_.array();
    sample += origin_;
  };

  // (0.0 - 1.0)
  double getUniRandNum()
  {
    return uniform_rand_(gen_);
  }

protected:
  Eigen::Vector3d range_, origin_;
  std::mt19937_64 gen_;
  std::uniform_real_distribution<double> uniform_rand_;
  std::normal_distribution<double> normal_rand_;
};

class InformedSampler : public BiasSampler
{
public:
  InformedSampler(const Eigen::Vector3d & start_pt, const Eigen::Vector3d & goal_pt): BiasSampler()
  {
    Eigen::Matrix<double, 1, 3> eZ_T;
    eZ_T << 0, 1, 0;
    centre_ = (start_pt + goal_pt) * 0.5;
    Eigen::Vector3d direction;
    direction = goal_pt - start_pt;
    dist_start_goal_ = direction.norm();
    direction = direction * (1 / dist_start_goal_);
    rotation_matrix_ = direction * eZ_T;
    radius_l_ = 4*dist_start_goal_;
    radius_s_ = 4*dist_start_goal_;
    range_ << radius_s_, radius_l_, radius_s_;
    origin_ = centre_;

  };

  // update the long and the short radius
  void updateRadius(double curr_best_length)
  {
    radius_l_ = 0.5 * curr_best_length;
    radius_s_ = 0.5 * sqrt( pow(curr_best_length, 2) - pow(dist_start_goal_, 2) );
  }

  void samplingInformed(Eigen::Vector3d &sample, const Eigen::Vector3d & start_point, const Eigen::Vector3d & goal_point, double dist_start_goal, double length_curr_path)
  {
    // sample in a sphere set at origin
    Eigen::Vector3d p;
    p << normal_rand_(gen_), normal_rand_(gen_), normal_rand_(gen_);
    double r = pow(uniform_rand_(gen_), 0.33333);
    sample = r * p.normalized();
    // apply transform
    sample.array() *= range_.array();
    sample = rotation_matrix_ * sample + centre_;
    
  }

private:
  Eigen::Vector3d centre_;
  Eigen::Matrix3d rotation_matrix_;
  double radius_l_;
  double radius_s_;
  double dist_start_goal_;
 
};

#endif
