#ifndef _HW_TOOL_H_
#define _HW_TOOL_H_

#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include "backward.hpp"
#include <cmath>
#include <State.h>
#include <vector>
#include <algorithm>

#define e_ 0.0001
#define MAX_ITR 200
#define T_start 100
#define SEARCH_STEP 0.1

class Homeworktool
{	
	private:

	protected:
		uint8_t * data;

		int GLX_SIZE, GLY_SIZE, GLZ_SIZE;
		int GLXYZ_SIZE, GLYZ_SIZE;

		double resolution, inv_resolution;
		double gl_xl, gl_yl, gl_zl;
		double gl_xu, gl_yu, gl_zu;	

		Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i & index);
		Eigen::Vector3i coord2gridIndex(const Eigen::Vector3d & pt);

		double sum_v_2;
		double sum_p_2;
		double sum_v_p;

	public:
		Homeworktool(){};
		~Homeworktool(){};

		void initGridMap(double _resolution, Eigen::Vector3d global_xyz_l, Eigen::Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id);
		void setObs(const double coord_x, const double coord_y, const double coord_z);
		bool isObsFree(const double coord_x, const double coord_y, const double coord_z);
				
		Eigen::Vector3d coordRounding(const Eigen::Vector3d & coord);
		double OptimalBVP(Eigen::Vector3d _start_position,Eigen::Vector3d _start_velocity,Eigen::Vector3d _target_position);

		double OptimalT_numerical(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v);
		double OptimalT_analytic(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v);
		double Calculate_J(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v, double T);
		double Calculate_dJ(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v, double T);
		double Calculate_dJ2(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v, double T);
};

#endif