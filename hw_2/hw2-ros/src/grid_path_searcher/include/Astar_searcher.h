#ifndef _ASTART_SEARCHER_H
#define _ASTART_SEARCHER_H

#include <iostream>
#include <cmath>
#include <algorithm>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include "backward.hpp"
#include "node.h"

#define _use_tie_breaker 0


enum HeuFunc {
	Manhattan,
	Euclidean,
	DiagonalHeu,
	Dijkstra
};


struct HeuValue{
	double gValue;
	double hValue;
	double fValue;

	HeuValue(){}
	HeuValue(double g, double h) : gValue(g), hValue(h) {fValue = g + h;}
	HeuValue(const GridNodePtr & node){
		gValue = node -> gScore;
		hValue = node -> fScore;
		fValue = gValue + hValue;
	}

	bool operator== (HeuValue & rhs) const{
		return ( (gValue == rhs.gValue) && (fValue == rhs.fValue) );
	}
};

struct Cmp{
	bool operator () (const HeuValue & a, const HeuValue & b) const{
		if ((a.fValue) != (b.fValue)) {return a.fValue < b.fValue;}
		else return a.hValue < b.hValue;
	}
};


class AstarPathFinder
{	
	private:

	protected:
		uint8_t * data;
		GridNodePtr *** GridNodeMap;
		Eigen::Vector3i goalIdx;
		int GLX_SIZE, GLY_SIZE, GLZ_SIZE;
		int GLXYZ_SIZE, GLYZ_SIZE;

		double resolution, inv_resolution;
		double gl_xl, gl_yl, gl_zl;
		double gl_xu, gl_yu, gl_zu;

		GridNodePtr terminatePtr;

#if _use_tie_breaker
		std::multimap<HeuValue, GridNodePtr, Cmp> openSet;
#else
		std::multimap<double, GridNodePtr> openSet;
#endif
		double getHeu(GridNodePtr node1, GridNodePtr node2, HeuFunc heu=Euclidean);
		void AstarGetSucc(GridNodePtr currentPtr, std::vector<GridNodePtr> & neighborPtrSets, std::vector<double> & edgeCostSets);		

    	bool isOccupied(const int & idx_x, const int & idx_y, const int & idx_z) const;
		bool isOccupied(const Eigen::Vector3i & index) const;
		bool isFree(const int & idx_x, const int & idx_y, const int & idx_z) const;
		bool isFree(const Eigen::Vector3i & index) const;
		
		Eigen::Vector3d gridIndex2coord(const Eigen::Vector3i & index);
		// 坐标转化为索引
		Eigen::Vector3i coord2gridIndex(const Eigen::Vector3d & pt);

	public:
		AstarPathFinder(){};
		~AstarPathFinder(){};
		void AstarGraphSearch(Eigen::Vector3d start_pt, Eigen::Vector3d end_pt);
		void resetGrid(GridNodePtr ptr);
		void resetUsedGrids();

		void initGridMap(double _resolution, Eigen::Vector3d global_xyz_l, Eigen::Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id);
		void setObs(const double coord_x, const double coord_y, const double coord_z);

		Eigen::Vector3d coordRounding(const Eigen::Vector3d & coord);
		std::vector<Eigen::Vector3d> getPath();
		std::vector<Eigen::Vector3d> getVisitedNodes();
};



#endif