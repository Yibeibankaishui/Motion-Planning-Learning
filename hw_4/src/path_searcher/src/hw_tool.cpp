#include <hw_tool.h>

using namespace std;
using namespace Eigen;

void Homeworktool::initGridMap(double _resolution, Vector3d global_xyz_l, Vector3d global_xyz_u, int max_x_id, int max_y_id, int max_z_id)
{   
    gl_xl = global_xyz_l(0);
    gl_yl = global_xyz_l(1);
    gl_zl = global_xyz_l(2);

    gl_xu = global_xyz_u(0);
    gl_yu = global_xyz_u(1);
    gl_zu = global_xyz_u(2);
    
    GLX_SIZE = max_x_id;
    GLY_SIZE = max_y_id;
    GLZ_SIZE = max_z_id;
    GLYZ_SIZE  = GLY_SIZE * GLZ_SIZE;
    GLXYZ_SIZE = GLX_SIZE * GLYZ_SIZE;

    resolution = _resolution;
    inv_resolution = 1.0 / _resolution;    

    data = new uint8_t[GLXYZ_SIZE];
    memset(data, 0, GLXYZ_SIZE * sizeof(uint8_t));
}

void Homeworktool::setObs(const double coord_x, const double coord_y, const double coord_z)
{   
    if( coord_x < gl_xl  || coord_y < gl_yl  || coord_z <  gl_zl || 
        coord_x >= gl_xu || coord_y >= gl_yu || coord_z >= gl_zu )
        return;

    int idx_x = static_cast<int>( (coord_x - gl_xl) * inv_resolution);
    int idx_y = static_cast<int>( (coord_y - gl_yl) * inv_resolution);
    int idx_z = static_cast<int>( (coord_z - gl_zl) * inv_resolution);      
    
    data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] = 1;
}

bool Homeworktool::isObsFree(const double coord_x, const double coord_y, const double coord_z)
{
    Vector3d pt;
    Vector3i idx;
    
    pt(0) = coord_x;
    pt(1) = coord_y;
    pt(2) = coord_z;
    idx = coord2gridIndex(pt);

    int idx_x = idx(0);
    int idx_y = idx(1);
    int idx_z = idx(2);

    return (idx_x >= 0 && idx_x < GLX_SIZE && idx_y >= 0 && idx_y < GLY_SIZE && idx_z >= 0 && idx_z < GLZ_SIZE && 
           (data[idx_x * GLYZ_SIZE + idx_y * GLZ_SIZE + idx_z] < 1));
}

Vector3d Homeworktool::gridIndex2coord(const Vector3i & index) 
{
    Vector3d pt;

    pt(0) = ((double)index(0) + 0.5) * resolution + gl_xl;
    pt(1) = ((double)index(1) + 0.5) * resolution + gl_yl;
    pt(2) = ((double)index(2) + 0.5) * resolution + gl_zl;

    return pt;
}

Vector3i Homeworktool::coord2gridIndex(const Vector3d & pt) 
{
    Vector3i idx;
    idx <<  min( max( int( (pt(0) - gl_xl) * inv_resolution), 0), GLX_SIZE - 1),
            min( max( int( (pt(1) - gl_yl) * inv_resolution), 0), GLY_SIZE - 1),
            min( max( int( (pt(2) - gl_zl) * inv_resolution), 0), GLZ_SIZE - 1);                  
  
    return idx;
}

Eigen::Vector3d Homeworktool::coordRounding(const Eigen::Vector3d & coord)
{
    return gridIndex2coord(coord2gridIndex(coord));
}

double Homeworktool::OptimalBVP(Eigen::Vector3d _start_position,Eigen::Vector3d _start_velocity,Eigen::Vector3d _target_position)
{
    double optimal_cost = 100000; // this just to initial the optimal_cost, you can delete it 
    /*
                    



    STEP 2: go to the hw_tool.cpp and finish the function Homeworktool::OptimalBVP
    the solving process has been given in the document

    because the final point of trajectory is the start point of OBVP, so we input the pos,vel to the OBVP

    after finish Homeworktool::OptimalBVP, the Trajctory_Cost will record the optimal cost of this trajectory


    */

    Eigen::Vector3d d_p = _target_position - _start_position;
    Eigen::Vector3d d_v = 0 - _start_velocity;

    double optimal_time = ;
    optimal_time = 
    optimal_cost = 0;


    return optimal_cost;
}


// numerical method
double OptimalT_numerical(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v){

}


double OptimalT_analytic(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v){
    double a = 16 * (d_v.array().pow(2).sum());
    double b = -48 * (d_p.transpose() * d_v);
    double c = 36 * (d_p.array().pow(2).sum()); 
}


double Calculate_J(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v, double T){
    double p3 = pow((3333*((12*d_p(0))/(pow(T, 3))) - (6*d_v(0))/(pow(T, 2))), 2)/10000 + pow((3333*((12*d_p(1))/(pow(T, 3))) - (6*d_v(1))/(pow(T, 2))), 2)/10000 + pow((3333*((12*d_p(2))/(pow(T, 3))) - (6*d_v(2))/(pow(T, 2))), 2)/10000;
    double p2 = -((6*d_p(0))/(pow(T, 2)) - (2*d_v(0))/T)*((12*d_p(0))/(pow(T, 3))) - (6*d_v(0))/(pow(T, 2))) - ((6*d_p(1))/(pow(T, 2)) - (2*d_v(1))/T)*((12*d_p(1))/(pow(T, 3))) - (6*d_v(1))/(pow(T, 2))) - ((6*d_p(2))/(pow(T, 2)) - (2*d_v(2))/T)*((12*d_p(2))/(pow(T, 3))) - (6*d_v(2))/(pow(T, 2)));
    double p1 = pow(((6*d_p(0))/(pow(T, 2)) - (2*d_v(0))/T), 2) + pow(((6*d_p(1))/(pow(T, 2)) - (2*d_v(1))/T),2) + pow(((6*d_p(2))/(pow(T, 2)) - (2*d_v(2))/T),2) + 1;
    double p0 = 0;
    double J = p3 * T * T * T + p2 * T * T + p1 * T;

}


double Calculate_dJ(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v, double T){
    double p4 = 1;
    double p3 = 0;
    double p2 = 16 * (d_v.array().pow(2).sum());
    double p1 = -48 * (d_p.transpose() * d_v);
    double p0 = 36 * (d_p.array().pow(2).sum()); 
    double dJ = p4 * pow(T, 4) + p3 * pow(T, 3) + p2 * pow(T, 2) + p1 * T + p0;
}


double Calculate_dJ2(const Eigen::Vector3d d_p, const Eigen::Vector3d d_v, double T){
    double p3 = 4;
    double p2 = 0;
    double p1 = 32 * (d_v.array().pow(2).sum());
    double p0 = -96 * (d_p.transpose() * d_v);
    double dJ2 = p3 * pow(T, 3) + p2 * pow(T, 2) + p1 * T + p0;
}

 
(3333*((12*d_p(0))/(pow(T, 3))) - (6*d_v(0))/(pow(T, 2)))^2)/10000 + (3333*((12*d_p(1))/(pow(T, 3))) - (6*d_v(1))/(pow(T, 2)))^2)/10000 + (3333*((12*d_p(2))/(pow(T, 3))) - (6*d_v(2))/(pow(T, 2)))^2)/10000,
 - ((6*d_p(0))/(pow(T, 2)) - (2*d_v(0))/T)*((12*d_p(0))/(pow(T, 3))) - (6*d_v(0))/(pow(T, 2))) - ((6*d_p(1))/(pow(T, 2)) - (2*d_v(1))/T)*((12*d_p(1))/(pow(T, 3))) - (6*d_v(1))/(pow(T, 2))) - ((6*d_p(2))/(pow(T, 2)) - (2*d_v(2))/T)*((12*d_p(2))/(pow(T, 3))) - (6*d_v(2))/(pow(T, 2))),
  ((6*d_p(0))/(pow(T, 2)) - (2*d_v(0))/T)^2 + ((6*d_p(1))/(pow(T, 2)) - (2*d_v(1))/T)^2 + ((6*d_p(2))/(pow(T, 2)) - (2*d_v(2))/T)^2 + 1,
