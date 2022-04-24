#include "lec5_hw/visualizer.hpp"
#include "lec5_hw/trajectory.hpp"

#include <ros/ros.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>

#include <cmath>
#include <iostream>
#include <vector>

struct Config
{
    std::string targetTopic;
    double clickHeight;
    std::vector<double> initialVel;
    std::vector<double> initialAcc;
    std::vector<double> terminalVel;
    std::vector<double> terminalAcc;
    double allocationSpeed;
    double allocationAcc;
    int maxPieceNum;

    Config(const ros::NodeHandle &nh_priv)
    {
        nh_priv.getParam("TargetTopic", targetTopic);
        nh_priv.getParam("ClickHeight", clickHeight);
        nh_priv.getParam("InitialVel", initialVel);
        nh_priv.getParam("InitialAcc", initialAcc);
        nh_priv.getParam("TerminalVel", terminalVel);
        nh_priv.getParam("TerminalAcc", terminalAcc);
        nh_priv.getParam("AllocationSpeed", allocationSpeed);
        nh_priv.getParam("AllocationAcc", allocationAcc);
        nh_priv.getParam("MaxPieceNum", maxPieceNum);
    }
};

double timeTrapzVel(const double dist,
                    const double vel,
                    const double acc)
{
    const double t = vel / acc;
    const double d = 0.5 * acc * t * t;

    if (dist < d + d)
    {
        return 2.0 * sqrt(dist / acc);
    }
    else
    {
        return 2.0 * t + (dist - 2.0 * d) / vel;
    }
}

void minimumJerkTrajGen(
    // Inputs:
    const int pieceNum,
    const Eigen::Vector3d &initialPos,
    const Eigen::Vector3d &initialVel,
    const Eigen::Vector3d &initialAcc,
    const Eigen::Vector3d &terminalPos,
    const Eigen::Vector3d &terminalVel,
    const Eigen::Vector3d &terminalAcc,
    const Eigen::Matrix3Xd &intermediatePositions,
    const Eigen::VectorXd &timeAllocationVector,    //Time allocated for each piece
    // Outputs:
    Eigen::MatrixX3d &coefficientMatrix)
{
    // coefficientMatrix is a matrix with 6*piece num rows and 3 columes
    // As for a polynomial c0+c1*t+c2*t^2+c3*t^3+c4*t^4+c5*t^5,
    // each 6*3 sub-block of coefficientMatrix is
    // --              --
    // | c0_x c0_y c0_z |
    // | c1_x c1_y c1_z |
    // | c2_x c2_y c2_z |
    // | c3_x c3_y c3_z |
    // | c4_x c4_y c4_z |
    // | c5_x c5_y c5_z |
    // --              --
    // Please computed coefficientMatrix of the minimum-jerk trajectory
    // in this function

    // ------------------------ Put your solution below ------------------------

    // T, T^1, T^2, T^3, T^4, T^5
    const Eigen::VectorXd T1 = timeAllocationVector;
    const Eigen::VectorXd T2 = timeAllocationVector.array().pow(2);
    const Eigen::VectorXd T3 = timeAllocationVector.array().pow(3);
    const Eigen::VectorXd T4 = timeAllocationVector.array().pow(4);
    const Eigen::VectorXd T5 = timeAllocationVector.array().pow(5);

    const int N = pieceNum;
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(6 * N, 6 * N);
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(6 * N, 3);

    b.row(0) = initialPos.transpose();
    b.row(1) = initialVel.transpose();
    b.row(2) = initialAcc.transpose();

    b.row(6 * N - 3) = terminalPos.transpose();
    b.row(6 * N - 2) = terminalVel.transpose();
    b.row(6 * N - 1) = terminalAcc.transpose();

    // initial constraint
    M(0, 0) = 1.0;
    M(1, 1) = 1.0;
    M(2, 2) = 2.0;
    // terminal constraint
    // pos
    M(6 * N - 3, 6 * N - 6) = 1.0;
    M(6 * N - 3, 6 * N - 5) = T1(N - 1);
    M(6 * N - 3, 6 * N - 4) = T2(N - 1);
    M(6 * N - 3, 6 * N - 3) = T3(N - 1);
    M(6 * N - 3, 6 * N - 2) = T4(N - 1);
    M(6 * N - 3, 6 * N - 1) = T5(N - 1);
    // velocity
    M(6 * N - 2, 6 * N - 5) = 1.0;
    M(6 * N - 2, 6 * N - 4) = 2.0 * T1(N - 1);
    M(6 * N - 2, 6 * N - 3) = 3.0 * T2(N - 1);
    M(6 * N - 2, 6 * N - 2) = 4.0 * T3(N - 1);
    M(6 * N - 2, 6 * N - 1) = 5.0 * T4(N - 1);
    // acceleration
    M(6 * N - 1, 6 * N - 4) = 2.0;
    M(6 * N - 1, 6 * N - 3) = 6.0 * T1(N - 1);
    M(6 * N - 1, 6 * N - 2) = 12.0 * T2(N - 1);
    M(6 * N - 1, 6 * N - 1) = 20.0 * T3(N - 1);    

    // construct the matrix M
    for (int i = 0; i < N - 1; i++){

        // continuity constraints
        // f_j(T) == f_{j+1}(0)

        // Jerk
        // jerk = 6 c3 + 24 t c4 + 60 t2 c5
        M(6 * i + 3, 6 * i + 3) = 6.0;
        M(6 * i + 3, 6 * i + 4) = 24.0 * T1(i);
        M(6 * i + 3, 6 * i + 5) = 60.0 * T2(i);
        // jerk of next piece at time 0
        M(6 * i + 3, 6 * i + 6) = -6.0;

        // Snap
        // snap = 24 c4 + 120 t c5
        M(6 * i + 4, 6 * i + 4) = 24.0;
        M(6 * i + 4, 6 * i + 5) = 120.0 * T1(i);
        // snap of next piece at time 0
        M(6 * i + 4, 6 * i + 6) = -24.0;

        // Position
        // Pos(T) = waypoints(i+1)
        M(6 * i + 5, 6 * i) = 1.0;
        M(6 * i + 5, 6 * i + 1) = T1(i);
        M(6 * i + 5, 6 * i + 2) = T2(i);
        M(6 * i + 5, 6 * i + 3) = T3(i);
        M(6 * i + 5, 6 * i + 4) = T4(i);
        M(6 * i + 5, 6 * i + 5) = T5(i);
        // next waypoint
        b.row(6 * i + 5) = intermediatePositions.col(i).transpose();

        // Velocity
        M(6 * i + 6, 6 * i + 1) = 1.0;
        M(6 * i + 6, 6 * i + 2) = 2.0 * T1(i);
        M(6 * i + 6, 6 * i + 3) = 3.0 * T2(i);
        M(6 * i + 6, 6 * i + 4) = 4.0 * T3(i);
        M(6 * i + 6, 6 * i + 5) = 5.0 * T4(i);
        // velocity of next piece at time 0
        M(6 * i + 6, 6 * i + 6) = -1.0;

        // Acceleration
        M(6 * i + 7, 6 * i + 2) = 2.0;
        M(6 * i + 7, 6 * i + 3) = 6.0 * T1(i);
        M(6 * i + 7, 6 * i + 4) = 12.0 * T2(i);
        M(6 * i + 7, 6 * i + 5) = 20.0 * T3(i);
        // acceleration of next piece at time 0
        M(6 * i + 7, 6 * i + 6) = -2.0;

        // Position
        // Pos_j(T) = Pos_j+1(0)
        M(6 * i + 8, 6 * i) = 1.0;
        M(6 * i + 8, 6 * i + 1) = T1(i);
        M(6 * i + 8, 6 * i + 2) = T2(i);
        M(6 * i + 8, 6 * i + 3) = T3(i);
        M(6 * i + 8, 6 * i + 4) = T4(i);
        M(6 * i + 8, 6 * i + 5) = T5(i);
        M(6 * i + 8, 6 * i + 6) = -1.0;
    }


    // Solve the equation
    coefficientMatrix = M.fullPivLu().solve(b);

    // ------------------------ Put your solution above ------------------------
}

class ClickGen
{
private:
    Config config;

    ros::NodeHandle nh;
    ros::Subscriber targetSub;

    Visualizer visualizer;

    Eigen::Matrix3Xd positions;
    Eigen::VectorXd times;
    int positionNum;
    Trajectory<5> traj;

public:
    ClickGen(const Config &conf,
             ros::NodeHandle &nh_)
        : config(conf),
          nh(nh_),
          visualizer(nh),
          positions(3, config.maxPieceNum + 1),
          times(config.maxPieceNum),
          positionNum(0)
    {
        targetSub = nh.subscribe(config.targetTopic, 1,
                                 &ClickGen::targetCallBack, this,
                                 ros::TransportHints().tcpNoDelay());
    }

    void targetCallBack(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        if (positionNum > config.maxPieceNum)
        {
            positionNum = 0;
            traj.clear();
        }

        positions(0, positionNum) = msg->pose.position.x;
        positions(1, positionNum) = msg->pose.position.y;
        positions(2, positionNum) = std::fabs(msg->pose.orientation.z) * config.clickHeight;

        if (positionNum > 0)
        {
            const double dist = (positions.col(positionNum) - positions.col(positionNum - 1)).norm();
            times(positionNum - 1) = timeTrapzVel(dist, config.allocationSpeed, config.allocationAcc);
        }

        ++positionNum;

        if (positionNum > 1)
        {
            const int pieceNum = positionNum - 1;
            const Eigen::Vector3d initialPos = positions.col(0);
            const Eigen::Vector3d initialVel(config.initialVel[0], config.initialVel[1], config.initialVel[2]);
            const Eigen::Vector3d initialAcc(config.initialAcc[0], config.initialAcc[1], config.initialAcc[2]);
            const Eigen::Vector3d terminalPos = positions.col(pieceNum);
            const Eigen::Vector3d terminalVel(config.terminalVel[0], config.terminalVel[1], config.terminalVel[2]);
            const Eigen::Vector3d terminalAcc(config.terminalAcc[0], config.terminalAcc[1], config.terminalAcc[2]);
            const Eigen::Matrix3Xd intermediatePositions = positions.middleCols(1, pieceNum - 1);
            const Eigen::VectorXd timeAllocationVector = times.head(pieceNum);

            Eigen::MatrixX3d coefficientMatrix = Eigen::MatrixXd::Zero(6 * pieceNum, 3);

            minimumJerkTrajGen(pieceNum,
                               initialPos, initialVel, initialAcc,
                               terminalPos, terminalVel, terminalAcc,
                               intermediatePositions,
                               timeAllocationVector,
                               coefficientMatrix);

            traj.clear();
            traj.reserve(pieceNum);
            for (int i = 0; i < pieceNum; i++)
            {
                traj.emplace_back(timeAllocationVector(i),
                                  coefficientMatrix.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
            }
        }

        visualizer.visualize(traj, positions.leftCols(positionNum));

        return;
    }
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "click_gen_node");
    ros::NodeHandle nh_;
    ClickGen clickGen(Config(ros::NodeHandle("~")), nh_);
    ros::spin();
    return 0;
}
