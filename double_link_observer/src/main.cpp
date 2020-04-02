#include "ros/ros.h"
// Generate "external" torques for the double link manipulator.
// Use
//   rqt_plot /ext/torque/tau1:tau2 
// for visualization

#include "double_link_observer/ExtTorque2.h"

#include "double_link.h"
#include "momentum_observer.h"

#define OMEGA1 1.3
#define OMEGA2 0.8 

int main(int argc, char** argv)
{
  // ROS settings
  ros::init(argc, argv, "dl_observer");
  ros::NodeHandle n;

  ros::Publisher torque_pub = n.advertise<double_link_observer::ExtTorque2>("/ext/torque",1);
  double_link_observer::ExtTorque2 msg;
  // Robot and observer settings
  DoubleLink robot;
  Vector k(2);
  k << 50,50;
  MomentumObserver m_observer(&robot, k);

  // Trajectory parameters
  Vector q(2), qd(2), q2d(2), tau(2), ext(2);
  bool addTorque = false;

  ros::Duration(1).sleep();  // wait a little to establish connection
  double tbegin = ros::Time::now().toSec();
  double point = 0, tprev = 0;

  ROS_INFO("double_link_observer started");
  ros::Rate rate(100);
  while(ros::ok()) {
    double t = ros::Time::now().toSec() - tbegin;
    double c1 = cos(OMEGA1*t), c2 = cos(OMEGA2*t);
    double s1 = sin(OMEGA1*t), s2 = sin(OMEGA2*t); 
    // state
    q(0) = s1; q(1) = s2;
    qd(0) = OMEGA1*c1; qd(1) = OMEGA2*c2;
    q2d(0) = -OMEGA1*OMEGA1*s1; q2d(1) = -OMEGA2*OMEGA2*s2;
    tau = robot.getM(q)*q2d + robot.getC(q,qd)*qd + robot.getG(q); 
    // "external" torque
    if(addTorque) {
      tau(0) += 0.5;
      tau(1) += 0.4;
    }
    if((t-point) >= 2) {
      addTorque = !addTorque;
      point = t;
    }
    // evaluate
    ext = m_observer.getExternalTorque(q,qd,tau,t-tprev);
    tprev = t;
    // publish 
    msg.time = t;
    msg.tau1 = ext(0);
    msg.tau2 = ext(1);
    torque_pub.publish(msg);

    rate.sleep();
  }

  return 0;
}
