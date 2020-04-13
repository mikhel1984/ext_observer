// Generate "external" torques for the double link manipulator.
// Use
//   rqt_plot /ext/torque/tau1:tau2:tau3 etc.
// for visualization

#include <string>

#include "ros/ros.h"
#include "ur_observer/ExtTorque6.h"  // messages

#include "include/momentum_observer.h" 
#include "include/disturbance_observer.h"
#include "include/sliding_mode_observer.h"

#include "include/ur_robot.h"
#include "include/ur_listener.h"

#define MO_OBSERVER  "MO"
#define DIS_OBSERVER "DIS"
#define SM_OBSERVER  "SM"
#define UR_PORT      30003


int main(int argc, char** argv)
{
  // ROS settings
  ros::init(argc, argv, "dl_observer");
  ros::NodeHandle n;
  ros::Publisher torque_pub = n.advertise<ur_observer::ExtTorque6>("/ext/torque",1);
  ur_observer::ExtTorque6 msg;
  
  // Robot
  UrRobot robot; 
  double K_i[UR_JOINTS] = {10.0,10.6956,8.4566,9.0029,9.48,10.1232};

  // Read parameters
  std::string s, addr;
  if(!n.getParam("ur_type", s))  s = MO_OBSERVER;
  if(!n.getParam("ur_ip", addr)) addr = "127.0.0.1";

  // Connection
  UrListener connection(addr.c_str(), UR_PORT);
  if(!connection.isConnected()) ROS_ERROR("Connection failed");

  // Choose observer
  MomentumObserver    *m_observer = 0;
  DisturbanceObserver *d_observer = 0;
  SlidingModeObserver *sm_observer = 0;
  
  ExternalObserver    *observer = 0;

  if(s == MO_OBSERVER) {
    Vector k(UR_JOINTS);
    k << 90,50,50,90,90,40;
    m_observer = new MomentumObserver(&robot,k);
    observer = m_observer;
  } else if(s == DIS_OBSERVER) {
    double sigma = 21, xeta = 18, beta = 50;
    d_observer = new DisturbanceObserver(&robot,sigma,xeta,beta);
    observer = d_observer;
  } else if(s == SM_OBSERVER) {
    Vector T1(UR_JOINTS), S1(UR_JOINTS), T2(UR_JOINTS), S2(UR_JOINTS);
    S1 << 20,30,20,30,20,30;
    for(int i = 0; i < UR_JOINTS; i++) T1(i) = 2*sqrt(S1(i));
    S2 << 10,10,10,10,10,10;
    for(int i = 0; i < UR_JOINTS; i++) T2(i) = 2*sqrt(S2(i));
    sm_observer = new SlidingModeObserver(&robot,T1,S1,T2,S2);
    observer = sm_observer;
  } else {
    ROS_ERROR("Unknown parameter");
  }

  ROS_INFO("Observer type: %s", s.c_str());

  // Trajectory parameters
  Vector q(UR_JOINTS), qd(UR_JOINTS), ii(UR_JOINTS), ext(UR_JOINTS);
  double tprev = -1, tm = 0;

  ros::Duration(0.5).sleep();  // wait a little to establish connection

  while(ros::ok()) {
    
    if(connection.listen(tm,q,qd,ii)) {
      // find torques 
      for(int j = 0; j < UR_JOINTS; j++) ii(j) *= K_i[j];
      if(tprev < 0) tprev = tm;
      
      // find external torques
      ext = observer->getExternalTorque(q,qd,ii,tm-tprev);
      tprev = tm;
      // publish 
      msg.time = tm;
      msg.tau1 = ext(0);
      msg.tau2 = ext(1);
      msg.tau3 = ext(2);
      msg.tau4 = ext(3);
      msg.tau5 = ext(4);
      msg.tau6 = ext(5);
      torque_pub.publish(msg);
    }
  }

  if(m_observer)  delete m_observer;
  if(d_observer)  delete d_observer;
  if(sm_observer) delete sm_observer;

  return 0;
}
