
#include <iostream> 
#include <fstream>
#include "double_link_rnea.h"

//
// Uncomment the desirable observer
//
//#include "../lib/momentum_observer_rnea.h" 
//#include "../lib/disturbance_observer_rnea.h"
#include "../lib/sliding_mode_observer_rnea.h"

#define OMEGA1 1.3
#define OMEGA2 0.8 
#define FNAME "force.csv"

// use it to emulate external torque
#define SET_TORQUE

int main(int argc, char** argv)
{
  
  DoubleLinkRnea robot; 
  Vector q(2), qd(2), q2d(2), tau(2), ext(2);

#ifdef MOMENTUM_OBSERVER_RNEA_H
  Vector k(2);
  k << 50,50; 
  MomentumObserverRnea m_observer(&robot,k);
#endif 
#ifdef DISTURBANCE_OBSERVER_RNEA_H
  double sigma = 21, xeta = 18, beta = 50;
  DisturbanceObserverRnea d_observer(&robot,sigma,xeta,beta); 
#endif
#ifdef SLIDING_MODE_OBSERVER_RNEA_H
  Vector T1(2), S1(2), T2(2), S2(2);
  S1 << 20,30; 
  T1(0) = 2*sqrt(S1(0)); T1(1) = 2*sqrt(S1(1));
  S2 << 10,10; // set 0 to exclude linear part
  T2(0) = 2*sqrt(S2(0)); T2(1) = 2*sqrt(S2(1));
  SlidingModeObserverRnea sm_observer(&robot,T1,S1,T2,S2);
#endif
  
  // save to file 
  std::ofstream file;
  file.open(FNAME); 

  double dt = 0.01;
  for(double t = 0; t < 3; t += dt) {
    double c1 = cos(OMEGA1*t), c2 = cos(OMEGA2*t);
    double s1 = sin(OMEGA1*t), s2 = sin(OMEGA2*t); 
    // state
    q(0) = s1; q(1) = s2;
    qd(0) = OMEGA1*c1; qd(1) = OMEGA2*c2;
    q2d(0) = -OMEGA1*OMEGA1*s1; q2d(1) = -OMEGA2*OMEGA2*s2;
    
    // torque 
    tau = robot.rnea(q,qd,q2d,GRAVITY); 
#ifdef SET_TORQUE    
    if(t > 1 && t < 2) {
      tau(0) -= 0.5;
      tau(1) -= 0.5;
    }
#endif // SET_TORQUE
    // estimate torque and save 
    
#ifdef MOMENTUM_OBSERVER_RNEA_H
    ext = m_observer.getExternalTorque(q,qd,tau,dt);
#endif 
#ifdef DISTURBANCE_OBSERVER_RNEA_H
    ext = d_observer.getExternalTorque(q,qd,tau,dt);
#endif
#ifdef SLIDING_MODE_OBSERVER_RNEA_H
    ext = sm_observer.getExternalTorque(q,qd,tau,dt);
#endif 
    
    file << t << "," << ext(0) << "," << ext(1) << std::endl;
  }
      
  file.close();
  
  return 0;
}

