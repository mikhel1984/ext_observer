// Test observers with UR10 expreimental results

#include <iostream> 
#include "ur_robot.h"
#include "ur_listener.h"

//
//   Uncomment desired observer
//
#include "../src/momentum_observer.h" 
//#include "../src/disturbance_observer.h"
//#include "../src/sliding_mode_observer.h"

// settings
//#define OUT_NAME "external.csv"

// torque coefficients
double K_i[] = {10.0,10.6956,8.4566,9.0029,9.48,10.1232};

// access to objects
ExternalObserver* observer;

void printExt(double t, Vector& q, Vector& qd, Vector& ii);

int main(int argc, char** argv)
{  

  UrRobot robot; 
  int N = robot.jointNo();

#ifdef MOMENTUM_OBSERVER_H
  Vector k(N);
  k << 90,50,50,90,90,40; 
  MomentumObserver m_observer(&robot,k);
  observer = &m_observer;
#endif 
#ifdef DISTURBANCE_OBSERVER_H
  double sigma = 21, xeta = 18, beta = 50;
  DisturbanceObserver d_observer(&robot,sigma,xeta,beta); 
  observer = &d_observer;
#endif
#ifdef SLIDING_MODE_OBSERVER_H
  Vector T1(N), S1(N), T2(N), S2(N);
  S1 << 20,30,20,30,20,30; 
  for(int i = 0; i < N; i++) T1(i) = 2*sqrt(S1(i));
  S2 << 10,10,10,10,10,10; // set 0 to exclude linear part
  for(int i = 0; i < N; i++) T2(i) = 2*sqrt(S2(i));
  SlidingModeObserver sm_observer(&robot,T1,S1,T2,S2);
  observer = &d_observer;
#endif
  
  // save to file 
  //std::ofstream oFile(OUT_NAME);
  //oFile.setf(std::ios::fixed);
  //oFile.precision(4);

  // TCP connection
  UrListener connection("127.0.0.1",30003);
  if(!connection.isConnected()) {
    std::cout << "Connection error" << std::endl;
    return 1;
  }

  connection.listen(printExt);
      
  //oFile.close();
  
  return 0;
}

// data processing
void printExt(double t, Vector& q, Vector& qd, Vector& ii)
{
  static double tprev = -1;
  double dt = 0;

  // time step
  if(tprev >= 0) dt = t - tprev;
  tprev = t;
  
  // find torque
  for(int j = 0; j < UR_JOINTS; j++) ii(j) *= K_i[j];

  // evaluate and show
  std::cout << t << " ";
  std::cout << observer->getExternalTorque(q,qd,ii,dt).transpose() << std::endl;
}
