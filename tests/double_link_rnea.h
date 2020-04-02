#ifndef DOUBLE_LINK_RNEA_H
#define DOUBLE_LINK_RNEA_H

#include "../lib/external_observer.h" 
#include <cmath>

class DoubleLinkRnea : public RobotDynamicsRnea {
public:
  DoubleLinkRnea();
  
  Vector rnea(Vector& q, Vector& qd, Vector& q2d, double g = 0);
  // friction model 
  Vector getFriction(Vector& qd) { return fric; }
  // number of joints
  int jointNo() { return 2; }
      
private:
  Vector tau, fric;
  double m1,m2;
  double l1, l2;
  double lc1, lc2; 
  double I1, I2;
}; // DoubleLinkRnea

DoubleLinkRnea::DoubleLinkRnea() 
           : RobotDynamicsRnea()
           , tau(Vector(2))    
           , fric(Vector(2)) 
{
  m1 = 1; m2 = 1;
  l1 = 0.5; l2 = 0.5;
  lc1 = 0.25; lc2 = 0.25;
  I1 = 0.3; I2 = 0.2;
  fric.setZero();
}

Vector DoubleLinkRnea::rnea(Vector& q, Vector& qd, Vector& q2d, double g)
{
  double M00 = m1*lc1*lc1 + m2*(l1*l1 + lc2*lc2 + 2*l1*lc2*cos(q(1))) + I1 + I2;
  double M01 = m2*(lc2*lc2 + l1*lc2*cos(q(1))) + I2;
  double M10 = M01;
  double M11 = m2*lc2*lc2 + I2;
  double h = -m2*l1*lc2*sin(q(1)); 
  double C00 = h*qd(1);
  double C01 = h*(qd(0)+qd(1));
  double C10 = -h*(qd(0)); 
  double c12 = cos(q(0)+q(1));
  double G0 = (m1*lc1 + m2*l1)*g*cos(q(0)) + m2*lc2*g*c12;
  double G1 = m2*lc2*g*c12; 
  
  tau(0) = M00*q2d(0) + M01*q2d(1) + C00*qd(0) + C01*qd(1) + G0; 
  tau(1) = M10*q2d(0) + M11*q2d(1) + C10*qd(0)             + G1;

  return tau;
}

/*
Matrix DoubleLink::getJacobian(Vector& q)
{
  double s12 = sin(q(0)+q(1)), c12 = cos(q(0)+q(1));
  J(0,0) = -l1*sin(q(0))-l2*s12;
  J(1,0) = l1*cos(q(0)) + l2*c12;
  J(5,0) = 1;
  J(0,1) = -l2*s12;
  J(1,1) = l2*c12;
  J(5,1) = 1; 
  return J;
}
*/

#endif // DOUBLE_LINK_RNEA_H
