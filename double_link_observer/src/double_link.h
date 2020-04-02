#ifndef DOUBLE_LINK_H
#define DOUBLE_LINK_H

#include "external_observer.h" 
#include <cmath>

class DoubleLink : public RobotDynamics {
public:
  DoubleLink();
  
  // get inertia
  Matrix getM(Vector& q);
  // get Coriolis/centrifugal matrix
  Matrix getC(Vector& q, Vector& qd);
  // get GRAVITYity
  Vector getG(Vector& q);
  // friction model 
  Vector getFriction(Vector& qd) { return fric; }
  // number of joints
  int jointNo() { return 2; }
      
private:
  Matrix M, C, J;
  Vector G, fric;
  double m1,m2;
  double l1, l2;
  double lc1, lc2; 
  double I1, I2;
};

DoubleLink::DoubleLink() 
           : RobotDynamics()
           , M(Matrix(2,2))
           , C(Matrix(2,2))
           , J(Matrix(6,2))
           , G(Vector(2))    
           , fric(Vector(2)) 
{
  m1 = 1; m2 = 1;
  l1 = 0.5; l2 = 0.5;
  lc1 = 0.25; lc2 = 0.25;
  I1 = 0.3; I2 = 0.2;
  M.setZero(); 
  C.setZero(); 
  J.setZero(); 
  G.setZero(); 
  fric.setZero();
}

Matrix DoubleLink::getM(Vector& q)
{
  M(0,0) = m1*lc1*lc1 + m2*(l1*l1 + lc2*lc2 + 2*l1*lc2*cos(q(1))) + I1 + I2;
  M(0,1) = m2*(lc2*lc2 + l1*lc2*cos(q(1))) + I2;
  M(1,0) = M(0,1);
  M(1,1) = m2*lc2*lc2 + I2; 
  return M;
}

Matrix DoubleLink::getC(Vector& q, Vector& qd)
{
  double h = -m2*l1*lc2*sin(q(1)); 
  C(0,0) = h*qd(1);
  C(0,1) = h*(qd(0)+qd(1));
  C(1,0) = -h*(qd(0)); 
  return C;
}

Vector DoubleLink::getG(Vector& q)
{
  double c12 = cos(q(0)+q(1));
  G(0) = (m1*lc1 + m2*l1)*GRAVITY*cos(q(0)) + m2*lc2*GRAVITY*c12;
  G(1) = m2*lc2*GRAVITY*c12; 
  return G;
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

#endif // DOUBLE_LINK_H
