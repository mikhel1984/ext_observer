/**
 * @file disturbance_kalman_filter.h
 *
 * @brief Disturbance kalman filter observer.
 *
 * Expected explicit robot dynamics matrices.
 */
#ifndef DISTURBANCE_KALMAN_FILTER_EXP_H
#define DISTURBANCE_KALMAN_FILTER_EXP_H

#include "external_observer.h"
#include "kalman_filter_continous.h"

#define ID_DKalmanObserverExp 7

//#include <iostream>

/**
 * @brief Disturbance Kalman filter from Hu et. al.
 *
 * Use matrix exponent for discretizaion.
 */
class DKalmanObserverExp : public ExternalObserver {
public:
  /**
   * @brief Object constructor
   * @param rd pointer to the robot object.
   * @param s disturbance dynamics matrix.
   * @param h disturbance observation matrix.
   * @param q covariance of the process noise.
   * @param r covariance of the observation noise.
   */
  DKalmanObserverExp(RobotDynamics* rd, Matrix& s, Matrix& h, Matrix& q, Matrix& r);
  /**
   * @brief Object destructor.
   */
  ~DKalmanObserverExp();
  /**
   * @brief External torque estimation.
   * @param q joint angle vector.
   * @param qd joint velocity vector.
   * @param tau joint torque vector.
   * @param dt time step.
   * @return external torque vector.
   */
  Vector getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt);
  /**
   * @brief Observer settings.
   * @param q covariance of the process noise.
   * @param r covariance of the observation noise.
   */
  void settings(Matrix& q, Matrix& r); 

private:
  // disturbance observation
  Matrix H, M;
  // state, input and momentum
  Vector X, u, p;
  // Kalman filter object
  KalmanFilterContinous *filter;

}; // DKalmanObserverExp

// Initialization
DKalmanObserverExp::DKalmanObserverExp(RobotDynamics* rd, Matrix& s, Matrix& h, Matrix& q, Matrix& r)
  : ExternalObserver(rd,ID_DKalmanObserverExp)
  , H(h)
  , M(Matrix(jointNo,jointNo))
  , X(Vector(2*jointNo))
  , u(Vector(jointNo))
  , p(Vector(jointNo))
  , filter(0)  
{
  // prepare matrices
  Matrix A(2*jointNo,2*jointNo), B(2*jointNo,jointNo), C(jointNo,2*jointNo); 
  B.setZero();
  C.setZero();
  // B & C matrices
  for(int i = 0; i < jointNo; i++) {
    B(i,i) = 1;
    C(i,i) = 1;
  }
  // A matrix
  A.setZero(); 
  A.block(      0,jointNo,jointNo,jointNo) = h;
  A.block(jointNo,jointNo,jointNo,jointNo) = s;
  // make filter
  filter = new KalmanFilterContinous(A,B,C);
  filter->setCovariance(q,r);
}

// Clear dynamically allocated memory
DKalmanObserverExp::~DKalmanObserverExp()
{
  delete filter;
}

// Torque estimation
Vector DKalmanObserverExp::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  M = dyn->getM(q);
  p = M * qd;
  u = tau - dyn->getG(q) - dyn->getFriction(qd);
  u += dyn->getC(q,qd).transpose() * qd;

  filter->updateR(M);
  
  if(isRun) {
    X = filter->step(u,p,dt);
  } else {
    // prepare X0
    X.setZero();
    for(int i = 0; i < jointNo; i++) X(i) = p(i);
    // reset
    filter->reset(X);
    isRun = true;
  }
  // disturbance = H * omega   (reuse variable)
  p = H * X.block(jointNo,0,jointNo,1);

  return p;
}

// settings
void DKalmanObserverExp::settings(Matrix& Q, Matrix& R)
{
  filter->setCovariance(Q,R);
}

#endif // DISTURBANCE_KALMAN_FILTER_EXP_H
