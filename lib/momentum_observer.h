/**
 * @file momentum_observer.h
 *
 * @brief Momentum observer class.
 * 
 * Expected explicit robot dynamics matrices. 
 */
#ifndef MOMENTUM_OBSERVER_H
#define MOMENTUM_OBSERVER_H

#include "external_observer.h"

#define ID_MomentumObserver 5

/**
 * @brief Momentum observer from De Luca et al. 
 */
class MomentumObserver : public ExternalObserver {
public:
  /**
   * @brief Object constructor. 
   * @param rd pointer to robot object.
   * @param k vector of joint gains.
   */
  MomentumObserver(RobotDynamics *rd, Vector& k);
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
   * @param k vector of the joint gains.
   */
  void settings(Vector& k);

private:
  // accumulators
  Vector sum, r;
  // intermediate variables
  Vector p, beta, torque, tprev;
  // coefficients
  Vector ko;
  
}; // MomentumObserver 

// Initialization
MomentumObserver::MomentumObserver(RobotDynamics *rd, Vector& k) 
  : ExternalObserver(rd,ID_MomentumObserver)
  , sum(Vector(jointNo))
  , r(Vector(jointNo))
  , p(Vector(jointNo))
  , beta(Vector(jointNo))
  , torque(Vector(jointNo))
  , tprev(Vector(jointNo))
  , ko(k)
{ 
  settings(k);
}

// Torque estimation
Vector MomentumObserver::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  p = dyn->getM(q) * qd;     // M * qd
  beta = dyn->getG(q) - dyn->getC(q,qd).transpose() * qd;  // G - C' * qd 
  
  torque = tau - dyn->getFriction(qd); // exclude friction   
  if(isRun) {
    torque += r - beta;      // tau + r - beta
    sum += 0.5 * dt * (torque + tprev);
  } else {
    torque -= beta;
    r.setZero();
    sum = p;
    isRun = true;
  }
  tprev = torque;

  p -= sum;                 // p - integral - p0

  // elementwise product
  for(int i = 0; i < jointNo; i++) {
    r(i) = ko(i) * p(i);
    torque(i) = r(i);
  }
  
  return torque;
}

// Set gains
void MomentumObserver::settings(Vector& k)
{
  ko = k;
}

#endif // MOMENTUM_OBSERVER_H
