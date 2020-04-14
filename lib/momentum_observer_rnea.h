/**
 * @file momentum_observer_rnea.h
 *
 * @brief Momentum observer class.
 * 
 * Expected robot dynamics in form of RNEA.
 */
#ifndef MOMENTUM_OBSERVER_RNEA_H
#define MOMENTUM_OBSERVER_RNEA_H

#include "external_observer.h"

#define ID_MomentumObserverRnea 25

/**
 * @brief Momentum observer from De Luca et al. 
 *
 * Find dynamics using RNEA technique.
 */
class MomentumObserverRnea : public ExternalObserverRnea {
public:
  /**
   * @brief Object constructor. 
   * @param rd pointer to robot object.
   * @param k vector of joint gains.
   */
  MomentumObserverRnea(RobotDynamicsRnea *rd, Vector& k);
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
  Vector sum, r, zero;
  // intermediate variables
  Vector p, beta, torque, tprev;
  // coefficients
  Vector ko, ko_rat;
  
}; // MomentumObserverRnea 

// Initialization
MomentumObserverRnea::MomentumObserverRnea(RobotDynamicsRnea *rd, Vector& k) 
  : ExternalObserverRnea(rd,ID_MomentumObserverRnea)
  , sum(Vector(jointNo))
  , r(Vector(jointNo))
  , zero(Vector::Zero(jointNo))
  , p(Vector(jointNo))
  , beta(Vector(jointNo))
  , torque(Vector(jointNo))
  , tprev(Vector(jointNo))
  , ko(k)
  , ko_rat(Vector(jointNo))
{ 
  settings(k);
}

// Torque estimation
Vector MomentumObserverRnea::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  p = dyn->rnea(q,zero,qd);     // M * qd
  beta = dyn->rnea(q,zero,zero,GRAVITY)- dyn->tranCqd(q,qd);  // G - C' * qd 
  
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
    torque(i) = ko_rat(i) * r(i);   // reuse variable
  }
  
  return torque;
}

// Set gains
void MomentumObserverRnea::settings(Vector& k)
{
  ko = k;
  for(int i = 0; i < jointNo; i++) ko_rat(i) = k(i)/(1+k(i));
}

#endif // MOMENTUM_OBSERVER_RNEA_H
