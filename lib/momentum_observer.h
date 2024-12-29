// Copyright 2020-2024 Stanislav Mikhel

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
class MomentumObserver final : public ExternalObserver {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to robot object.
   * @param k vector of joint gains.
   */
  MomentumObserver(RobotDynamics *rd, VectorJ& k);

  /**
   * @brief External torque estimation.
   * @param q joint angle vector.
   * @param qd joint velocity vector.
   * @param tau joint torque vector.
   * @param dt time step.
   * @return external torque vector.
   */
  VectorJ getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt) override;

  /**
   * @brief Observer settings.
   * @param k vector of the joint gains.
   */
  inline void settings(VectorJ& k) { ko = k; }

private:
  // accumulators
  VectorJ sum, r;
  // intermediate variables
  VectorJ p, beta, torque, tprev;
  // coefficients
  VectorJ ko;
};  // MomentumObserver


// Initialization
MomentumObserver::MomentumObserver(RobotDynamics *rd, VectorJ& k)
  : ExternalObserver(rd, ID_MomentumObserver)
  , sum(VectorJ(jointNo))
  , r(VectorJ(jointNo))
  , p(VectorJ(jointNo))
  , beta(VectorJ(jointNo))
  , torque(VectorJ(jointNo))
  , tprev(VectorJ(jointNo))
  , ko(k)
{
  settings(k);
}

// Torque estimation
VectorJ MomentumObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;     // M * qd
  beta = dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;  // G - C' * qd

  torque = tau - dyn->getFriction(qd);  // exclude friction
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

#endif  // MOMENTUM_OBSERVER_H
