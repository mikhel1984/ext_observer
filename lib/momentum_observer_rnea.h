// Copyright 2020-2024 Stanislav Mikhel

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
class MomentumObserverRnea final : public ExternalObserverRnea {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to robot object.
   * @param k vector of joint gains.
   */
  MomentumObserverRnea(RobotDynamicsRnea *rd, VectorJ& k);

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
  VectorJ sum, r, zero;
  // intermediate variables
  VectorJ p, beta, torque, tprev;
  // coefficients
  VectorJ ko;
};  // MomentumObserverRnea


// Initialization
MomentumObserverRnea::MomentumObserverRnea(RobotDynamicsRnea *rd, VectorJ& k)
  : ExternalObserverRnea(rd, ID_MomentumObserverRnea)
  , sum(VectorJ(jointNo))
  , r(VectorJ(jointNo))
  , zero(VectorJ::Zero(jointNo))
  , p(VectorJ(jointNo))
  , beta(VectorJ(jointNo))
  , torque(VectorJ(jointNo))
  , tprev(VectorJ(jointNo))
  , ko(k)
{
  settings(k);
}

// Torque estimation
VectorJ MomentumObserverRnea::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->rnea(q, zero, qd);     // M * qd
  beta = dyn->rnea(q, zero, zero, GRAVITY)- dyn->tranCqd(q, qd);  // G - C' * qd

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
    torque(i) = r(i);   // reuse variable
  }

  return torque;
}

#endif  // MOMENTUM_OBSERVER_RNEA_H
