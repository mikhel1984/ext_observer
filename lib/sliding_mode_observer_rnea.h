// Copyright 2020-2024 Stanislav Mikhel

/**
 * @file sliding_mode_observer_rnea.h
 *
 * @brief Sliding mode observer class.
 *
 * Expected robot dynamics in form of RNEA.
 */
#ifndef SLIDING_MODE_OBSERVER_RNEA_H
#define SLIDING_MODE_OBSERVER_RNEA_H

#include "external_observer.h"

#define ID_SlidingModeObserverRnea 23

/**
 * @brief Sliding mode observer from ...
 *
 * Use RNEA for dynamics.
 */
class SlidingModeObserverRnea final : public ExternalObserverRnea {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to the robot object.
   * @param t1 ...
   * @param s1 ...
   * @param t2 ...
   * @param s2 ...
   */
  SlidingModeObserverRnea(
    RobotDynamicsRnea* rd, VectorJ& t1, VectorJ& s1, VectorJ& t2, VectorJ& s2);

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
   * @param t1 ...
   * @param s1 ...
   * @param t2 ...
   * @param s2 ...
   */
  void settings(VectorJ& t1, VectorJ& s1, VectorJ& t2, VectorJ& s2);

private:
  const double BIG = 50.0;  /**< Map tanh to sign. */

  // Temporary objects
  VectorJ sigma, p_hat;
  VectorJ p, spp, dp_hat, torque, dsigma, zeros;
  VectorJ T1, S1, T2, S2;
};  // SlidingModeObserverRnea


// Initialization
SlidingModeObserverRnea::SlidingModeObserverRnea(
  RobotDynamicsRnea* rd, VectorJ& t1, VectorJ& s1, VectorJ& t2, VectorJ& s2)
  : ExternalObserverRnea(rd, ID_SlidingModeObserverRnea)
  , sigma(VectorJ(jointNo))
  , p_hat(VectorJ(jointNo))
  , p(VectorJ(jointNo))
  , spp(VectorJ(jointNo))
  , dp_hat(VectorJ(jointNo))
  , torque(VectorJ(jointNo))
  , dsigma(VectorJ(jointNo))
  , zeros(VectorJ::Zero(jointNo))
  , T1(VectorJ(jointNo)), S1(VectorJ(jointNo))
  , T2(VectorJ(jointNo)), S2(VectorJ(jointNo))
{
  settings(t1, s1, t2, s2);
}

// Update settings
void SlidingModeObserverRnea::settings(VectorJ& t1, VectorJ& s1, VectorJ& t2, VectorJ& s2)
{
  T1 = t1;
  S1 = s1;
  T2 = t2;
  S2 = s2;
}

// External torque
VectorJ SlidingModeObserverRnea::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->rnea(q, zeros, qd);
  torque = tau - dyn->getFriction(qd);

  if(!isRun) {
    sigma.setZero();
    p_hat = p;
    isRun = true;
  }

  p = p_hat - p;  // reuse p
  for(int i = 0; i < jointNo; i++) {
    spp(i) = tanh(p(i)*BIG);
  }

  dp_hat = torque + dyn->tranCqd(q, qd) - dyn->rnea(q, zeros, zeros, GRAVITY) + sigma;
  for(int i = 0; i < jointNo; i++) {
    // - T2*p
    dp_hat(i) -= T2(i)*p(i);
    // - sqrt(abs(p)).*T1*spp
    dp_hat(i) -= sqrt( fabs(p(i)) ) * T1(i) * spp(i);
  }
  // dsigma
  for(int i = 0; i < jointNo; i++) {
    dsigma(i) = -S1(i)*spp(i) - S2(i)*p(i);
  }

  p = sigma;     // reuse to save result
  p_hat += dt*dp_hat;
  sigma += dt*dsigma;

  return p;
}

#endif  // SLIDING_MODE_OBSERVER_RNEA_H
