// Copyright 2020-2024 Stanislav Mikhel

/**
 * @file sliding_mode_observer.h
 *
 * @brief Sliding mode observer class.
 *
 * Expected explicit robot dynamics matrices.
 */
#ifndef SLIDING_MODE_OBSERVER_H
#define SLIDING_MODE_OBSERVER_H

#include "external_observer.h"

#define ID_SlidingModeObserver 3

/**
 * @brief Sliding mode observer from Garofalo et. al.
 *
 * "Sliding mode momentum observers for estimation of external torques
 *  and joint acceleration", 2019
 */
class SlidingModeObserver final : public ExternalObserver {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to the robot object.
   * @param t1 ...
   * @param s1 ...
   * @param t2 ...
   * @param s2 ...
   */
  SlidingModeObserver(RobotDynamics* rd, VectorJ& t1, VectorJ& s1, VectorJ& t2, VectorJ& s2);

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
  VectorJ p, spp, dp_hat, torque, dsigma;
  VectorJ T1, S1, T2, S2;
};  // SlidingModeObserver


// Initialization
SlidingModeObserver::SlidingModeObserver(RobotDynamics* rd,
                                VectorJ& t1, VectorJ& s1, VectorJ& t2, VectorJ& s2)
  : ExternalObserver(rd, ID_SlidingModeObserver)
  , sigma(VectorJ(jointNo))
  , p_hat(VectorJ(jointNo))
  , p(VectorJ(jointNo))
  , spp(VectorJ(jointNo))
  , dp_hat(VectorJ(jointNo))
  , torque(VectorJ(jointNo))
  , dsigma(VectorJ(jointNo))
  , T1(VectorJ(jointNo)), S1(VectorJ(jointNo))
  , T2(VectorJ(jointNo)), S2(VectorJ(jointNo))
{
  settings(t1, s1, t2, s2);
}

// Update parameters
void SlidingModeObserver::settings(VectorJ& t1, VectorJ& s1, VectorJ& t2, VectorJ& s2)
{
  T1 = t1;
  S1 = s1;
  T2 = t2;
  S2 = s2;
}

// External torque
VectorJ SlidingModeObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;
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
  // for(int i = 0; i < jointNo; i++) spp(i) = (p(i) > 0 ? 1 : (p(i) < 0 ? -1 : 0));

  dp_hat = torque + dyn->getC(q, qd).transpose()*qd - dyn->getG(q) + sigma;
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

#endif  // SLIDING_MODE_OBSERVER_H
