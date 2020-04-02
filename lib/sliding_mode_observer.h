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

#define BIG 50   /**< Map tanh to sign. */

/**
 * @brief Sliding mode observer from ...
 */
class SlidingModeObserver : public ExternalObserver {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to the robot object.
   * @param t1 ...
   * @param s1 ...
   * @param t2 ...
   * @param s2 ...
   */
  SlidingModeObserver(RobotDynamics* rd, Vector& t1, Vector& s1, Vector& t2, Vector& s2);
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
   * @param t1 ...
   * @param s1 ...
   * @param t2 ...
   * @param s2 ...
   */
  void settings(Vector& t1, Vector& s1, Vector& t2, Vector& s2);

private:
 // Temporary objects
 Vector sigma, p_hat;
 Vector p, spp, dp_hat, torque, dsigma;
 Vector T1, S1, T2, S2;

}; // SlidingModeObserver

// Initialization
SlidingModeObserver::SlidingModeObserver(RobotDynamics* rd,
                                Vector& t1, Vector& s1, Vector& t2, Vector& s2)
  : ExternalObserver(rd)
  , sigma(Vector(jointNo))
  , p_hat(Vector(jointNo))
  , p(Vector(jointNo))
  , spp(Vector(jointNo))
  , dp_hat(Vector(jointNo))
  , torque(Vector(jointNo))
  , dsigma(Vector(jointNo))
  , T1(Vector(jointNo)), S1(Vector(jointNo))
  , T2(Vector(jointNo)), S2(Vector(jointNo))
{
  settings(t1,s1,t2,s2);
}

// Update parameters
void SlidingModeObserver::settings(Vector& t1, Vector& s1, Vector& t2, Vector& s2)
{
  T1 = t1;
  S1 = s1;
  T2 = t2;
  S2 = s2;
}

// External torque
Vector SlidingModeObserver::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  p = dyn->getM(q) * qd;  
  torque = tau - dyn->getFriction(qd);

  if(!isRun) {
    sigma.setZero();
    p_hat = p;
    isRun = true;
  }

  p = p_hat - p;  // reuse p
  for(int i = 0; i < jointNo; i++) spp(i) = tanh(p(i)*BIG);

  dp_hat = torque + dyn->getC(q,qd).transpose()*qd - dyn->getG(q) + sigma;
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

#endif // SLIDING_MODE_OBSERVER_H
