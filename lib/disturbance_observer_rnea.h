// Copyright 2020-2024 Stanislav Mikhel

/**
 * @file disturbance_observer_rnea.h
 *
 * @brief Disturbance observer class.
 *
 * Expected robot dynamics in form of RNEA.
 */
#ifndef DISTURBANCE_OBSERVER_RNEA_H
#define DISTURBANCE_OBSERVER_RNEA_H

#include "external_observer.h"

#define ID_DisturbanceObserverRnea 24

/**
 * @brief Disturbance observer from Mohammadi et. al.
 *
 * Use RNEA technique for dynamics.
 */
class DisturbanceObserverRnea final : public ExternalObserverRnea {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to the robot object.
   * @param sigma ...
   * @param xeta ...
   * @param beta ...
   */
  DisturbanceObserverRnea(RobotDynamicsRnea *rd, double sigma, double xeta, double beta);

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
   * @param sigma ...
   * @param xeta ...
   * @param beta ...
   */
  void settings(double sigma, double xeta, double beta);

private:
  // temporary objects
  MatrixJ Y, L, I;
  MatrixJ lft, rht;
  VectorJ p, z, torque, zeros;
};  // DisturbanceObserverRnea


// Initialization
DisturbanceObserverRnea::DisturbanceObserverRnea(
  RobotDynamicsRnea *rd, double sigma, double xeta, double beta)
  : ExternalObserverRnea(rd, ID_DisturbanceObserverRnea)
  , Y(MatrixJ(jointNo, jointNo))
  , L(MatrixJ(jointNo, jointNo))
  , I(MatrixJ::Identity(jointNo, jointNo))
  , lft(MatrixJ(jointNo, jointNo))
  , rht(MatrixJ(jointNo, jointNo))
  , p(VectorJ(jointNo))
  , z(VectorJ(jointNo))
  , torque(VectorJ(jointNo))
  , zeros(VectorJ::Zero(jointNo))
{
  settings(sigma, xeta, beta);
}

// Get torque
VectorJ DisturbanceObserverRnea::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  L = Y * dyn->getM(q).inverse();
  L *= dt;
  p = Y * qd;

  if(isRun) {
    torque = tau - dyn->getFriction(qd);
    lft = I + L;
    rht = z + L*(dyn->rnea(q, qd, zeros, GRAVITY) - torque - p);
    z = lft.inverse() * rht;
  } else {
    z = -p;
    isRun = true;
  }

  p += z;

  return p;
}

// Update parameters
void DisturbanceObserverRnea::settings(double sigma, double xeta, double beta)
{
  double k = 0.5*(xeta + 2*beta*sigma);
  Y = k * MatrixJ::Identity(jointNo, jointNo);
}

#endif  // DISTURBANCE_OBSERVER_RNEA_H
