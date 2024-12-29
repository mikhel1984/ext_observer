// Copyright 2020-2024 Stanislav Mikhel

/**
 * @file disturbance_observer.h
 *
 * @brief Disturbance observer class.
 *
 * Expected explicit robot dynamics matrices.
 */
#ifndef DISTURBANCE_OBSERVER_H
#define DISTURBANCE_OBSERVER_H

#include "external_observer.h"

#define ID_DisturbanceObserver 4

/**
 * @brief Disturbance observer from Mohammadi et. al.
 */
class DisturbanceObserver final : public ExternalObserver {
public:
  /**
   * @brief Object constructor.
   * @param rd pointer to the robot object.
   * @param sigma ...
   * @param xeta ...
   * @param beta ...
   */
  DisturbanceObserver(RobotDynamics *rd, double sigma, double xeta, double beta);

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
  VectorJ p, z, torque;
};  // DisturbanceObserver


// Initialization
DisturbanceObserver::DisturbanceObserver(RobotDynamics *rd, double sigma, double xeta, double beta)
  : ExternalObserver(rd, ID_DisturbanceObserver)
  , Y(MatrixJ(jointNo, jointNo))
  , L(MatrixJ(jointNo, jointNo))
  , I(MatrixJ::Identity(jointNo, jointNo))
  , lft(MatrixJ(jointNo, jointNo))
  , rht(MatrixJ(jointNo, jointNo))
  , p(VectorJ(jointNo))
  , z(VectorJ(jointNo))
  , torque(VectorJ(jointNo))
{
  settings(sigma, xeta, beta);
}

// Get torque
VectorJ DisturbanceObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  L = Y * dyn->getM(q).inverse();
  L *= dt;
  p = Y * qd;

  if(isRun) {
    torque = tau - dyn->getFriction(qd);
    lft = I + L;
    rht = z + L*(dyn->getC(q, qd)*qd + dyn->getG(q) - torque - p);
    z = lft.inverse() * rht;
  } else {
    z = -p;
    isRun = true;
  }

  p += z;

  return p;
}

// Update parameters
void DisturbanceObserver::settings(double sigma, double xeta, double beta)
{
  double k = 0.5*(xeta + 2*beta*sigma);
  Y = k * MatrixJ::Identity(jointNo, jointNo);
}

#endif  // DISTURBANCE_OBSERVER_H
