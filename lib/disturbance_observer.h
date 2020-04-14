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
class DisturbanceObserver : public ExternalObserver {
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
  Vector getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt);
  /**
   * @brief Observer settings.
   * @param sigma ...
   * @param xeta ...
   * @param beta ...
   */
  void settings(double sigma, double xeta, double beta);

private:
  // temporary objects
  Matrix Y, L, I;
  Matrix lft, rht;
  Vector p, z, torque;  
  
}; // DisturbanceObserver 

// Initialization
DisturbanceObserver::DisturbanceObserver(RobotDynamics *rd, double sigma, double xeta, double beta)
  : ExternalObserver(rd,ID_DisturbanceObserver)
  , Y(Matrix(jointNo,jointNo))
  , L(Matrix(jointNo,jointNo))
  , I(Matrix::Identity(jointNo,jointNo))
  , lft(Matrix(jointNo,jointNo))
  , rht(Matrix(jointNo,jointNo))
  , p(Vector(jointNo))
  , z(Vector(jointNo))
  , torque(Vector(jointNo))
{
  settings(sigma,xeta,beta); 
}

// Get torque
Vector DisturbanceObserver::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  L = Y * dyn->getM(q).inverse();
  L *= dt; 
  p = Y * qd;
  
  if(isRun) {
    torque = tau - dyn->getFriction(qd);
    lft = I + L;
    rht = z + L*(dyn->getC(q,qd)*qd + dyn->getG(q) - torque - p);
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
  Y = k * Matrix::Identity(jointNo,jointNo); 
}

#endif // DISTURBANCE_OBSERVER_H
