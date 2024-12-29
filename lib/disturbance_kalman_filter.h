// Copyright 2020-2024 Stanislav Mikhel

/**
 * @file disturbance_kalman_filter.h
 *
 * @brief Disturbance kalman filter observer.
 *
 * Expected explicit robot dynamics matrices.
 */
#ifndef DISTURBANCE_KALMAN_FILTER_H
#define DISTURBANCE_KALMAN_FILTER_H

#include "external_observer.h"
#include "kalman_filter.h"

#define ID_DKalmanObserver 1

/**
 * @brief Disturbance Kalman filter from Hu et. al.
 */
class DKalmanObserver final : public ExternalObserver {
public:
  /**
   * @brief Object constructor
   * @param rd pointer to the robot object.
   * @param s disturbance dynamics matrix.
   * @param h disturbance observation matrix.
   * @param q covariance of the process noise.
   * @param r covariance of the observation noise.
   */
  DKalmanObserver(RobotDynamics* rd, MatrixJ& s, MatrixJ& h, MatrixJ& q, MatrixJ& r);

  /**
   * @brief Object destructor.
   */
  virtual ~DKalmanObserver();

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
   * @param q covariance of the process noise.
   * @param r covariance of the observation noise.
   */
  void settings(MatrixJ& q, MatrixJ& r);

private:
  // disturbance observation
  MatrixJ H;
  // state, input and momentum
  VectorJ X, u, p;
  // Kalman filter object
  KalmanFilter *filter;
};  // DKalmanObserver


// Initialization
DKalmanObserver::DKalmanObserver(RobotDynamics* rd, MatrixJ& s, MatrixJ& h, MatrixJ& q, MatrixJ& r)
  : ExternalObserver(rd, ID_DKalmanObserver)
  , H(h)
  , X(VectorJ(2*jointNo))
  , u(VectorJ(jointNo))
  , p(VectorJ(jointNo))
  , filter(0)
{
  // prepare matrices
  MatrixJ A(2*jointNo, 2*jointNo), B(2*jointNo, jointNo), C(jointNo, 2*jointNo);
  B.setZero();
  C.setZero();
  // B & C matrices
  for(int i = 0; i < jointNo; i++) {
    B(i, i) = 1;
    C(i, i) = 1;
  }
  // A matrix
  A.setZero();
  A.block(      0, jointNo, jointNo, jointNo) = h;
  A.block(jointNo, jointNo, jointNo, jointNo) = s;
  // make filter
  filter = new KalmanFilter(A, B, C);
  filter->setCovariance(q, r);
}

// Clear dynamically allocated memory
DKalmanObserver::~DKalmanObserver()
{
  delete filter;
}

// Torque estimation
VectorJ DKalmanObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;
  u = tau - dyn->getG(q) - dyn->getFriction(qd);
  u += dyn->getC(q, qd).transpose() * qd;

  if(isRun) {
    X = filter->step(u, p, dt);
  } else {
    // prepare X0
    X.setZero();
    for(int i = 0; i < jointNo; i++) X(i) = p(i);
    // reset
    filter->reset(X);
    isRun = true;
  }
  // disturbance = H * omega   (reuse variable)
  p = H * X.block(jointNo, 0, jointNo, 1);

  return p;
}

// settings
void DKalmanObserver::settings(MatrixJ& Q, MatrixJ& R)
{
  filter->setCovariance(Q, R);
}

#endif  // DISTURBANCE_KALMAN_FILTER_H
