// Copyright 2020-2024 Stanislav Mikhel

/**
 * @file disturbance_kalman_filter.h
 *
 * @brief Disturbance kalman filter observer.
 *
 * Expected explicit robot dynamics matrices.
 */
#ifndef DISTURBANCE_KALMAN_FILTER_RNEA_H
#define DISTURBANCE_KALMAN_FILTER_RNEA_H

#include "external_observer.h"
#include "kalman_filter.h"

#define ID_DKalmanObserverRnea 21

/**
 * @brief Disturbance Kalman filter from Hu et. al.
 *
 * Use RNEA technique for dynamics.
 */
class DKalmanObserverRnea final : public ExternalObserverRnea {
public:
  /**
   * @brief Object constructor
   * @param rd pointer to the robot object.
   * @param s disturbance dynamics matrix.
   * @param h disturbance observation matrix.
   * @param q covariance of the process noise.
   * @param r covariance of the observation noise.
   */
  DKalmanObserverRnea(RobotDynamicsRnea* rd, MatrixJ& s, MatrixJ& h, MatrixJ& q, MatrixJ& r);

  /**
   * @brief Object destructor.
   */
  virtual ~DKalmanObserverRnea();

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
  VectorJ X, u, p, zero;
  // Kalman filter object
  KalmanFilter *filter;
};  // DKalmanObserverRnea

// Initialization
DKalmanObserverRnea::DKalmanObserverRnea(
  RobotDynamicsRnea* rd, MatrixJ& s, MatrixJ& h, MatrixJ& q, MatrixJ& r)
  : ExternalObserverRnea(rd, ID_DKalmanObserverRnea)
  , H(h)
  , X(VectorJ(2*jointNo))
  , u(VectorJ(jointNo))
  , p(VectorJ(jointNo))
  , zero(VectorJ::Zero(jointNo))
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
DKalmanObserverRnea::~DKalmanObserverRnea()
{
  delete filter;
}

// Torque estimation
VectorJ DKalmanObserverRnea::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->rnea(q, zero, qd);     // M * qd
  u = tau - dyn->rnea(q, zero, zero, GRAVITY) - dyn->getFriction(qd);
  u += dyn->tranCqd(q, qd);      // C' * qd

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
void DKalmanObserverRnea::settings(MatrixJ& Q, MatrixJ& R)
{
  filter->setCovariance(Q, R);
}

#endif  // DISTURBANCE_KALMAN_FILTER_RNEA_H
