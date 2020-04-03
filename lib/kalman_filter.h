/**
 * @file kalman_filter.h
 *
 * @brief Discrete Kalman filter implementation. 
 */

#ifndef KALMAN_FILTER_H
#define KALMAN_FILTER_H

#include <eigen3/Eigen/Geometry>
#include <iostream>

/**
 * @brief Discrete Kalma filter implementation.
 */
class KalmanFilter {
public:
  /** 
   * @brief Object constructor.
   * @param a ...
   * @param b ...
   * @param c ...
   */
  KalmanFilter(Eigen::MatrixXd& a, Eigen::MatrixXd& b, Eigen::MatrixXd& c);
  /**
   * @brief Update covariance matrices.
   * 
   * Defauls are unit matrices.
   * @param q ...
   * @param r ...
   */
  void setCovariance(Eigen::MatrixXd& q, Eigen::MatrixXd& r);
  /** 
   * @brief Reset filter state.
   * 
   * Define initial system state.
   * @param x0 initial system state vector.
   */
  void reset(Eigen::VectorXd& x0);
  /**
   * @brief Estimate current system state.
   *
   * Time step is constant.
   * @param u control input vector.
   * @param y measured output.
   * @return expected system state.
   */
  Eigen::VectorXd step(Eigen::VectorXd& u, Eigen::VectorXd& y);
  /**
   * @brief Estimate current system state.
   * 
   * Time step is variable.
   * @param u control input vector.
   * @param y measured output.
   * @param dt current time step.
   * @return expected system state.
   */
  Eigen::VectorXd step(Eigen::VectorXd& u, Eigen::VectorXd& y, double dt);

private:
  // System description
  Eigen::MatrixXd A, B, C;
  // Covariance 
  Eigen::MatrixXd Q, R; 
  // Intermediate matrices
  Eigen::MatrixXd P, K, Y, I, At;
  // System state
  Eigen::VectorXd X;
  // matrix size
  int nx, ny;

}; // KalmanFilter 

// Initialization
KalmanFilter::KalmanFilter(Eigen::MatrixXd& a, Eigen::MatrixXd& b, Eigen::MatrixXd& c)
  : A(a)
  , B(b)
  , C(c)
  , Q(Eigen::MatrixXd(1,1))
  , R(Eigen::MatrixXd(1,1))
  , P(Eigen::MatrixXd(1,1))
  , K(Eigen::MatrixXd(1,1))
  , Y(Eigen::MatrixXd(1,1))
  , I(Eigen::MatrixXd(1,1))
  , At(a)
  , X(Eigen::VectorXd(1))
{
  nx = A.rows();
  ny = C.rows();
  P.resize(nx,nx);
  X.resize(nx);
  Q = Eigen::MatrixXd::Identity(nx,nx);
  R = Eigen::MatrixXd::Identity(ny,ny);
  K.resize(nx,ny);
  Y.resize(ny,ny);
  I = Eigen::MatrixXd::Identity(nx,nx);
}

// Filter settings
void KalmanFilter::setCovariance(Eigen::MatrixXd& q, Eigen::MatrixXd& r)
{
  Q = q; 
  R = r;
}

// Define new initial state
void KalmanFilter::reset(Eigen::VectorXd& x0)
{
  X = x0;
  P = Eigen::MatrixXd::Zero(nx,nx);
}

// State estimation for constant time step
Eigen::VectorXd KalmanFilter::step(Eigen::VectorXd& u, Eigen::VectorXd& y)
{
  // predict 
  X = A * X + B * u;                     // i | i-1
  P = A * P * A.transpose() + Q;         // i | i-1
  // update 
  Y = C * P * C.transpose() + R;         // i
  K = P * C.transpose() * Y.inverse();   // i

  X = X + K * (y - C * X);               // i | i
  P = (I - K * C) * P;                   // i | i
  
  return X;
}

// State estimation for variable time step
Eigen::VectorXd KalmanFilter::step(Eigen::VectorXd& u, Eigen::VectorXd& y, double dt)
{
  // find current A matrix
  At = I + dt*A;
  // predict 
  X = At * X + dt * B * u;               // i | i-1
  P = At * P * At.transpose() + Q;       // i | i-1
  // update 
  Y = C * P * C.transpose() + R;         // i
  K = P * C.transpose() * Y.inverse();   // i

  X = X + K * (y - C * X);               // i | i
  P = (I - K * C) * P;                   // i | i
  
  return X;
}

#endif // KALMAN_FILTER_H
