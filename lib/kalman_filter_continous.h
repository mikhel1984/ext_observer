// Copyright 2020-2024 Stanislav Mikhel

/**
 * @file kalman_filter_continous.h
 *
 * @brief Discrete Kalman filter implementation based on continous system description.
 */

#ifndef KALMAN_FILTER_CONTINOUS_H
#define KALMAN_FILTER_CONTINOUS_H

#include <eigen3/Eigen/Geometry>

/**
 * @brief Discrete Kalman filter implementation.
 *
 * Use matrices for continous system.
 */
class KalmanFilterContinous {
public:
  /**
   * @brief Object constructor.
   * @param a state transition matrix.
   * @param b control input matrix.
   * @param c observation matrix.
   */
  KalmanFilterContinous(Eigen::MatrixXd& a, Eigen::MatrixXd& b, Eigen::MatrixXd& c);

  /**
   * @brief Update covariance matrices.
   *
   * Defauls are unit matrices.
   * @param q covariance of the process noise.
   * @param r covariance of the observation noise.
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
   * Time step is variable.
   * @param u control input vector.
   * @param y measured output.
   * @param dt current time step.
   * @return expected system state.
   */
  Eigen::VectorXd step(Eigen::VectorXd& u, Eigen::VectorXd& y, double dt);

  void updateR(Eigen::MatrixXd& m);

  Eigen::MatrixXd exponential(Eigen::MatrixXd& m, double dt);

private:
  /**
   * @brief Find discrete-time versions of matrices.
   *
   * @param dt current time step.
   */
  void makeDiscrete(double dt);

  // System description
  Eigen::MatrixXd Ad, Bd, Cd;
  // Covariance
  Eigen::MatrixXd R, Rd, Qd, Rupd;
  // Groups
  Eigen::MatrixXd AB, AQ, ABd, AQd;
  // Intermediate matrices
  Eigen::MatrixXd P, K, Y, I;
  // System state
  Eigen::VectorXd X;
  // matrix size
  int na, nc, nb;
};  // KalmanFilterContinous

#endif  // KALMAN_FILTER_CONTINOUS_H
