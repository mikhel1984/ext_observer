/**
 * @file kalman_filter_continous.h
 *
 * @brief Discrete Kalman filter implementation based on continous system description. 
 */

#ifndef KALMAN_FILTER_CONTINOUS_H
#define KALMAN_FILTER_CONTINOUS_H

#include <eigen3/unsupported/Eigen/MatrixFunctions>
//#include <iostream>

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
  Eigen::MatrixXd R, Rd, Qd; 
  // Groups
  Eigen::MatrixXd AB, AQ, ABd, AQd;
  // Intermediate matrices
  Eigen::MatrixXd P, K, Y, I;
  // System state
  Eigen::VectorXd X;
  // matrix size
  int na, nc, nb;

}; // KalmanFilterContinous 

// Initialization
KalmanFilterContinous::KalmanFilterContinous(Eigen::MatrixXd& a, Eigen::MatrixXd& b, Eigen::MatrixXd& c)
  : Ad(a)
  , Bd(b)
  , Cd(c)
  , R(Eigen::MatrixXd(1,1))
  , Rd(Eigen::MatrixXd(1,1))
  , Qd(Eigen::MatrixXd(1,1))
  , AB(Eigen::MatrixXd(1,1))
  , AQ(Eigen::MatrixXd(1,1))
  , ABd(Eigen::MatrixXd(1,1))
  , AQd(Eigen::MatrixXd(1,1))
  , P(Eigen::MatrixXd(1,1))
  , K(Eigen::MatrixXd(1,1))
  , Y(Eigen::MatrixXd(1,1))
  , I(Eigen::MatrixXd(1,1))
  , X(Eigen::VectorXd(1))
{
  na = a.rows();
  nc = c.rows();
  nb = b.cols();
  P.resize(na,na);
  X.resize(na);
  Qd = Eigen::MatrixXd::Identity(na,na);
  R  = Eigen::MatrixXd::Identity(nc,nc);
  Rd = R;
  K.resize(na,nc);
  Y.resize(nc,nc);
  I = Eigen::MatrixXd::Identity(na,na);
  // prepare matrix combinations
  // [A B;0 0]
  AB.resize(na+nb,na+nb);
  AB.setZero(); 
  AB.block(0,0,na,na) = a;
  AB.block(0,na,na,nb) = b;
  ABd.resize(na+nb,na+nb);
  // [-A Q; 0 A'] 
  AQ.resize(na+na,na+na);
  AQ.block(0,0,na,na) = -a;
  AQ.block(0,na,na,na) = Qd;
  AQ.block(na,na,na,na) = a.transpose();
  AQd.resize(na+na,na+na);
}

// Filter settings
void KalmanFilterContinous::setCovariance(Eigen::MatrixXd& q, Eigen::MatrixXd& r)
{
  AQ.block(0,na,na,na) = q;
  R = r;
}

// Define new initial state
void KalmanFilterContinous::reset(Eigen::VectorXd& x0)
{
  X = x0;
  P = Eigen::MatrixXd::Zero(na,na);  
}

// State estimation for constant time step
Eigen::VectorXd KalmanFilterContinous::step(Eigen::VectorXd& u, Eigen::VectorXd& y, double dt)
{
  makeDiscrete(dt);

  // predict 
  X = Ad * X + Bd * u;                       // i | i-1
  P = Ad * P * Ad.transpose() + Qd;           // i | i-1
  // update 
  Y = Cd * (P * Cd.transpose()) + Rd;         // i
  K = P * (Cd.transpose() * Y.inverse());   // i

  X += K * (y - Cd * X);                    // i | i
  P = (I - K * Cd) * P;                     // i | i
  
  return X;
}

void KalmanFilterContinous::makeDiscrete(double dt)
{
  // find At and Bt
  ABd = (AB*dt).exp(); 
  Ad = ABd.block(0,0,na,na);
  Bd = ABd.block(0,na,na,nb);
  // find Qt and Rt
  Rd = R / dt;
  AQd = (AQ*dt).exp();
  Qd = AQd.block(na,na,na,na).transpose() * AQd.block(0,na,na,na);
}

#endif // KALMAN_FILTER_CONTINOUS_H
