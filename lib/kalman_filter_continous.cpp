
//#include <eigen3/unsupported/Eigen/MatrixFunctions>
#include "kalman_filter_continous.h"

#define EXP_TERMS 5

// Initialization
KalmanFilterContinous::KalmanFilterContinous(Eigen::MatrixXd& a, Eigen::MatrixXd& b, Eigen::MatrixXd& c)
  : Ad(a)
  , Bd(b)
  , Cd(c)
  , R(Eigen::MatrixXd(1,1))
  , Rd(Eigen::MatrixXd(1,1))
  , Qd(Eigen::MatrixXd(1,1))
  , Rupd(Eigen::MatrixXd(1,1))
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
  Rupd = R;
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
  //ABd = (AB*dt).exp(); 
  ABd = exponential(AB,dt);
  Ad = ABd.block(0,0,na,na);
  Bd = ABd.block(0,na,na,nb);
  // find Qt and Rt
  Rd = Rupd / dt;
  //AQd = (AQ*dt).exp();
  AQd = exponential(AQ,dt);
  Qd = AQd.block(na,na,na,na).transpose() * AQd.block(0,na,na,na);
}

void KalmanFilterContinous::updateR(Eigen::MatrixXd& m)
{
  Rupd = m * R * m.transpose();
}

Eigen::MatrixXd KalmanFilterContinous::exponential(Eigen::MatrixXd& m, double dt)
{
  Eigen::MatrixXd res = Eigen::MatrixXd::Identity(m.rows(),m.cols());
  Eigen::MatrixXd acc = res; 
  
  for(int i = 1; i <= EXP_TERMS; i++) {
    acc *= m * (dt/i);
    res += acc;
  }

  return res;  
}