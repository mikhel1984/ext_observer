
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#include "../src/kalman_filter.h"

#define RANDOM ((rand() % 1000)/1000.0)
#define FNAME "estimation.csv"
#define DELTA 0.02

using namespace Eigen;

int main(int argc, char **argv)
{
  // system 
  MatrixXd A(2,2), B(2,1), C(1,2);
  A << -5, 1, 1, -4;
  B << 0.5, 0.1;
  C << 2, 1;

  // to discret system
  A = MatrixXd::Identity(2,2) + DELTA * A;
  B = DELTA * B;

  // control impact 
  VectorXd u(1);

  // initial state 
  VectorXd x(2,1), y(1), yest(1); 
  x.setZero();

  // filter
  KalmanFilter filter(A,B,C);
  filter.reset(x);

  MatrixXd Q(2,2), R(1,1);
  Q << 0.003, 0, 0, 0.003;
  R << 1;
  filter.setCovariance(Q,R);

  // save to file 
  std::ofstream file;
  file.open(FNAME); 

  for(double t = 0; t < 3; t += DELTA) {
    // next state
    u(0) = sin(2*t);
    x = A * x + B * u; 
    y = C * x;
    y(0) += 0.02*(2*RANDOM - 1);
    // estimation
    yest = C * filter.step(u,y);
    // write
    file << t << "," << y(0) << "," << yest(0) << std::endl;
  }

  file.close();

  return 0;
}
