#ifndef IIR_FILTER_H
#define IIR_FILTER_H

#include <cmath>
#include <eigen3/Eigen/Geometry>
#include <iostream>

class FilterIIR {
public:
  FilterIIR(double cutOffHz, double sampHz) 
  { update(cutOffHz,sampHz); } 

  virtual Eigen::VectorXd filt(Eigen::VectorXd& x) = 0;

  virtual void clear() = 0;

  void update(double cutOffHz, double sampHz) 
  { omega = tan(3.1415926*cutOffHz / sampHz); }  // skip 2 / 2 

  double getOmega() { return omega; }

protected:
  double omega;
}; // FilterIIR

// H(s) = w / (s + w)
class FilterF1 : public FilterIIR {
public:
  FilterF1(double cutOffHz, double sampHz, int N)
  : FilterIIR(cutOffHz, sampHz) 
  , x1(Eigen::VectorXd(N))
  , y1(Eigen::VectorXd(N))
  { }

  Eigen::VectorXd filt(Eigen::VectorXd& x);

  void clear() { x1.setZero(); y1.setZero(); }


private:
  Eigen::VectorXd x1, y1;

}; // FilterF1

// H(s) = -w^2 / (s + w)
class FilterF2 : public FilterIIR {
public:
  FilterF2(double cutOffHz, double sampHz, int N)
  : FilterIIR(cutOffHz, sampHz) 
  , x1(Eigen::VectorXd(N))
  , y1(Eigen::VectorXd(N))
  { }

  Eigen::VectorXd filt(Eigen::VectorXd& x);

  void clear() { x1.setZero(); y1.setZero(); }

private:
  Eigen::VectorXd x1, y1;

}; // FilterF2

class FilterButterworth : public FilterIIR {
public:
  FilterButterworth(double cutOffHz, double sampHz, int N)
  : FilterIIR(cutOffHz, sampHz) 
  , x1(Eigen::VectorXd(N))
  , x2(Eigen::VectorXd(N))
  , y1(Eigen::VectorXd(N))
  , y2(Eigen::VectorXd(N))
  , res(Eigen::VectorXd(N))
  { }

  Eigen::VectorXd filt(Eigen::VectorXd& x);

  void clear() { x1.setZero(); x2.setZero(); y1.setZero(); y2.setZero(); }

private:
  Eigen::VectorXd x1, x2, y1, y2, res;
}; // FilterButterworth

Eigen::VectorXd FilterF1::filt(Eigen::VectorXd& x)
{
  double k1 = (1-omega)/(1+omega);
  double k2 = omega / (1+omega);

  //std::cout << omega << " " << k1 << " " << k2 << std::endl;

  y1 = k1*y1 + k2*(x + x1);
  x1 = x;

  return y1;
}

Eigen::VectorXd FilterF2::filt(Eigen::VectorXd& x)
{
  double k1 = (1-omega)/(1+omega);
  double k2 = -omega*omega/(1+omega);

  y1 = k1*y1 + k2*(x + x1);
  x1 = x;

  return y1;
}

Eigen::VectorXd FilterButterworth::filt(Eigen::VectorXd& x)
{
  double denom = 1 + sqrt(2)*omega + omega*omega;
  double k1 = (2*(1-omega*omega))/denom;
  double k2 = (sqrt(2)*omega - 1 + omega*omega)/denom;
  double k3 = omega*omega/denom; // == k5
  double k4 = 2*k3/denom;

  res = k1*y1 + k2*y2 + k3*(x + x2) + k4*x1;

  y2 = y1; y1 = res;
  x2 = x1; x1 = x;

  return res;
}

#endif
