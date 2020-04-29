#ifndef IIR_FILTER_H
#define IIR_FILTER_H

#include <cmath>
#include <eigen3/Eigen/Dense>
//#include <iostream>

#define SQ2 1.4142135623731

class FilterIIR {
public:

  virtual Eigen::VectorXd filt(Eigen::VectorXd& x) = 0;

  virtual void update(double cutOff, double sampTime) = 0;

}; // FilterIIR

// H(s) = w / (s + w)
class FilterF1 : public FilterIIR {
public:
  FilterF1(double cutOff, double sampTime, int N)
  : FilterIIR()
  , x1(Eigen::VectorXd(N))
  , y1(Eigen::VectorXd(N))
  { update(cutOff,sampTime); }

  Eigen::VectorXd filt(Eigen::VectorXd& x);
  
  Eigen::VectorXd filt(Eigen::VectorXd& x, double dt);

  void set(Eigen::VectorXd& x0) { x1 = x0; y1 = x0; }

  void update(double cutOff, double sampTime);
private:
  Eigen::VectorXd x1, y1;
  double k1, k2, cut;

}; // FilterF1

// H(s) = -w^2 / (s + w)
class FilterF2 : public FilterIIR {
public:
  FilterF2(double cutOff, double sampTime, int N)
  : FilterIIR()
  , x1(Eigen::VectorXd(N))
  , y1(Eigen::VectorXd(N))
  { update(cutOff,sampTime); }

  Eigen::VectorXd filt(Eigen::VectorXd& x);
  
  Eigen::VectorXd filt(Eigen::VectorXd& x, double dt);

  void set(Eigen::VectorXd& x0) { x1 = x0; y1 = -f2*omega*x0; }

  double getOmega() { return cut; }

  void update(double cutOff, double sampTime);
private:
  Eigen::VectorXd x1, y1;
  double f2, omega, k1, k2, cut;

}; // FilterF2


class FilterButterworth : public FilterIIR {
public:
  FilterButterworth(double cutOff, double sampTime, int N)
  : FilterIIR() 
  , x1(Eigen::VectorXd(N))
  , x2(Eigen::VectorXd(N))
  , y1(Eigen::VectorXd(N))
  , y2(Eigen::VectorXd(N))
  , res(Eigen::VectorXd(N))
  { update(cutOff,sampTime); }

  Eigen::VectorXd filt(Eigen::VectorXd& x);
  
  void update(double cutOff, double sampTime);

  
private:
  Eigen::VectorXd x1, x2, y1, y2, res;
  double kx0, kx1, ky1, ky2;
  
}; // FilterButterworth


class FilterLowPass : public FilterIIR {
public: 
  FilterLowPass(double cutOff, double sampTime, int N)
  : FilterIIR()
  , x1(Eigen::VectorXd(N))
  , y1(Eigen::VectorXd(N))
  { update(cutOff,sampTime); }
  
  void update(double cutOff, double sampTime);
  
  Eigen::VectorXd filt(Eigen::VectorXd& x);

private:
  Eigen::VectorXd x1, y1;
  double k1, k2;
  
}; // FilterLowPass

class FilterHighPass : public FilterIIR {
public:
  FilterHighPass(double cutOff, double sampTime, int N)
  : FilterIIR()
  , x1(Eigen::VectorXd(N))
  , y1(Eigen::VectorXd(N))
  { update(cutOff,sampTime); }
  
  void update(double cutOff, double sampTime);
  
  Eigen::VectorXd filt(Eigen::VectorXd& x);

private:
  Eigen::VectorXd x1, y1;
  double k1, k2;

}; // FilterHighPass

Eigen::VectorXd FilterF1::filt(Eigen::VectorXd& x)
{

  y1 = k1*y1 + k2*(x + x1);
  x1 = x;

  return y1;
}

void FilterF1::update(double cutOff, double sampTime) 
{
  double omega = tan(cutOff*sampTime*0.5);  
  k1 = (1-omega) / (1+omega);
  k2 = omega / (1 + omega);
  cut = cutOff;
}

Eigen::VectorXd FilterF1::filt(Eigen::VectorXd& x, double dt)
{
  double omega = tan(cut*dt*0.5);  
  k1 = (1-omega) / (1+omega);
  k2 = omega / (1 + omega);
  
  return filt(x);
}

Eigen::VectorXd FilterF2::filt(Eigen::VectorXd& x)
{
  y1 = k1*y1 + k2*(x + x1);
  x1 = x;

  return y1;
}


void FilterF2::update(double cutOff, double sampTime) 
{
  cut = cutOff;
  f2 = 2/sampTime; 
  omega = tan(cutOff*sampTime*0.5);
  k1 = (1-omega)/(1+omega);
  k2 = -f2*omega*omega/(1+omega);
}

Eigen::VectorXd FilterF2::filt(Eigen::VectorXd& x, double dt)
{
  f2 = 2/dt; 
  omega = tan(cut*dt*0.5);
  k1 = (1-omega)/(1+omega);
  k2 = -f2*omega*omega/(1+omega);
  
  return filt(x);
}

void FilterButterworth::update(double cutOff, double sampTime)
{
  double omega = tan(cutOff * sampTime * 0.5); 
  double omega2 = omega * omega;
  double denom = 1 + SQ2*omega + omega2; 
  
  ky1 = 2*(1-omega2)/denom;
  ky2 = (SQ2*omega-1-omega2)/denom;
  kx0 = omega2 / denom;            // kx2 == kx0
  kx1 = 2*kx0;  
}

Eigen::VectorXd FilterButterworth::filt(Eigen::VectorXd& x)
{
  res = ky1*y1 + ky2*y2 + kx0*(x + x2) + kx1*x1;
  x2 = x1; x1 = x;
  y2 = y1; y1 = res;
  
  return res;
}

void FilterLowPass::update(double cutOff, double sampTime) 
{ 
  double omega = tan(cutOff*sampTime*0.5); 
  k1 = (1-omega)/(1+omega); 
  k2 = omega/(1+omega); 
}

Eigen::VectorXd FilterLowPass::filt(Eigen::VectorXd& x)
{
  y1 *= k1;
  y1 += k2*(x + x1);
  x1 = x;
  return y1;
}

void FilterHighPass::update(double cutOff, double sampTime)
{
  double omega = tan(cutOff*sampTime*0.5); 
  k1 = (1-omega) / (1+omega);
  k2 = 1 / (1+omega); 
}

Eigen::VectorXd FilterHighPass::filt(Eigen::VectorXd& x)
{
  y1 *= k1;
  y1 += k2*(x - x1);
  x1 = x;
  return y1;
}


#endif // IIR_FILTER_H
