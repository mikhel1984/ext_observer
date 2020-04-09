
#ifndef FILTERED_DYNAMIC_OBSERVER_H
#define FILTERED_DYNAMIC_OBSERVER_H

#include "external_observer.h"
#include "iir_filter.h"

#include <iostream>

class FDynObserver : public ExternalObserver {
public:
  FDynObserver(RobotDynamics *rd, double cutOffHz, double sampHz);

  Vector getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt);

  void settings(double cutOffHz, double sampHz);

private:
  FilterF1 f1_1, f1_2;
  FilterF2 f2;
  Vector p, res, p0;
  double time;

}; // FDynObserver

FDynObserver::FDynObserver(RobotDynamics *rd, double cutOffHz, double sampHz)
  : ExternalObserver(rd)
  , f1_1(FilterF1(cutOffHz,sampHz,jointNo))
  , f1_2(FilterF1(cutOffHz,sampHz,jointNo))
  , f2(FilterF2(cutOffHz,sampHz,jointNo))
  , p(Vector(jointNo))
  , res(Vector(jointNo))
  , p0(Vector(jointNo))
{
}

void FDynObserver::settings(double cutOffHz, double sampHz)
{
  f1_1.update(cutOffHz,sampHz);
  f1_2.update(cutOffHz,sampHz);
  f2.update(cutOffHz,sampHz);
}

Vector FDynObserver::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  p = dyn->getM(q) * qd;
  double w = f2.getOmega();

  if(isRun) {
    time += dt;
    res = f2.filt(p) + w * p ; //- w*exp(-w*time)*p0;
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q,qd).transpose() * qd;  // reuse 
    p -= tau;
    res += f1_1.filt(p);
    //res -= f1_2.filt(tau);
  } else {
    //f1_1.clear();
    //f1_2.clear();
    //f2.clear();
    f2.set(p);
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q,qd).transpose() * qd;  // reuse 
    p -= tau;
    f1_1.set(p);
    //f1_2.set(tau);
    res.setZero();
    p0 = p;
    time = 0;
    isRun = true;
  }

  return res;
}

#endif // FILTERED_DYNAMIC_OBSERVER_H
