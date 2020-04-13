
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
  FilterF1 f1;
  FilterF2 f2;
  Vector p, res;

}; // FDynObserver

FDynObserver::FDynObserver(RobotDynamics *rd, double cutOffHz, double sampHz)
  : ExternalObserver(rd)
  , f1(FilterF1(cutOffHz,sampHz,jointNo))
  , f2(FilterF2(cutOffHz,sampHz,jointNo))
  , p(Vector(jointNo))
  , res(Vector(jointNo))
{
}

void FDynObserver::settings(double cutOffHz, double sampHz)
{
  f1.update(cutOffHz,sampHz);
  f2.update(cutOffHz,sampHz);
}

Vector FDynObserver::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  p = dyn->getM(q) * qd;

  if(isRun) {
    res = f2.filt(p,dt) + f2.getOmega() * p ;
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q,qd).transpose() * qd;  // reuse 
    p -= tau;
    res += f1.filt(p,dt);
  } else {
    f2.set(p);
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q,qd).transpose() * qd;  // reuse 
    p -= tau;
    f1.set(p);
    res.setZero();
    isRun = true;
  }

  return res;
}

#endif // FILTERED_DYNAMIC_OBSERVER_H
