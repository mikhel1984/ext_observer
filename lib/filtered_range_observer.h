
#ifndef FILTERED_RANGE_OBSERVER_H
#define FILTERED_RANGE_OBSERVER_H

#include "external_observer.h"
#include "iir_filter.h"

#include <iostream>

#define ID_FRangeObserver 6

class FRangeObserver : public ExternalObserver {
public:
  FRangeObserver(RobotDynamics *rd, double cutOffHz, double sampHz, double k);

  Vector getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt);

  void settings(double cutOffHz, double sampHz, double k);

private:
  FilterF1 f1;
  FilterF2 f2;
  Vector p, res;
  double shift;

}; // FRangeObserver

FRangeObserver::FRangeObserver(RobotDynamics *rd, double cutOffHz, double sampHz, double k)
  : ExternalObserver(rd,ID_FRangeObserver)
  , f1(FilterF1(cutOffHz,sampHz,jointNo))
  , f2(FilterF2(cutOffHz,sampHz,jointNo))
  , p(Vector(jointNo))
  , res(Vector(jointNo))
  , shift(k)
{
}

void FRangeObserver::settings(double cutOffHz, double sampHz, double k)
{
  f1.update(cutOffHz,sampHz);
  f2.update(cutOffHz,sampHz);
  shift = k;
}

Vector FRangeObserver::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  p = dyn->varM(q) * qd;
  //p *= shift;

  if(isRun) {
    res = f2.filt(p,dt).array().abs() + f2.getOmega() * p.array().abs() ;
    p = dyn->varFriction(qd) + dyn->varG(q) 
        + (dyn->varC(q,qd).transpose() * qd);  // reuse 
    //p *= shift;
    p = f1.filt(p,dt).array().abs();
    res += p;
  } else {
    f2.set(p);
    p = dyn->varFriction(qd) + dyn->varG(q) 
        + (dyn->varC(q,qd).transpose() * qd);  // reuse 
    //p *= shift;
    f1.set(p);
    res.setZero();
    isRun = true;
  }

  return res;
}

#endif // FILTERED_RANGE_OBSERVER_H
