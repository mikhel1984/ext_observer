
#ifndef FILTERED_DYNAMIC_OBSERVER_RNEA_H
#define FILTERED_DYNAMIC_OBSERVER_RNEA_H

#include "external_observer.h"
#include "iir_filter.h"

#include <iostream>

#define ID_FDynObserverRnea 22

class FDynObserverRnea : public ExternalObserverRnea {
public:
  FDynObserverRnea(RobotDynamicsRnea *rd, double cutOffHz, double sampHz);

  Vector getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt);

  void settings(double cutOffHz, double sampHz);

private:
  FilterF1 f1;
  FilterF2 f2;
  Vector p, res, zero;

}; // FDynObserver

FDynObserverRnea::FDynObserverRnea(RobotDynamicsRnea *rd, double cutOffHz, double sampHz)
  : ExternalObserverRnea(rd,ID_FDynObserverRnea)
  , f1(FilterF1(cutOffHz,sampHz,jointNo))
  , f2(FilterF2(cutOffHz,sampHz,jointNo))
  , p(Vector(jointNo))
  , res(Vector(jointNo))
  , zero(Vector::Zero(jointNo))
{
}

void FDynObserverRnea::settings(double cutOffHz, double sampHz)
{
  f1.update(cutOffHz,sampHz);
  f2.update(cutOffHz,sampHz);
}

Vector FDynObserverRnea::getExternalTorque(Vector& q, Vector& qd, Vector& tau, double dt)
{
  p = dyn->rnea(q,zero,qd);     // M * qd

  if(isRun) {
    res = f2.filt(p,dt) + f2.getOmega() * p ;
    p = dyn->getFriction(qd) + dyn->rnea(q,zero,zero,GRAVITY) - dyn->tranCqd(q,qd);  // reuse 
    p -= tau;
    res += f1.filt(p,dt);
  } else {
    f2.set(p);
    p = dyn->getFriction(qd) + dyn->rnea(q,zero,zero,GRAVITY) - dyn->tranCqd(q,qd);  // reuse 
    p -= tau;
    f1.set(p);
    res.setZero();
    isRun = true;
  }

  return res;
}

#endif // FILTERED_DYNAMIC_OBSERVER_RNEA_H
