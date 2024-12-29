// Copyright 2020-2024 Stanislav Mikhel

#ifndef FILTERED_RANGE_OBSERVER_H
#define FILTERED_RANGE_OBSERVER_H

#include "external_observer.h"
#include "iir_filter.h"

#define ID_FRangeObserver 6

class FRangeObserver final : public ExternalObserver {
public:
  FRangeObserver(RobotDynamics *rd, double cutOffHz, double sampHz, double k);

  VectorJ getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt) override;

  void settings(double cutOffHz, double sampHz, double k);

private:
  FilterF1 f1;
  FilterF2 f2;
  VectorJ p, res;
  double shift;
};  // FRangeObserver


FRangeObserver::FRangeObserver(RobotDynamics *rd, double cutOffHz, double sampHz, double k)
  : ExternalObserver(rd, ID_FRangeObserver)
  , f1(FilterF1(cutOffHz, sampHz, jointNo))
  , f2(FilterF2(cutOffHz, sampHz, jointNo))
  , p(VectorJ(jointNo))
  , res(VectorJ(jointNo))
  , shift(k)
{
}

void FRangeObserver::settings(double cutOffHz, double sampHz, double k)
{
  f1.update(cutOffHz, sampHz);
  f2.update(cutOffHz, sampHz);
  shift = k;
}

VectorJ FRangeObserver::getExternalTorque(VectorJ& q, VectorJ& qd, VectorJ& tau, double dt)
{
  p = dyn->getM(q) * qd;
  p *= shift;

  if(isRun) {
    res = f2.filt(p, dt) + f2.getOmega() * p ;
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;  // reuse
    p *= shift;
    res += f1.filt(p, dt);
  } else {
    f2.set(p);
    p = dyn->getFriction(qd) + dyn->getG(q) - dyn->getC(q, qd).transpose() * qd;  // reuse
    p *= shift;
    f1.set(p);
    res.setZero();
    isRun = true;
  }

  return res;
}

#endif  // FILTERED_RANGE_OBSERVER_H
