// Copyright 2020-2024 Stanislav Mikhel

#include "external_observer.h"

RobotDynamicsRnea::RobotDynamicsRnea()
  : RobotDynamicsBase()
  , _qext(VectorJ(1))
  , _p0(VectorJ(1))
  , _zero(VectorJ(1))
  , _sum(VectorJ(1))
  , _M(MatrixJ(1, 1))
{
}

// Find C^T * qd
VectorJ RobotDynamicsRnea::tranCqd(VectorJ& q, VectorJ& qd)
{
  // TODO: call it once
  int N = jointNo();
  _zero.resize(N);
  _zero.setZero();
  _p0 = rnea(q, _zero, qd);  // M * qd
  _sum.setZero();

  for(int i = 0; i < N; i++) {
    _qext = q;
    _qext(i) += delta;
    _sum += (rnea(_qext, _zero, qd) - _p0) * (qd(i) / delta);
  }
  _sum -= rnea(q, qd, _zero);  // M'*qd - C*qd

  return _sum;
}

// Create M matrix from sequence of RNEA calls
MatrixJ RobotDynamicsRnea::getM(VectorJ& q)
{
  // TODO: call it once
  int N = jointNo();
  _zero.resize(N);
  _zero.setZero();
  _M.resize(N, N);
  _qext.resize(N);

  for(int i = 0; i < N; i++) {
    _qext.setZero();  // use _qext to choose column
    _qext(i) = 1;
    _p0 = rnea(q, _zero, _qext);  // use _p0 to save matrix column
    for(int j = 0; j < N; j++)
      _M(j, i) = _p0(j);
  }

  return _M;
}
