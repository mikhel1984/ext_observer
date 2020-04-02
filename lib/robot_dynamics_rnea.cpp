
#include "external_observer.h" 

#define DELTA_INIT 1E-7

RobotDynamicsRnea::RobotDynamicsRnea()
  : RobotDynamicsBase()
  , _qext(Vector(1))
  , _p0(Vector(1))
  , _zero(Vector(1))
  , _sum(Vector(1))
  , _M(Matrix(1,1))
  , delta(DELTA_INIT) 
{
}

// Find C^T * qd
Vector RobotDynamicsRnea::tranCqd(Vector& q, Vector& qd)
{
  // TODO: call it once
  int N = jointNo();
  _zero.resize(N); _zero.setZero();
  _p0 = rnea(q,_zero,qd);  // M * qd
  _sum.setZero();

  for(int i = 0; i < N; i++) {
    _qext = q;
    _qext(i) += delta; 
    _sum += (rnea(_qext,_zero,qd) - _p0) * (qd(i) / delta);
  }
  _sum -= rnea(q,qd,_zero); // M'*qd - C*qd

  return _sum;
}

// Create M matrix from sequence of RNEA calls
Matrix RobotDynamicsRnea::getM(Vector& q)
{
  // TODO: call it once
  int N = jointNo();
  _zero.resize(N); _zero.setZero();
  _M.resize(N,N); 
  _qext.resize(N);
  
  for(int i = 0; i < N; i++) {
    _qext.setZero();  // use _qext to choose column
    _qext(i) = 1;   
    _p0 = rnea(q,_zero,_qext); // use _p0 to save matrix column
    for(int j = 0; j < N; j++)
      _M(j,i) = _p0(j);
  }
  
  return _M;
}
