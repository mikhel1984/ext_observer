#include "../tests/double_link.h"
#include "../lib/momentum_observer.h"
#include "../lib/disturbance_observer.h"
#include "../lib/sliding_mode_observer.h"
#include "../lib/disturbance_kalman_filter.h"
#include "../lib/filtered_dyn_observer.h"
#include "../lib/filtered_range_observer.h"
#include "../lib/disturbance_kalman_filter_exp.h"

#include "observers.h"

#define ARRAY_LEN 30  // max number of observers

// robot dynamics
static DoubleLink robot;
// observers
static ExternalObserver *observer[ARRAY_LEN] = {0};
// position to add new element
static int _nextIndex = 0;

void reset(int ind) 
{
  if(ind >= 0 && ind < _nextIndex) 
    observer[ind]->reset();
}

void freeAll()
{
  MomentumObserver *mo;
  DisturbanceObserver *dis;
  SlidingModeObserver *sm;
  DKalmanObserver *dk;
  FDynObserver *fd;
  FRangeObserver *fr;
  DKalmanObserverExp *dke;
  for(int i = 0; i < _nextIndex; i++) {
    switch (observer[i]->type() ) {
    case ID_MomentumObserver: 
      mo = (MomentumObserver*) observer[i];
      delete mo;
      break;
    case ID_DisturbanceObserver:
      dis = (DisturbanceObserver*) observer[i];
      delete dis;
      break;
    case ID_SlidingModeObserver:
      sm = (SlidingModeObserver*) observer[i];
      delete sm;
      break;
    case ID_DKalmanObserver:
      dk = (DKalmanObserver*) observer[i];
      delete dk;
      break;
    case ID_FDynObserver:
      fd = (FDynObserver*) observer[i];
      delete fd;
      break;
    case ID_FRangeObserver:
      fr = (FRangeObserver*) observer[i];
      delete fr;
      break;
    case ID_DKalmanObserverExp:
      dke = (DKalmanObserverExp*) observer[i];
      delete dke;
      break;
    }
    observer[i] = 0;
  }
  _nextIndex = 0;
}

int getExternalTorque(int ind, double* ext, double *q, double *qd, double *tau, double dt) 
{
  if(ind < 0 || ind >= _nextIndex) return ERR_WRONG_INDEX;
  ExternalObserver* ob = observer[ind];

  Vector vext(JOINT_NO), vq(JOINT_NO), vqd(JOINT_NO), vtau(JOINT_NO);

  for(int i = 0; i < JOINT_NO; i++) {
    vq(i) = q[i];
    vqd(i) = qd[i];
    vtau(i) = tau[i];
  }

  vext = ob->getExternalTorque(vq,vqd,vtau,dt);
  
  for(int i = 0; i < JOINT_NO; i++) {
    ext[i] = vext(i);
  }

  return 0;
}

void getRobotTorque(double *tau, double *q, double *qd, double *q2d)
{
  Vector vq(JOINT_NO), vqd(JOINT_NO), vq2d(JOINT_NO), res(JOINT_NO);
  for(int i = 0; i < JOINT_NO; i++) {
    vq(i) = q[i];
    vqd(i) = qd[i];
    vq2d(i) = q2d[i];
  }
  
  res = robot.getM(vq) * vq2d + robot.getC(vq,vqd) * vqd + robot.getG(vq) + robot.getFriction(vqd); 
  
  for(int i = 0; i < JOINT_NO; i++) 
    tau[i] = res(i);
}

int configMomentumObserver(int ind, double *k)
{
  // prepare
  MomentumObserver* ptr;
  Vector vk(JOINT_NO);
  for(int i = 0; i < JOINT_NO; i++) 
    vk(i) = k[i];

  if(ind == ADD_NEW) {
    if(_nextIndex == ARRAY_LEN) return ERR_NO_SLOTS;  // can't add new
    // create new
    ptr = new MomentumObserver(&robot, vk);
    observer[_nextIndex] = ptr;
    return _nextIndex++;    // new element _nextIndex
  } else if(ind < ADD_NEW || ind >= ARRAY_LEN) {
    return ERR_WRONG_INDEX;
  }
  
  // update state
  ptr = (MomentumObserver*) observer[ind];
  if(ptr->type() != ID_MomentumObserver)
    return ERR_WRONG_TYPE;
  ptr->settings(vk);

  return ind; // ok
}

int configDisturbanceObserver(int ind, double sigma, double xeta, double beta)
{
  // prepare
  DisturbanceObserver* ptr;

  if(ind == ADD_NEW) {
    if(_nextIndex == ARRAY_LEN) return ERR_NO_SLOTS;  // can't add new
    // create new
    ptr = new DisturbanceObserver(&robot,sigma,xeta,beta);
    observer[_nextIndex] = ptr;
    return _nextIndex++;    // new element _nextIndex
  } else if(ind < ADD_NEW || ind >= ARRAY_LEN) {
    return ERR_WRONG_INDEX;
  }
  // update state
  ptr = (DisturbanceObserver*) observer[ind];
  if(ptr->type() != ID_DisturbanceObserver)
    return ERR_WRONG_TYPE;
  ptr->settings(sigma,xeta,beta);

  return ind; // ok
}

int configSlidingModeObserver(int ind, double *T1, double *S1, double *T2, double *S2)
{
  // prepare
  SlidingModeObserver* ptr;
  Vector vT1(JOINT_NO), vS1(JOINT_NO), vT2(JOINT_NO), vS2(JOINT_NO);
  for(int i = 0; i < JOINT_NO; i++) {
    vT1(i) = T1[i];
    vS1(i) = S1[i];
    vT2(i) = T2[i];
    vS2(i) = S2[i];
  }

  if(ind == ADD_NEW) {
    if(_nextIndex == ARRAY_LEN) return ERR_NO_SLOTS;  // can't add new
    // create new
    ptr = new SlidingModeObserver(&robot,vT1,vS1,vT2,vS2);
    observer[_nextIndex] = ptr;
    return _nextIndex++;    // new element _nextIndex
  } else if(ind < ADD_NEW || ind >= ARRAY_LEN) {
    return ERR_WRONG_INDEX;
  }

  ptr = (SlidingModeObserver*) observer[ind];
  if(ptr->type() != ID_SlidingModeObserver)
    return ERR_WRONG_TYPE;
  ptr->settings(vT1,vS1,vT2,vS2);

  return ind;
}

int configDistKalmanObserver(int ind, double *S, double *H, double *Q, double *R)
{
  DKalmanObserver* ptr;
  Matrix mS(JOINT_NO,JOINT_NO), mH(JOINT_NO,JOINT_NO), mQ(2*JOINT_NO,2*JOINT_NO), mR(JOINT_NO,JOINT_NO);
  int k = 0;
  for(int c = 0; c < JOINT_NO; c++) {
    for(int r = 0; r < JOINT_NO; r++,k++) {
      mS(r,c) = S[k];
      mH(r,c) = H[k];
      mR(r,c) = R[k];
    }
  }
  k = 0;
  for(int c = 0; c < 2*JOINT_NO; c++) {
    for(int r = 0; r < 2*JOINT_NO; r++,k++) 
      mQ(r,c) = Q[k];
  }

  if(ind == ADD_NEW) {
    if(_nextIndex == ARRAY_LEN) return ERR_NO_SLOTS;  // can't add new
    // create new
    ptr = new DKalmanObserver(&robot,mS,mH,mQ,mR);
    observer[_nextIndex] = ptr;
    return _nextIndex++;    // new element _nextIndex
  } else if(ind < ADD_NEW || ind >= ARRAY_LEN) {
    return ERR_WRONG_INDEX;
  }

  ptr = (DKalmanObserver*) observer[ind];
  if(ptr->type() != ID_DKalmanObserver)
    return ERR_WRONG_TYPE;
  ptr->settings(mQ,mR);

  return ind;
}

int configDistKalmanObserverExp(int ind, double *S, double *H, double *Q, double *R)
{
  DKalmanObserverExp* ptr;
  Matrix mS(JOINT_NO,JOINT_NO), mH(JOINT_NO,JOINT_NO), mQ(2*JOINT_NO,2*JOINT_NO), mR(JOINT_NO,JOINT_NO);
  int k = 0;
  for(int c = 0; c < JOINT_NO; c++) {
    for(int r = 0; r < JOINT_NO; r++,k++) {
      mS(r,c) = S[k];
      mH(r,c) = H[k];
      mR(r,c) = R[k];
    }
  }
  k = 0;
  for(int c = 0; c < 2*JOINT_NO; c++) {
    for(int r = 0; r < 2*JOINT_NO; r++,k++) 
      mQ(r,c) = Q[k];
  }

  if(ind == ADD_NEW) {
    if(_nextIndex == ARRAY_LEN) return ERR_NO_SLOTS;  // can't add new
    // create new
    ptr = new DKalmanObserverExp(&robot,mS,mH,mQ,mR);
    observer[_nextIndex] = ptr;
    return _nextIndex++;    // new element _nextIndex
  } else if(ind < ADD_NEW || ind >= ARRAY_LEN) {
    return ERR_WRONG_INDEX;
  }

  ptr = (DKalmanObserverExp*) observer[ind];
  if(ptr->type() != ID_DKalmanObserverExp)
    return ERR_WRONG_TYPE;
  ptr->settings(mQ,mR);

  return ind;
}


int configFilterDynObserver(int ind, double cutOff, double dt)
{
  FDynObserver* ptr;

  if(ind == ADD_NEW) {
    if(_nextIndex == ARRAY_LEN) return ERR_NO_SLOTS;  // can't add new
    // create new
    ptr = new FDynObserver(&robot,cutOff,dt);
    observer[_nextIndex] = ptr;
    return _nextIndex++;    // new element _nextIndex
  } else if(ind < ADD_NEW || ind >= ARRAY_LEN) {
    return ERR_WRONG_INDEX;
  }

  ptr = (FDynObserver*) observer[ind];
  if(ptr->type() != ID_FDynObserver)
    return ERR_WRONG_TYPE;
  ptr->settings(cutOff,dt);

  return ind;
}

int configFilterRangeObserver(int ind, double cutOff, double dt, double k)
{
  FRangeObserver* ptr;
  
  if(ind == ADD_NEW) {
    if(_nextIndex == ARRAY_LEN) return ERR_NO_SLOTS; 
    ptr = new FRangeObserver(&robot,cutOff,dt,k);
    observer[_nextIndex] = ptr;
    return _nextIndex++;
  } else if(ind < ADD_NEW || ind >= ARRAY_LEN) {
    return ERR_WRONG_INDEX;
  }
  
  ptr = (FRangeObserver*) observer[ind];
  if(ptr->type() != ID_FRangeObserver)
    return ERR_WRONG_TYPE;
  ptr->settings(cutOff,dt,k);
  
  return ind;
}
