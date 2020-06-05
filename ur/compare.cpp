// Test observers with UR10 expreimental results

#include <iostream> 
#include <fstream>
#include <vector>
#include <string> 
#include <cstdlib>
#include <iomanip>
#include "ur_robot.h"
#include <ctime>


#include "../lib/momentum_observer.h" 
#include "../lib/disturbance_observer.h"
#include "../lib/sliding_mode_observer.h"
#include "../lib/disturbance_kalman_filter.h"
#include "../lib/filtered_dyn_observer.h"


#define BEG_Q    1
#define BEG_QD   7
#define BEG_TAU  13
#define BEG_Q2D  19

// Reag csv line of numbers
bool csvRead(std::istream& str, std::vector<double>& v);

Vector acceleration(Vector& qd3, double dt);

int main(int argc, char** argv)
{  
  if(argc != 2) {
    std::cout << "Usage: ./test file.csv" << std::endl;
    return 0;
  }
  UrRobot robot; 
  int N = robot.jointNo();
  Vector q(N), qd(N), tau(N), ext(N);

  // torque coefficients
  double K[] = {10.0,10.6956,8.4566,9.0029,9.48,10.1232};

  Vector k(N);
  k << 90,50,50,90,90,40; 
  MomentumObserver m_observer(&robot,k);

  double sigma = 21, xeta = 18, beta = 50;
  DisturbanceObserver d_observer(&robot,sigma,xeta,beta); 

  Vector T1(N), S1(N), T2(N), S2(N);
  S1 << 20,30,20,30,20,30; 
  for(int i = 0; i < N; i++) T1(i) = 2*sqrt(S1(i));
  S2 << 10,10,10,10,10,10; // set 0 to exclude linear part
  for(int i = 0; i < N; i++) T2(i) = 2*sqrt(S2(i));
  SlidingModeObserver sm_observer(&robot,T1,S1,T2,S2);

  Matrix S = Matrix::Zero(N,N);
  Matrix H = Matrix::Identity(N,N);
  Matrix Q = Matrix::Identity(2*N,2*N);
  for(int i = 0; i < N; i++) Q(i,i) = 0.002;
  for(int i = N; i < 2*N; i++) Q(i,i) = 0.3;
  Matrix R = Matrix::Identity(N,N);  
  R *= 0.05;
  DKalmanObserver dkm_observer(&robot,S,H,Q,R);

  FDynObserver fd_observer(&robot, 8, 1/100.0);  // assume frequency is 100 Hz

  // read from file
  std::ifstream iFile(argv[1]);
  
  std::vector<double> val;
  val.reserve(N);

  double prev = -1, dt = 0;

  std::vector<Vector> qs, qds, taus;
  std::vector<double> steps;

  while(csvRead(iFile,val)) {

    double curr = val[0];
    if(prev > 0) {
      dt = curr - prev;
    } else {
      dt = 0;
    }
    prev = curr;

    // read values
    for(int i = 0; i < N; i++) {
      q(i)   = val[BEG_Q + i];
      qd(i)  = val[BEG_QD + i];
      tau(i) = val[BEG_TAU + i] * K[i];
    }
    qs.push_back(q);
    qds.push_back(qd);
    taus.push_back(tau);
    steps.push_back(dt);
  }
  iFile.close();

  int nSum = qs.size(); 
  std::cout << "Samples: " << nSum << std::endl; 

  clock_t tBeg; 

  tBeg = clock();
  for(int i = 0; i < nSum; i++) {
    ext = m_observer.getExternalTorque(qs[i],qds[i],taus[i],steps[i]);
  }
  double timeSec = double(clock() - tBeg) / CLOCKS_PER_SEC;
  std::cout << "Momentum observer, total (s): " << timeSec << ", average (ms): " << (timeSec * 1000/ nSum) << std::endl;

  tBeg = clock();
  for(int i = 0; i < nSum; i++) {
    ext = d_observer.getExternalTorque(qs[i],qds[i],taus[i],steps[i]);
  }
  timeSec = double(clock() - tBeg) / CLOCKS_PER_SEC;
  std::cout << "Disturbance observer, total (s): " << timeSec << ", average (ms): " << (1000* timeSec / nSum) << std::endl;

  tBeg = clock();
  for(int i = 0; i < nSum; i++) {
    ext = sm_observer.getExternalTorque(qs[i],qds[i],taus[i],steps[i]);
  }
  timeSec = double(clock() - tBeg) / CLOCKS_PER_SEC;
  std::cout << "Sliding mode observer, total (s): " << timeSec << ", average (ms): " << (timeSec * 1000/ nSum) << std::endl;
    
  tBeg = clock();
  for(int i = 0; i < nSum; i++) {
    ext = dkm_observer.getExternalTorque(qs[i],qds[i],taus[i],steps[i]);
  }
  timeSec = double(clock() - tBeg) / CLOCKS_PER_SEC;
  std::cout << "Kalman filter observer, total (s): " << timeSec << ", average (ms): " << (timeSec * 1000 / nSum) << std::endl;

  tBeg = clock();
  for(int i = 0; i < nSum; i++) {
    ext = fd_observer.getExternalTorque(qs[i],qds[i],taus[i],steps[i]);
  }
  timeSec = double(clock() - tBeg) / CLOCKS_PER_SEC;
  std::cout << "Filtered dynamics observer, total (s): " << timeSec << ", average (ms): " << (timeSec * 1000 / nSum) << std::endl;
      
  
  return 0;
}

bool csvRead(std::istream& str, std::vector<double>& result)
{
  std::string                line;
  std::getline(str,line);
  std::stringstream          lineStream(line);
  std::string                cell;

  result.clear();

  while(std::getline(lineStream,cell, ','))  {
    result.push_back((double) atof(cell.c_str()) );
  }

  return !result.empty();
}

Vector acceleration(Vector& qd3, double dt)
{
  static int count = 0;
  static Vector q2d(6);
  static Vector qd1(6), qd2(6);

  if(count == 0) {
    q2d.setZero();
  } 
  else if(count == 1) {
    q2d = qd3 - qd2;
    q2d /= dt;
  } 
  else {
    q2d = qd1 - 4*qd2 + 3*qd3;
    q2d /= 2*dt;
  }
  qd1 = qd2;
  qd2 = qd3;
  count++;
  return q2d;
}
