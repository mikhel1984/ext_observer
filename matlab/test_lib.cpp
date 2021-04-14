// Test observers with expreimental results

#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <string> 
#include <cstdlib>
#include <iomanip>

#include "observers.h"

// settings
#define IN_NAME  "input.csv"       // CSV: time, positions, velocities, torques
#define OUT_NAME "external.csv"    // CSV: time, external torques 

#define BEG_Q    1
#define BEG_QD   (BEG_Q+JOINT_NO)
#define BEG_TAU  (BEG_QD+JOINT_NO)
#define BEG_Q2D  (BEG_TAU+JOINT_NO)

// Reag csv line of numbers
bool csvRead(std::istream& str, std::vector<double>& v);


int main(int argc, char** argv)
{  
  int N = JOINT_NO;
  double q[JOINT_NO], qd[JOINT_NO], tau[JOINT_NO], ext[JOINT_NO];

  // torque coefficients
  double K[] = {10.0,10.6956};

  double S[JOINT_NO*JOINT_NO] = {0};
  double H[] = {1,0,
                0,1};   
  double Q[(2*JOINT_NO)*(2*JOINT_NO)] = {0};
  for(int i = 0; i < (2*JOINT_NO)*(2*JOINT_NO); i += JOINT_NO+1) {
    Q[i] = (i < 2*JOINT_NO*JOINT_NO) ? 0.2 : 30;    
  }
  double R[] = {5E-4,0,
                0,5E-4};
  int id_kfe = configDistKalmanObserverExp(-1,S,H,Q,R);
  
  
  // save to file 
  std::ofstream oFile(OUT_NAME);
  oFile.setf(std::ios::fixed);
  oFile.precision(4);

  // read from file
  std::ifstream iFile(IN_NAME);

  std::vector<double> val;
  val.reserve(JOINT_NO);

  double prev = -1, zero = 0, dt = 0;
  // range
  //double tBeg = -1, tEnd = 1E100;
  
  while(csvRead(iFile,val)) {

    double curr = val[0];
    if(prev > 0) {
      dt = curr - prev;
    } else {
      dt = 0;
      zero = curr;
    }
    prev = curr;
    curr -= zero;  // start from 0

    //if(curr < tBeg || curr > tEnd) continue;

    // read values
    for(int i = 0; i < N; i++) {
      q[i]   = val[BEG_Q + i];
      qd[i]  = val[BEG_QD + i];
      tau[i] = val[BEG_TAU + i] * K[i];
      //q2d(i) = val[BEG_Q2D + i];
    }
    
    getExternalTorque(id_kfe,ext,q,qd,tau,dt);

    // save result 
    oFile <<  curr;
    for(int i = 0; i < N; i++) {
      oFile << "," << ext[i];
    }
    oFile << std::endl;
  }
      
  oFile.close();
  iFile.close();
  
  freeAll();
  
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
