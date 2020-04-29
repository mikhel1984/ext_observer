// Test observers with UR10 expreimental results

#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector>
#include <string> 
#include <cstdlib>
#include <iomanip>

#include "observers.h"


// settings
#define OUT_NAME "external.csv"

#define BEG_Q    1
#define BEG_QD   7
#define BEG_TAU  13
#define BEG_Q2D  19

// Reag csv line of numbers
bool csvRead(std::istream& str, std::vector<double>& v);


int main(int argc, char** argv)
{  
  int N = 6;
  double q[N], qd[N], tau[N], ext[N];

  // torque coefficients
  double K[] = {10.0,10.6956,8.4566,9.0029,9.48,10.1232};

  double S[6*6] = {0};
  double H[] = {1,0,0,0,0,0,
                0,1,0,0,0,0,
                0,0,1,0,0,0,
                0,0,0,1,0,0,
                0,0,0,0,1,0,
                0,0,0,0,0,1}; 
  double Q[] = {0.2,0,0,0,0,0,
                0,0.2,0,0,0,0,
                0,0,0.2,0,0,0,
                0,0,0,30,0,0,
                0,0,0,0,30,0,
                0,0,0,0,0,30};
  double R[] = {5E-4,0,0,0,0,0,
                0,5E-4,0,0,0,0,
                0,0,5E-4,0,0,0,
                0,0,0,5E-4,0,0,
                0,0,0,0,5E-4,0,
                0,0,0,0,0,5E-4};
  int id_kfe = configDistKalmanObserverExp(-1,S,H,Q,R);
  
  
  // save to file 
  std::ofstream oFile(OUT_NAME);
  oFile.setf(std::ios::fixed);
  oFile.precision(4);

  // read from file
  std::ifstream iFile("link4.csv");

  std::vector<double> val;
  val.reserve(N);

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
    curr -= zero;  // plot from 0

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
