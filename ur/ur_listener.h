#ifndef UR_LISTENER_H
#define UR_LISTENER_H

#include <sys/socket.h> // socket
//#include <sys/types.h>
#include <arpa/inet.h>  // inet_addr, htonl
#include <unistd.h>     // usleep

#include <eigen3/Eigen/Geometry>

#include <iostream>

#define TCP_BUF_SIZE (1024*5)
#define SLEEP_USEC   1000
#define COUNT_MAX    1000

#define CHAR_TIME    4
#define CHAR_Q       (4+8+5*48)
#define UR_JOINTS    6

class UrListener {
public:
  UrListener(const char*, unsigned short p);
  ~UrListener();
  
  bool isConnected() { return connected; }

  void listen(void (*fn)(double,
                 Eigen::VectorXd&,Eigen::VectorXd&,Eigen::VectorXd&));

private:
  double swap(double val);

  struct sockaddr_in sender;
  char* buf;
  int sock;
  bool connected;

}; // UrListener

UrListener::UrListener(const char* addr, unsigned short port)
  : buf(0)
  , connected(false)
{
  // create socket
  if((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
    std::cout << "Can't create socket" << std::endl;
    return;
  }
  // define address
  if(inet_addr(addr) == INADDR_NONE) {
    std::cout << "Wrong address " << addr << std::endl;
    return;
  }
  sender.sin_addr.s_addr = inet_addr(addr);
  sender.sin_family = AF_INET;
  sender.sin_port = htons(port);
  sender.sin_zero[0] = '\0';
  // connect
  if(connect(sock, (struct sockaddr*) &sender, sizeof(sender)) < 0) {
    std::cout << "Connection failed!" << std::endl;
    return;
  }
  connected = true;
  std::cout << "Connected" << std::endl;

  buf = new char[TCP_BUF_SIZE];
}

UrListener::~UrListener()
{
  //if(sock > -1) close(sock);
  if(buf) delete[] buf;
}

void UrListener::listen(void (*fn)(double,
                 Eigen::VectorXd&,Eigen::VectorXd&,Eigen::VectorXd&))
{
  Eigen::VectorXd vq(UR_JOINTS), vqd(UR_JOINTS), vi(UR_JOINTS);
  
  int bytes, i;
  unsigned int count = 0;
  double *arr, tm;
  while(count < COUNT_MAX) {
    bytes = recv(sock, buf, TCP_BUF_SIZE, MSG_DONTWAIT);
    if(bytes > 0) {
      // time 
      arr = (double*) (buf + CHAR_TIME);
      tm = swap(*arr);
      // angles
      arr = (double*) (buf + CHAR_Q);
      for(i = 0; i < UR_JOINTS; i++) vq(i) = swap(arr[i]);
      // velocities
      arr += UR_JOINTS;
      for(i = 0; i < UR_JOINTS; i++) vqd(i) = swap(arr[i]);
      // currents
      arr += UR_JOINTS;
      for(i = 0; i < UR_JOINTS; i++) vi(i) = swap(arr[i]);
      // call function
      if(fn) 
        fn(tm,vq,vqd,vi);
      count = 0;
    } else {
      usleep(SLEEP_USEC);
      ++count;
    }
  }
  if(count >= COUNT_MAX) 
    std::cout << "No data" << std::endl;
}

double UrListener::swap(double val)
{
  uint32_t *src, *dst;
  double res;

  src = (uint32_t*) &val;
  dst = (uint32_t*) &res;
  dst[0] = htonl(src[1]);
  dst[1] = htonl(src[0]);

  return res;
}

#endif // UR_LISTENER_H
