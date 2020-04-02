cmake_minimum_required(VERSION 2.8.3)
project(ur_observer)

## Compile as C++11, supported in ROS Kinetic and newer
# add_compile_options(-std=c++11)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  message_generation
)

## Generate messages in the 'msg' folder
add_message_files(
   DIRECTORY msg
)

## Generate added messages and services with any dependencies listed here
generate_messages(DEPENDENCIES std_msgs)


catkin_package(
  CATKIN_DEPENDS message_runtime
)

include_directories(
 include
  ${catkin_INCLUDE_DIRS}
)

set(UR_SOURCES
  # observer & robot
  src/ur_robot.cpp
  # dynamics
  src/M_mtrx_fcn/M_mtrx_fcn.c
  src/M_mtrx_fcn/M_mtrx_fcn_initialize.c
  src/M_mtrx_fcn/M_mtrx_fcn_terminate.c
  src/C_mtrx_fcn/C_mtrx_fcn.c
  src/C_mtrx_fcn/C_mtrx_fcn_initialize.c
  src/C_mtrx_fcn/C_mtrx_fcn_terminate.c
  src/G_vctr_fcn/G_vctr_fcn.c
  src/G_vctr_fcn/G_vctr_fcn_initialize.c
  src/G_vctr_fcn/G_vctr_fcn_terminate.c
  src/F_vctr_fcn/F_vctr_fcn.c
  src/F_vctr_fcn/F_vctr_fcn_initialize.c
  src/F_vctr_fcn/F_vctr_fcn_terminate.c
)

add_executable(ur_observer src/main.cpp ${UR_SOURCES}) 
target_link_libraries(ur_observer ${catkin_LIBRARIES}) 