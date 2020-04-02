cmake_minimum_required(VERSION 2.8.3)
project(double_link_observer)

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
# include
  ${catkin_INCLUDE_DIRS}
)

add_executable(observer src/main.cpp)
#add_dependencies(observer double_link_observer_generate_messages_cpp)
target_link_libraries(observer ${catkin_LIBRARIES})
