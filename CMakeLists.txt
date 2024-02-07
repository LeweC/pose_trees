cmake_minimum_required(VERSION 3.8)
project(pose_trees)

set(CMAKE_CXX_STANDARD 20)


find_package(Eigen3 REQUIRED)


include_directories(include)

set(SOURCES src/pose_test_3d.cpp)

add_library(${PROJECT_NAME} SHARED ${SOURCES})
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)

add_executable(PoseTest3D src/pose_test_3d.cpp)
target_link_libraries(PoseTest3D Eigen3::Eigen)
target_link_libraries(PoseTest3D
    ${PROJECT_NAME}
)

install(TARGETS PoseTest3D
        DESTINATION lib/${PROJECT_NAME})

set_target_properties(PoseTest3D PROPERTIES LINKER_LANGUAGE CXX)

install(
  DIRECTORY include/
  DESTINATION include
)

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)