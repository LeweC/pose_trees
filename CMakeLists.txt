cmake_minimum_required(VERSION 3.8)
project(pose_trees)

set(CMAKE_CXX_STANDARD 20)

find_package(ament_cmake REQUIRED)
find_package(TBB REQUIRED)
find_package(Eigen3 REQUIRED)
list(APPEND CMAKE_MODULE_PATH "/usr/include/tbb/")

find_package(OpenMP)
if(OPENMP_FOUND)
  if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
  endif()
endif()

include_directories(include)

set(SOURCES src/pose_test_3d.cpp)
set(SOURCES src/pose_test_2d.cpp)
set(SOURCES src/pose_test_1d.cpp)

add_library(${PROJECT_NAME} SHARED ${SOURCES})

target_link_libraries(${PROJECT_NAME} TBB::tbb OpenMP::OpenMP_CXX)

ament_target_dependencies(${PROJECT_NAME}
    ament_cmake
)

add_executable(PoseTest3D src/pose_test_3d.cpp)
target_link_libraries(PoseTest3D Eigen3::Eigen)
target_link_libraries(PoseTest3D
    ${PROJECT_NAME}
)

install(TARGETS PoseTest3D
        DESTINATION lib/${PROJECT_NAME})

set_target_properties(PoseTest3D PROPERTIES LINKER_LANGUAGE CXX)

add_executable(PoseTest2D src/pose_test_2d.cpp)
target_link_libraries(PoseTest2D Eigen3::Eigen)
target_link_libraries(PoseTest2D
    ${PROJECT_NAME}
)
install(TARGETS PoseTest2D
        DESTINATION lib/${PROJECT_NAME})

set_target_properties(PoseTest2D PROPERTIES LINKER_LANGUAGE CXX)


add_executable(PoseTest1D src/pose_test_1d.cpp)
target_link_libraries(PoseTest1D Eigen3::Eigen)
target_link_libraries(PoseTest1D
    ${PROJECT_NAME}
)
install(TARGETS PoseTest1D
        DESTINATION lib/${PROJECT_NAME})

set_target_properties(PoseTest1D PROPERTIES LINKER_LANGUAGE CXX)

install(
  DIRECTORY include/
  DESTINATION include
)

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include)

ament_export_dependencies(TBB)
ament_export_libraries(${PROJECT_NAME})
ament_export_include_directories(include)
ament_package()
