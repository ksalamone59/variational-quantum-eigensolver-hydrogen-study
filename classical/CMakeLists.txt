cmake_minimum_required(VERSION 3.12)
project(classical_eigensolver)
set(CMAKE_CXX_STANDARD 17)

add_executable(classical_eigensolver src/classical_eigensolver.cpp)

include_directories(
        ${ONNXRUNTIME_DIR}/build/native/include
        ${CMAKE_SOURCE_DIR}/include
        ${EIGEN3_INCLUDE_DIRS}
    )