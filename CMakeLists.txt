PROJECT(mesh)
#Required in DGtal
cmake_minimum_required(VERSION 3.1)
set(ROOT "${CMAKE_CURRENT_LIST_DIR}")
message(${ROOT})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# create executable
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

if (NOT CMAKE_VERSION VERSION_LESS "3.11")
  cmake_policy(SET CMP0072 NEW)
endif()

find_package(DGtal REQUIRED)
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)
find_package(LIBIGL REQUIRED QUIET)

include_directories(${DGTAL_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})

link_directories(${DGTAL_LIBRARY_DIRS})

add_executable(${PROJECT_NAME} shortcuts.cpp)

target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/eigen/)
target_link_libraries(${PROJECT_NAME} ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} igl::core)
