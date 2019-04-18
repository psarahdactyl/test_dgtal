PROJECT(cc)
cmake_minimum_required(VERSION 3.1)
set(ROOT "${CMAKE_CURRENT_LIST_DIR}")
message(${ROOT})
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED)

# DGtal and Boost
find_package(DGtal REQUIRED)
find_package(Boost 1.40 COMPONENTS program_options REQUIRED)

include_directories(${DGTAL_INCLUDE_DIRS})
include_directories(${Boost_INCLUDE_DIR})

link_directories(${DGTAL_LIBRARY_DIRS})

# executable
add_executable(${PROJECT_NAME} shortcuts.cpp)
#target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${OPENGL_INCLUDE_DIR})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/eigen/)
#target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/glfw/include)
#target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/glad/include)
target_link_libraries(${PROJECT_NAME} ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} igl::opengl igl::opengl_glfw igl::core)