PROJECT(cc)
#Required in DGtal
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

set(ROOT "${CMAKE_CURRENT_LIST_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

FIND_PACKAGE(DGtal REQUIRED)
FIND_PACKAGE( Boost 1.40 COMPONENTS program_options REQUIRED )
find_package(LIBIGL REQUIRED QUIET)
find_package(Qt5Core)
find_package(Qt5Widgets)
find_package(Qt53DCore)
find_package(Qt53DExtras)
find_package(Qt53DInput)
find_package(Qt53DLogic)
find_package(Qt53DRender)
#find_package(QGLViewer)
#set(QGLVIEWER_INCLUDE_DIR /usr/local/lib/QGLViewer.framework/Headers)
#set(QGLVIEWER_LIBRARIES /usr/local/lib/QGLViewer.framework)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
#include_directories( ${QGLVIEWER_INCLUDE_DIR} )   
include_directories(${Qt5Core_INCLUDE_DIRS} 
   ${Qt5Widgets_INCLUDE_DIRS} 
   ${Qt53DCore_INCLUDE_DIRS}
   ${Qt53DExtras_INCLUDE_DIRS}
   ${Qt53DInput_INCLUDE_DIRS}
   ${Qt53DLogic_INCLUDE_DIRS}
   ${Qt53DRender_INCLUDE_DIRS})

set(UI_LIBS  Qt5::Core Qt5::Widgets Qt5::3DCore Qt5::3DExtras 
   Qt5::3DExtras Qt5::3DInput Qt5::3DLogic
   Qt5::3DRender CACHE INTERNAL "" FORCE)

LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

#ADD_EXECUTABLE(helloworld helloworld)
#ADD_EXECUTABLE(vol_viewer vol_viewer)
#ADD_EXECUTABLE(connected_comp connected_comp.cpp)
ADD_EXECUTABLE(shortcuts shortcuts.cpp)

#TARGET_LINK_LIBRARIES(helloworld ${DGTAL_LIBRARIES})
#TARGET_LINK_LIBRARIES(vol_viewer ${DGTAL_LIBRARIES} ${Boost_LIBRARIES})
#target_link_libraries(shortcuts ${UI_LIBS})
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${ROOT}/eigen/)
TARGET_LINK_LIBRARIES(shortcuts ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} igl::core)
