cmake_minimum_required(VERSION 2.8)

if(WIN32)
        if(MSVC)
                set(CMAKE_FIND_LIBRARY_PREFIXES "")
                ses(CMAKE_FIND_LIBRARY_SUFFIXES ".lib")
        endif()
else()
        set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
        set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
endif()

find_package( OpenGL REQUIRED)

find_package(Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED 4.5)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_COMPILER /usr/bin/g++)


project(Parameter)


message("${CMAKE_BUILD_TYPE}")
message("${CMAKE_CXX_COMPILER}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

add_subdirectory(src/Graphite)
add_subdirectory(src/Common)
add_subdirectory(src/Numerical)
add_subdirectory(src/ModelMesh)
add_subdirectory(src/MainWindow)
add_subdirectory(src/OpenGL)
add_subdirectory(src/Param)
add_subdirectory(src/UI)

     
        
       
