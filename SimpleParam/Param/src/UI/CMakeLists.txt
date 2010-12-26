
include_directories( $ENV{QTDIR}/include
                     $ENV{QTDIR}/include/QtOpenGL
                   )

file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(ui STATIC ${HEADERS} ${SOURCES})
