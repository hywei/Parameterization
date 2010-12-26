file(GLOB HEADERS *.h)
file(GLOB SOURCES *.cpp)

add_library(meshmodel STATIC ${HEADERS} ${SOURCES})
