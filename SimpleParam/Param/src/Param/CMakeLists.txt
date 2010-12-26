include_directories( ${PROJECT_SOURCE_DIR}/include )
include_directories( ${PROJECT_SOURCE_DIR}/include/hj_3rd)

set ( HEADERS Parameterization.h
              Barycentric.h
              ParamPatch.h
              ParamChart.h
              ChartCreator.h
              TransFunctor.h
              ParamDrawer.h
              Parameter.h )

set ( SOURCES Parameterization.cc
              Barycentric.cc
              TransFunctor.cc
              ChartCreator.cc
              ParamDrawer.cc
              Parameter.cc )
              

                  

add_library(param STATIC ${HEADERS} ${SOURCES})
