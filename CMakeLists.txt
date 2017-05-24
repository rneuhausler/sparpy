cmake_minimum_required(VERSION 2.8)

# Boost
find_package(Boost 1.50.0 COMPONENTS python REQUIRED)
list(APPEND sparpy_LIBRARIES ${Boost_LIBRARIES})
list(APPEND sparpy_INCLUDES ${Boost_INCLUDE_DIRS})

# VTK 
find_package(VTK REQUIRED)
if (VTK_FOUND)
    add_definitions(-DHAVE_VTK)
endif(VTK_FOUND)
list(APPEND sparpy_LIBRARIES ${VTK_LIBRARIES})
list(APPEND sparpy_INCLUDES ${VTK_INCLUDE_DIRS})

# Python
find_package(PythonLibs REQUIRED)
list(APPEND sparpy_LIBRARIES ${PYTHON_LIBRARIES})
list(APPEND sparpy_INCLUDES ${PYTHON_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wno-deprecated -std=c++14")
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

# Aboria
list(APPEND sparpy_INCLUDES Aboria/src)
list(APPEND sparpy_INCLUDES Aboria/third-party)

include_directories(src ${sparpy_INCLUDES})

set(sparpy_SOURCE
    src/python.cpp
    )

set(sparpy_HEADERS
    src/python.hpp
    src/sparpy.hpp
    src/simulation.hpp
    src/interactions.hpp
    src/timestepping.hpp
    )

PYTHON_ADD_MODULE (sparpy ${sparpy_SOURCE})
target_link_libraries(sparpy ${sparpy_LIBRARIES})

add_subdirectory(tests)
