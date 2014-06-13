cmake_minimum_required (VERSION 2.6)
project (ContourRelaxedSuperpixels)

set (CMAKE_BUILD_TYPE Release)

if (CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -ansi")
endif()

find_package (OpenCV REQUIRED)
find_package (Boost COMPONENTS filesystem REQUIRED)
find_package (Boost COMPONENTS system REQUIRED)

add_executable (contourRelaxedSuperpixels contourRelaxedSuperpixels.cpp)
target_link_libraries (contourRelaxedSuperpixels ${OpenCV_LIBS} ${Boost_LIBRARIES})

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)
