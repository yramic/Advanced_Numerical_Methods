cmake_minimum_required(VERSION 2.8)
#project(quad)
add_subdirectory(../../../../CppHilbert/Library ${CMAKE_SOURCE_DIR}/build/bin)

include_directories(../../../../CppHilbert/Library/source/)

add_library(quadrature_internal STATIC integral_gauss.cpp logweight_quadrature.cpp genLaguerreRule.cpp)

#target_link_libraries(quadrature CppHilbert)
#target_link_libraries(quadrature )

add_custom_target(quadrature ALL
  COMMAND ${CMAKE_AR} rc libquadrature.a $<TARGET_FILE:quadrature_internal> $<TARGET_FILE:CppHilbert>)
