CMAKE_MINIMUM_REQUIRED(VERSION 3.24.0)

PROJECT(small_mcftle_gpu CXX)

SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

IF(WIN32)
  ADD_DEFINITIONS(-DNOMINMAX)
ENDIF()

# target for the compute shader
set(SOURCES mcftle.hlsl flow.hlsli params.hlsli random.hlsli)
ADD_CUSTOM_TARGET(resources ALL SOURCES ${SOURCES})
ADD_CUSTOM_COMMAND(OUTPUT mcftle.cso
                     COMMAND fxc.exe /Fo ${CMAKE_BINARY_DIR}/mcftle.cso ${CMAKE_CURRENT_SOURCE_DIR}/mcftle.hlsl /nologo /O2 /Tcs_5_0
                     MAIN_DEPENDENCY mcftle.hlsl
                     COMMENT "Effect-compile mcftle.hlsl \n ${COMMAND}"
                     VERBATIM)

# executable
set(SOURCES main.cpp camera.hpp cbuffer.hpp d3d.hpp frame.hpp math.hpp timer.hpp)
ADD_EXECUTABLE(small_mcftle_gpu ${SOURCES})
TARGET_LINK_LIBRARIES(small_mcftle_gpu d3d11.lib dxgi.lib)
ADD_DEPENDENCIES(small_mcftle_gpu resources)