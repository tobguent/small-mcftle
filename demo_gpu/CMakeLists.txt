CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)

PROJECT(small_mcftle_gpu CXX)

SET (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)
FIND_PACKAGE(DirectX REQUIRED)

IF(WIN32)
  ADD_DEFINITIONS(-DNOMINMAX)
ENDIF()

# target for the compute shader
set(SOURCES mcftle.hlsl flow.hlsli params.hlsli random.hlsli)
ADD_CUSTOM_TARGET(resources ALL SOURCES ${SOURCES})
ADD_CUSTOM_COMMAND(OUTPUT mcftle.cso
                     COMMAND ${DirectX_FXC_EXECUTABLE} /Fo ${CMAKE_BINARY_DIR}/mcftle.cso ${CMAKE_CURRENT_SOURCE_DIR}/mcftle.hlsl /nologo /O2 /Tcs_5_0
                     MAIN_DEPENDENCY mcftle.hlsl
                     COMMENT "Effect-compile mcftle.hlsl \n ${COMMAND}"
                     VERBATIM)

# executable
set(SOURCES main.cpp camera.hpp cbuffer.hpp d3d.hpp frame.hpp math.hpp timer.hpp)
ADD_EXECUTABLE(small_mcftle_gpu ${SOURCES})
TARGET_LINK_LIBRARIES(small_mcftle_gpu ${DirectX_D3D11_LIBRARY} ${DirectX_DXGI_LIBRARY})
ADD_DEPENDENCIES(small_mcftle_gpu resources)