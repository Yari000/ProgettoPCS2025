# message("SUBFOLDER:")

# message(${CMAKE_SOURCE_DIR})
# message(${CMAKE_CURRENT_SOURCE_DIR})
# message(${CMAKE_BINARY_DIR})
# message(${CMAKE_CURRENT_BINARY_DIR})
# message(${PROJECT_NAME})
# message(${CMAKE_PROJECT_NAME})

list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/PolygonalMesh.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Triangle.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/GeometryUtils.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/GeodesicGenerator.hpp)
list(APPEND src_headers ${CMAKE_CURRENT_SOURCE_DIR}/ShortestPath.hpp)

list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/GeometryUtils.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/GeodesicGenerator.cpp)
list(APPEND src_sources ${CMAKE_CURRENT_SOURCE_DIR}/ShortestPath.cpp)

list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_headers ${src_headers} PARENT_SCOPE)
set(src_sources ${src_sources} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)
