cmake_minimum_required(VERSION 3.10.0)
project(FMMPotential)

set(CMAKE_DEBUG_POSTFIX -debug)
set(CMAKE_RELEASE_POSTFIX -release)
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

find_package(OpenMP)
find_package(Eigen3 3.4.0 REQUIRED)
find_package(MishMesh)
find_package(OpenMesh REQUIRED)

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

add_library(FMMPotential SHARED src/FMMPotential.cpp src/FMMPotential_impl.h src/quadrature.hpp include/FMMPotential/FMMPotential.h)
add_library(libFMMPotential STATIC src/FMMPotential.cpp src/FMMPotential_impl.h src/quadrature.hpp include/FMMPotential/FMMPotential.h)

target_include_directories(FMMPotential PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>
)

target_include_directories(libFMMPotential PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>
)
target_link_libraries(FMMPotential OpenMeshCore Eigen3::Eigen)

set(CLANG_LIBOMP "CLANG_LIBOMP-NOTFOUND" CACHE FILEPATH "Path the clang libomp.lib")
if(${CLANG_LIBOMP})
    target_link_libraries(FMMPotential ${CLANG_LIBOMP})
    target_link_libraries(libFMMPotential ${CLANG_LIBOMP})
endif()

if(TARGET MishMesh)
        target_link_libraries(FMMPotential MishMesh)
        target_link_libraries(libFMMPotential MishMesh)
        target_compile_definitions(FMMPotential PUBLIC WITH_MISHMESH)
        target_compile_definitions(libFMMPotential PUBLIC WITH_MISHMESH)
endif(TARGET MishMesh)


set(HEADER_DIRS ${CMAKE_CURRENT_BINARY_DIR}/include)
string(TOUPPER "${PROJECT_NAME}" PROJECT_NAME_UPPERCASE)
configure_file(FMMPotentialConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/FMMPotentialConfig.cmake @ONLY)

install(TARGETS FMMPotential libFMMPotential EXPORT FMMPotential-targets LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(EXPORT FMMPotential-targets FILE FMMPotentialTargets.cmake DESTINATION ${CMAKE_BINARY_DIR})
install(EXPORT FMMPotential-targets FILE FMMPotentialTargets.cmake DESTINATION ${CMAKE_INSTALL_PREFIX})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/FMMPotentialConfig.cmake DESTINATION ${CMAKE_INSTALL_PREFIX})
