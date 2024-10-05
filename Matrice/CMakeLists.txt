cmake_minimum_required(VERSION 3.10)

# Set project name and version
project(matrice # project name
    VERSION 1.10.05 # version
    LANGUAGES CXX # languages
    )

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

# Add the library source files
file(GLOB_RECURSE SRC_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cxx"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
)
add_library(matrice STATIC ${SRC_FILES})


# Specify include directories
target_include_directories(matrice
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Set the target properties for the library
set_target_properties(matrice PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 1
)

# Define MATRICE_MATH_KERNEL option
set(MATRICE_MATH_KERNEL "MATRICE_USE_NAT" CACHE STRING "Set math kernel (switch to MATRICE_USE_MKL for enabling MKL support)")

# Check if MATRICE_MATH_KERNEL is set to MATRICE_USE_MKL
if(MATRICE_MATH_KERNEL STREQUAL "MATRICE_USE_MKL")
    find_package(MKL REQUIRED)
    if(MKL_FOUND)
        message(STATUS "MATRICE_MATH_KERNEL set to MATRICE_USE_MKL. Enabling MKL support.")
        add_definitions(-DMATRICE_MATH_KERNEL)  # Define for C++ preprocessor
        target_link_libraries(matrice PUBLIC MKL::MKL)  # Link MKL
    else()
        message(WARNING "MKL not found, but MATRICE_MATH_KERNEL is set to MATRICE_USE_MKL.")
    endif()
else()
    message(STATUS "MATRICE_MATH_KERNEL is not set to MATRICE_USE_MKL. MKL disabled.")
endif()

# Install the library and headers
install(TARGETS matrice
    EXPORT MatriceTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

# Install the headers
install(DIRECTORY include/ DESTINATION include)

# Export the library for use in other CMake projects
install(EXPORT MatriceTargets
    FILE MatriceConfig.cmake
    NAMESPACE matrice::
    DESTINATION lib/cmake/Matrice
)