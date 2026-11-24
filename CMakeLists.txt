cmake_minimum_required(VERSION 3.16)
project(ParallelAnkle VERSION 0.1.0 LANGUAGES CXX)

set(parallel_ankle_top_level OFF)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(parallel_ankle_top_level ON)
endif()
option(PARALLEL_ANKLE_WITH_YAML "Build the YAML profile loader" ON)
option(PARALLEL_ANKLE_BUILD_EXAMPLES "Build examples" ${parallel_ankle_top_level})
option(PARALLEL_ANKLE_BUILD_TESTS "Build tests with an installed GoogleTest" ${parallel_ankle_top_level})
option(PARALLEL_ANKLE_BUILD_BENCHMARKS "Build bounded ankle kinematics microbenchmark" OFF)

include(GNUInstallDirs)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_library(parallel_ankle src/parameters.cpp src/parallel_ankle.cpp src/five_bar.cpp src/constraints.cpp
    src/transmissions.cpp src/second_order.cpp)
add_library(ParallelAnkle::core ALIAS parallel_ankle)
set_target_properties(parallel_ankle PROPERTIES EXPORT_NAME core POSITION_INDEPENDENT_CODE ON)
target_compile_features(parallel_ankle PUBLIC cxx_std_17)
target_include_directories(parallel_ankle PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(parallel_ankle PUBLIC Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(parallel_ankle PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(PARALLEL_ANKLE_WITH_YAML)
    find_package(yaml-cpp REQUIRED)
    add_library(parallel_ankle_config src/config.cpp)
    add_library(ParallelAnkle::config ALIAS parallel_ankle_config)
    set_target_properties(parallel_ankle_config PROPERTIES EXPORT_NAME config POSITION_INDEPENDENT_CODE ON)
    if(TARGET yaml-cpp::yaml-cpp)
        set(parallel_ankle_yaml_target yaml-cpp::yaml-cpp)
    else()
        set(parallel_ankle_yaml_target yaml-cpp)
    endif()
    target_link_libraries(parallel_ankle_config PUBLIC ParallelAnkle::core PRIVATE ${parallel_ankle_yaml_target})
    if(PARALLEL_ANKLE_BUILD_EXAMPLES)
        add_executable(ankle_demo examples/ankle_demo.cpp)
        target_link_libraries(ankle_demo PRIVATE ParallelAnkle::config)
    endif()
endif()

if(PARALLEL_ANKLE_BUILD_EXAMPLES)
    add_executable(five_bar_demo examples/five_bar_demo.cpp)
    target_link_libraries(five_bar_demo PRIVATE ParallelAnkle::core)
    add_executable(constraints_demo examples/constraints_demo.cpp)
    target_link_libraries(constraints_demo PRIVATE ParallelAnkle::core)
    add_executable(transmissions_demo examples/transmissions_demo.cpp)
    target_link_libraries(transmissions_demo PRIVATE ParallelAnkle::core)
    add_executable(loaded_stiffness_demo examples/loaded_stiffness_demo.cpp)
    target_link_libraries(loaded_stiffness_demo PRIVATE ParallelAnkle::core)
    add_executable(linear_universal_ankle_demo examples/linear_universal_ankle_demo.cpp)
    target_link_libraries(linear_universal_ankle_demo PRIVATE ParallelAnkle::core)
endif()

if(PARALLEL_ANKLE_BUILD_BENCHMARKS)
    add_executable(benchmark_jacobians EXCLUDE_FROM_ALL benchmarks/benchmark_jacobians.cpp)
    target_link_libraries(benchmark_jacobians PRIVATE ParallelAnkle::core)
endif()

if(PARALLEL_ANKLE_BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
    if(TARGET GTest::gtest_main)
        set(parallel_ankle_gtest_main GTest::gtest_main)
    else()
        set(parallel_ankle_gtest_main GTest::Main)
    endif()
    add_executable(parallel_ankle_tests tests/test_solver.cpp tests/test_edge_cases.cpp tests/test_geometry.cpp
        tests/test_five_bar.cpp tests/test_constraints.cpp tests/test_transmissions.cpp tests/test_second_order.cpp)
    target_link_libraries(parallel_ankle_tests PRIVATE ParallelAnkle::core ${parallel_ankle_gtest_main})
    if(PARALLEL_ANKLE_WITH_YAML)
        target_sources(parallel_ankle_tests PRIVATE tests/test_config.cpp)
        target_link_libraries(parallel_ankle_tests PRIVATE ParallelAnkle::config)
        target_compile_definitions(parallel_ankle_tests PRIVATE PARALLEL_ANKLE_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    endif()
    include(GoogleTest)
    gtest_discover_tests(parallel_ankle_tests DISCOVERY_MODE PRE_TEST)
    if(TARGET five_bar_demo)
        add_test(NAME Examples.FiveBar COMMAND five_bar_demo)
    endif()
    if(TARGET constraints_demo)
        add_test(NAME Examples.Constraints COMMAND constraints_demo)
    endif()
    if(TARGET transmissions_demo)
        add_test(NAME Examples.Transmissions COMMAND transmissions_demo)
    endif()
    if(TARGET loaded_stiffness_demo)
        add_test(NAME Examples.LoadedStiffness COMMAND loaded_stiffness_demo)
    endif()
    if(TARGET linear_universal_ankle_demo)
        add_test(NAME Examples.LinearUniversalAnkle COMMAND linear_universal_ankle_demo)
    endif()
endif()

include(CMakePackageConfigHelpers)
install(TARGETS parallel_ankle EXPORT ParallelAnkleCoreTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(EXPORT ParallelAnkleCoreTargets NAMESPACE ParallelAnkle::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ParallelAnkle)
if(PARALLEL_ANKLE_WITH_YAML)
    install(TARGETS parallel_ankle_config EXPORT ParallelAnkleConfigTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(EXPORT ParallelAnkleConfigTargets NAMESPACE ParallelAnkle::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ParallelAnkle)
endif()
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY configs/ DESTINATION ${CMAKE_INSTALL_DATADIR}/ParallelAnkle/configs)
configure_package_config_file(cmake/ParallelAnkleConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/ParallelAnkleConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ParallelAnkle
    PATH_VARS CMAKE_INSTALL_DATADIR)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/ParallelAnkleConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ParallelAnkleConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/ParallelAnkleConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ParallelAnkle)
