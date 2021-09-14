cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
   project(cca_inttime)

   set(CMAKE_BUILD_TYPE RELEASE)
   
   
    set(CMAKE_C_FLAGS "-Wall")
    set(CMAKE_C_FLAGS_DEBUG "-O0 -ggdb")
    set(CMAKE_C_FLAGS_RELEASE "-pipe -O3 -march=native -fomit-frame-pointer -fopt-info-vec-missed ")
    set(CMAKE_CXX_FLAGS "-Wall")
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb")
    set(CMAKE_CXX_FLAGS_RELEASE "-pipe -O3 -march=native -fomit-frame-pointer -fopt-info-vec-missed ")
   
   include_directories(${HDF5_INCLUDE_DIRS})

   find_package(HDF5 REQUIRED COMPONENTS C CXX)
   find_package(OpenMP REQUIRED)
   find_package(Protobuf REQUIRED)


   add_executable(cca_inttime cca_inttime ccset.pb.cc)
   set_property(TARGET cca_inttime PROPERTY CXX_STANDARD 20)
   target_link_libraries(cca_inttime ${HDF5_HL_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES} ${PROTOBUF_LIBRARY} OpenMP::OpenMP_CXX)

   include(CheckIPOSupported)
   check_ipo_supported(RESULT hasipo OUTPUT error)
   if( hasipo )
      message(STATUS "IPO / LTO enabled")
      set_property(TARGET cca_inttime PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
   else()
      message(STATUS "IPO / LTO not supported: <${error}>")
   endif()
      
   message(STATUS "INCLUDE LOCATION" ${HDF5_INCLUDE_DIRS})
   message(STATUS "version" ${HDF5_VERSION})
   message(STATUS "DEFINITIONS" ${HDF5_DEFINITIONS})
   message(STATUS "LIBRARIES" ${HDF5_LIBRARIES})
   message(STATUS "HL_LIBRARIES" ${HDF5_HL_LIBRARIES})
