cmake_minimum_required(VERSION 3.20.1)

if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(speed)

set (CMAKE_CXX_STANDARD 17)

# find the correct path for isysroot with `xcrun --show-sdk-path`
IF("${CMAKE_SYSTEM}" MATCHES "Linux")
   SET(CMAKE_CXX_FLAGS "-std=c++17 -I/usr/local/include -I/usr/include")
ELSEIF("${CMAKE_HOST_SYSTEM_PROCESSOR}" MATCHES "arm64")
   message(STATUS CMAKE_HOST_SYSTEM_PROCESSOR = ${CMAKE_HOST_SYSTEM_PROCESSOR})
ELSE()
   SET(CMAKE_CXX_FLAGS "-std=c++17 -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk -I/usr/local/include")
ENDIF()

IF(UNIX AND NOT APPLE)
  set(LINUX TRUE)
  message(STATUS ">>> Linux")
ELSE()
  message(STATUS ">>> MacOS")
ENDIF()

# Boost
set(boost_min_ver 1.50.0)
set(boost_libs system filesystem timer date_time program_options)
find_package(Boost CONFIG)
#find_package(Boost ${boost_min_ver})

if(Boost_FOUND)
  find_package(Boost ${boost_min_ver} COMPONENTS ${boost_libs})
endif()

find_package(TBB REQUIRED)

if(TBB_FOUND)
  MESSAGE(STATUS ">>> Found TBB")
ENDIF()

find_package(LIBZIP REQUIRED)
if (LIBZIP_FOUND)
   MESSAGE(STATUS ">>> Found LIBZIP")
ELSE(LIBZIP_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build app without libzip.")
ENDIF()

add_executable(prune prune.cxx)
target_include_directories(prune PUBLIC)
target_link_libraries(prune ${Boost_LIBRARIES} TBB::tbb)
