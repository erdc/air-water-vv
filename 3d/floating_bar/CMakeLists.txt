cmake_minimum_required(VERSION 2.8)

SET(MODELS
  BarApp
  )

#--------------------------------------------------------------
# Find the Chrono package with required components
#--------------------------------------------------------------

# Invoke find_package in CONFIG mode

find_package(Chrono
             COMPONENTS Irrlicht FEA
             CONFIG
	     )
message("data dir")
message(${CHRONO_DATA_DIR}) 
# If Chrono and/or the required component(s) were not found, return now.

if(NOT Chrono_FOUND)
  message("Could not find requirements for Chrono tutorials")
  return()
endif()

include_directories(${CHRONO_INCLUDE_DIRS})

#--------------------------------------------------------------
# Compilation flags
#--------------------------------------------------------------

set(COMPILE_FLAGS ${CHRONO_CXX_FLAGS})

# Disable some warnings triggered by Irrlicht (Windows only)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    SET(COMPILE_FLAGS "${COMPILE_FLAGS} /wd4275")
endif()


foreach(PROGRAM ${MODELS})

  message(STATUS "...add ${PROGRAM}")

  add_executable(${PROGRAM}  "${PROGRAM}.cpp")
  source_group(""  FILES "${PROGRAM}.cpp")

  set_target_properties(${PROGRAM} PROPERTIES
    FOLDER demos
    COMPILE_FLAGS "${COMPILE_FLAGS}"
    COMPILE_DEFINITIONS "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\""
    LINK_FLAGS "${COMPILE_FLAGS} ${CHRONO_LINKER_FLAGS}"
  )
  target_link_libraries(${PROGRAM} ${CHRONO_LIBRARIES})

endforeach(PROGRAM)
