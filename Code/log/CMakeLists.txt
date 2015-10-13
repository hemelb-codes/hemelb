configure_file("Logger.in.cc" "${CMAKE_CURRENT_BINARY_DIR}/Logger.cc")
add_library(hemelb_log ${CMAKE_CURRENT_BINARY_DIR}/Logger.cc ${logger})
