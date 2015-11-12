#------Capture build environment -------------
find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
  RESULT_VARIABLE rev_ok OUTPUT_VARIABLE HEMELB_REVISION_NUMBER OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT rev_ok EQUAL 0)
  message("Could not get revision number from mercurial, looking for revision_info.txt")
  execute_process(COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../revision_info.txt RESULT_VARIABLE rev_ok OUTPUT_VARIABLE HEMELB_REVISION_NUMBER OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
message("Revision number: ${HEMELB_REVISION_NUMBER}")

execute_process(COMMAND date OUTPUT_VARIABLE HEMELB_BUILD_TIME OUTPUT_STRIP_TRAILING_WHITESPACE)

