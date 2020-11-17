target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/ConcernTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MockConcern.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MockIteratedAction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/StepManagerTests.cc
)
