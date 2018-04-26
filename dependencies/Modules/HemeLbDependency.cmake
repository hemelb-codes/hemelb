
macro(hemelb_dependency NAME STEP)
  include(${HEMELB_DEPENDENCIES_PATH}/${NAME}/${STEP}.cmake)
endmacro()
