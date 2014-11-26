add_library(hemelb_unittests main.cc)

add_executable(tetrahedron redblood/tetrahedron.cc)
target_link_libraries(tetrahedron
  ${heme_libraries}
	${MPI_LIBRARIES}
	${PARMETIS_LIBRARIES}
	${TINYXML_LIBRARIES}
	${Boost_LIBRARIES}
	${CTEMPLATE_LIBRARIES}
	${ZLIB_LIBRARIES}
  ${MPWide_LIBRARIES}
)
