set(MESHHEXER_CLI_SOURCES
  meshhexer-cli.cpp
)

set(MESHHEXER_CLI_HEADERS
  timing.hpp
)

add_executable(meshhexer-cli ${MESHHEXER_CLI_SOURCES})
target_sources(meshhexer-cli
  PUBLIC FILE_SET headers
  TYPE HEADERS
  FILES ${MESHHEXER_CLI_HEADERS})
target_link_libraries(meshhexer-cli PRIVATE MeshHexer)

if(MESHHEXER_INSTALL)
  install(
    TARGETS meshhexer-cli
  )
endif(MESHHEXER_INSTALL)
