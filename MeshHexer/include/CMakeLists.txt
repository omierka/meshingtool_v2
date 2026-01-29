set(MESHHEXER_PUBLIC_HEADERS
  meshhexer/meshhexer.hpp
  meshhexer/types.hpp
)

target_sources(MeshHexer
  PUBLIC FILE_SET public_headers
  TYPE HEADERS
  FILES ${MESHHEXER_PUBLIC_HEADERS})