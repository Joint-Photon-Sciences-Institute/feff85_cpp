add_library(feffpath
  paths.cpp pathsd.cpp prcrit.cpp ccrit.cpp mcrith.cpp mcritk.cpp
  mrb.cpp mpprmd.cpp mpprmp.cpp timrep.cpp phash.cpp
  heap.cpp ipack.cpp sortix.cpp repath.cpp outcrt.cpp
  ffmod4.cpp
)
target_include_directories(feffpath PUBLIC ${CMAKE_SOURCE_DIR}/include PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(feffpath PUBLIC feffcommon feffmath feffpar nlohmann_json::nlohmann_json)
