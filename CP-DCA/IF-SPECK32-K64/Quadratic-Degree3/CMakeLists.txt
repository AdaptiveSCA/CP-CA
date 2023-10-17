CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(CPDCA)

include_directories(include)

aux_source_directory(src/m4ri CPDCA_SRC)

add_library(CPDCA_LIB ${CPDCA_SRC})

add_executable(CPDCA test/white_box_arx.c)

target_link_libraries(CPDCA CPDCA_LIB)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O2")