cmake_minimum_required(VERSION 3.10)
project("vecimg")

include_directories("include")

add_executable(vecimg
	src/main.c
	src/stbi.c
)

target_compile_options(vecimg PUBLIC -lm -g -Wall -O2)
target_link_options(vecimg PUBLIC -lm)
