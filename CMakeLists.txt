cmake_minimum_required(VERSION 3.17)
project(RTC)

set(CMAKE_CXX_STANDARD 14)

add_executable(RTC main.cpp)
target_link_libraries(RTC -lopengl32 -lglut32 -lGlaux -lGLU32 -lglui32)