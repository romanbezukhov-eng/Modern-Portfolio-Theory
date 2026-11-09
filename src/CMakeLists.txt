# Основная реализация.
add_library(portfolio portfolio.cpp)
target_include_directories(portfolio PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
