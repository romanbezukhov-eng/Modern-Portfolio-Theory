# Только заголовочный файл, поэтому INTERFACE-библиотека: у неё нет .cpp,
# она только говорит, где лежат заголовки.
add_library(data_io INTERFACE)
target_include_directories(data_io INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
