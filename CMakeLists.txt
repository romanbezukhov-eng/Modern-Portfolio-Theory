cmake_minimum_required(VERSION 3.16)
project(itlab_portfolio CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# По умолчанию собираем с оптимизацией (для профилировки/бенчмарков это важно).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra)
endif()

# Путь к папке с данными, чтобы main и тесты находили файлы из любой папки запуска.
# Переменная видна во всех подпапках.
set(DATA_DIR "${PROJECT_SOURCE_DIR}/data")

enable_testing()

add_subdirectory(common)   # чтение данных (только заголовок)
add_subdirectory(src)      # основная реализация
add_subdirectory(naive)    # наивная реализация, для сравнения
add_subdirectory(app)      # программа для запуска с числами
add_subdirectory(tests)    # тесты
# В data/ нет кода, только текстовые файлы, поэтому CMake там не нужен.
