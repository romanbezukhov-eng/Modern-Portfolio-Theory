# Наивная реализация, хранится для сравнения.
# Класс там называется так же (PortfolioOptimizer), поэтому в один исполняемый файл
# вместе с основной её линковать нельзя. Каждая программа берёт одну из двух.
add_library(portfolio_naive portfolio.cpp)
target_include_directories(portfolio_naive PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
