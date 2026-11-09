# Один и тот же файл тестов гоняем на обеих реализациях.
add_executable(tests test_portfolio.cpp)
target_link_libraries(tests PRIVATE portfolio data_io)
target_compile_definitions(tests PRIVATE PORTFOLIO_DATA_DIR="${DATA_DIR}")
add_test(NAME tests COMMAND tests)

add_executable(tests_naive test_portfolio.cpp)
target_link_libraries(tests_naive PRIVATE portfolio_naive data_io)
target_compile_definitions(tests_naive PRIVATE PORTFOLIO_DATA_DIR="${DATA_DIR}")
add_test(NAME tests_naive COMMAND tests_naive)
