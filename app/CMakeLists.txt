add_executable(portfolio_app main.cpp)
target_link_libraries(portfolio_app PRIVATE portfolio data_io)
target_compile_definitions(portfolio_app PRIVATE PORTFOLIO_DATA_DIR="${DATA_DIR}")

# make run_martin / cmake --build build --target run_martin
# запускает пример из статьи Мартина (6% и 9%).
add_custom_target(run_martin
    COMMAND portfolio_app ${DATA_DIR}/martin1955_e6.txt
    COMMAND portfolio_app ${DATA_DIR}/martin1955_e9.txt
    DEPENDS portfolio_app
    USES_TERMINAL)
