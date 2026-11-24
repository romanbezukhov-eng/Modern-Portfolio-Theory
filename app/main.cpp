// Запуск: portfolio_app
// Создаёт data/random.txt со случайными данными, читает его и считает портфель.
// Если передать имя готового файла (portfolio_app файл.txt), берётся он, а не случайный.
#include <iostream>
#include <random>
#include <string>

#include "data_io.hpp"
#include "portfolio.hpp"

int main(int argc, char** argv) {
    const std::size_t n = 10;  // сколько активов в случайной таблице
    std::string path = std::string(PORTFOLIO_DATA_DIR) + "/random.txt";

    if (argc > 1) path = argv[1];
    else writeRandomCase(path, n, std::random_device{}());  // новая таблица каждый раз

    DataCase c = loadCase(path);
    PortfolioOptimizer opt(c.mu, c.sigma);
    PortfolioResult res = opt.solve(c.E);

    for (std::size_t i = 0; i < res.x.size(); ++i) std::cout << "x[" << i << "] = " << res.x[i] << "\n";
    std::cout << "E = " << res.E << ", R = " << res.R << "\n";
}
