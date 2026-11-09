#pragma once

#include <cstddef>
#include <utility>
#include <vector>

// Результат расчёта портфеля.
struct PortfolioResult {
    std::vector<double> x;  // доли активов, длина n, сумма = 1 (у выброшенных активов 0)
    double V;               // дисперсия (риск) портфеля: x^T * sigma * x
    double E;               // фактическая ожидаемая доходность: mu^T * x
};

class PortfolioOptimizer {
public:
    // mu    - ожидаемые доходности, длина n
    // sigma - ковариационная матрица n x n, хранится плоско по строкам: элемент (i, j) лежит в sigma[i * n + j]
    // eps   - точность для сравнения чисел
    // Бросает std::invalid_argument при неверных размерах или несимметричной sigma.
    PortfolioOptimizer(std::vector<double> mu, std::vector<double> sigma, double eps = 1e-9); //поменять на ссылки


    // Бросает std::invalid_argument, если E вне диапазона доходностей активов, и std::runtime_error, если решение найти не удалось.
    PortfolioResult solve(double E) const;

    std::size_t size() const { return n; }

    // Активы, выброшенные как дубликаты хранятся как пары (выброшенный, оставленный).
    const std::vector<std::pair<std::size_t, std::size_t>>& duplicates() const { return dup; }

private:
    std::size_t n;
    std::vector<double> mu;
    std::vector<double> sigma;
    double eps;
    std::vector<std::size_t> active0;                        // активы без дубликатов
    std::vector<std::pair<std::size_t, std::size_t>> dup;    // (выброшенный, оставленный)

    double S(std::size_t i, std::size_t j) const { return sigma[i * n + j]; }

    // Решает систему Лагранжа для активов из active, ответ кладёт в xa.
    // Возвращает false, если система вырождена.
    bool solveSystem(const std::vector<std::size_t>& active, double E,
                     std::vector<double>& xa) const;
};
