#pragma once

#include <cstddef>
#include <utility>
#include <vector>

// Результат расчёта портфеля
struct PortfolioResult {
    std::vector<double> x;  // доли активов, длина n, сумма = 1 (у выброшенных активов 0)
    double R;               // дисперсия (риск) портфеля: x^T * sigma * x
    double E;               // фактическая ожидаемая доходность: mu^T * x
};


//Предполагается, что данные невырожденные, то есть sigma симметрична и положительно определена (это не проверяется, проверяется только симметрия)

class PortfolioOptimizer {
public:
    // mu    - ожидаемые доходности, длина n
    // sigma - ковариационная матрица n x n, хранится плоско по строкам: элемент (i, j) лежит в sigma[i * n + j]
    // eps   - точность для сравнения чисел
    // Бросает std::invalid_argument при неверных размерах или несимметричной sigma
    PortfolioOptimizer(std::vector<double>& mu, std::vector<double>& sigma, double eps = 1e-9); //изменены на ссылки


    // Ищет портфель минимального риска с ожидаемой доходностью E
    // Алгоритм (active - список активов, которые ещё участвуют в расчёте, сначала это все 0..n-1):
    //   1. Методом сопряжённых градиентов решаем две системы sigma * y1 = mu и sigma * y2 = 1
    //      (sigma берём только по активным активам)
    //   2. Считаем числа A = mu^T y1, B = mu^T y2, C = 1^T y2 и решаем систему 2 x 2
    //      A*l1 + B*l2 = E, B*l1 + C*l2 = 1, откуда получаем l1, l2
    //   3. Доли активов x = l1 * y1 + l2 * y2, риск R = l1 * E + l2
    //   4. Если есть отрицательная доля (в долг брать нельзя), выбрасываем этот актив
    //      и возвращаемся к шагу 1 уже для оставшихся. Иначе ответ готов
    // Бросает std::invalid_argument, если E вне диапазона доходностей активов, и std::runtime_error, если решение найти не удалось
    PortfolioResult solve(double E) const;


    std::size_t size() const { return n; }

private:
    std::size_t n;
    std::vector<double> mu;
    std::vector<double> sigma;
    double eps;

    double S(std::size_t i, std::size_t j) const { return sigma[i * n + j]; } //метод помощник, так как матрица хранится в строчку, так удобнее отриентироваться по эл-там

};
