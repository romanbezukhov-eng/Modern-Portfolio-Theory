#include "../portfolio.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace {

const double kPivotEps = 1e-14;  // главный элемент меньше этого считаем нулём

    // Метод Гаусса с выбором главного элемента по столбцу.
    // A - матрица m x m (плоско, по строкам), b - правая часть.
    // Ответ записывается в b. Возвращает false, если матрица вырождена.
    bool gaussSolve(std::vector<double>& A, std::vector<double>& b, std::size_t m) {
        for (std::size_t col = 0; col < m; ++col) {
            // ищем строку с максимальным по модулю элементом в столбце
            std::size_t pivot = col;
            double best = std::fabs(A[col * m + col]);
            for (std::size_t r = col + 1; r < m; ++r) {
                double v = std::fabs(A[r * m + col]);
                if (v > best) {
                    best = v;
                    pivot = r;
                }
            }
            if (best < kPivotEps) return false;

            if (pivot != col) {
                for (std::size_t c = 0; c < m; ++c) std::swap(A[col * m + c], A[pivot * m + c]);
                std::swap(b[col], b[pivot]);
            }

            // обнуляем столбец под диагональю
            for (std::size_t r = col + 1; r < m; ++r) {
                double f = A[r * m + col] / A[col * m + col];
                if (f == 0.0) continue;
                for (std::size_t c = col; c < m; ++c) A[r * m + c] -= f * A[col * m + c];
                b[r] -= f * b[col];
            }
        }

    // обратный ход
    for (std::size_t i = m; i-- > 0;) {
        double s = b[i];
        for (std::size_t c = i + 1; c < m; ++c) s -= A[i * m + c] * b[c];
        b[i] = s / A[i * m + i];
    }
    return true;
}

}  // конец namespace

PortfolioOptimizer::PortfolioOptimizer(std::vector<double> mu_, std::vector<double> sigma_,
                                       double eps_)
    : n(mu_.size()), mu(std::move(mu_)), sigma(std::move(sigma_)), eps(eps_) {
    if (n == 0) throw std::invalid_argument("mu пустой");
    if (sigma.size() != n * n)
        throw std::invalid_argument("размер sigma должен быть n * n");

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j)
            if (std::fabs(S(i, j) - S(j, i)) > 1e-12)
                throw std::invalid_argument("sigma должна быть симметричной");

    // Ищем дубликаты: активы с одинаковыми mu и одинаковыми строками sigma.
    std::vector<bool> removed(n, false);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < j; ++i) {
            if (removed[i]) continue;
            bool same = std::fabs(mu[i] - mu[j]) < eps;
            for (std::size_t k = 0; same && k < n; ++k)
                same = std::fabs(S(i, k) - S(j, k)) < eps;
            if (same) {
                removed[j] = true;
                dup.emplace_back(j, i);
                break;
            }
        }
        if (!removed[j]) active0.push_back(j);
    }
}

bool PortfolioOptimizer::solveSystem(const std::vector<std::size_t>& active, double E,
                                     std::vector<double>& A, std::vector<double>& b,
                                     std::vector<double>& xa) const {
    // Система размера m = k + 2 (k - число активных активов):
    //   sigma * x + h1 * mu + h2 * e = 0   (k строк)
    //   mu^T * x                     = E   (1 строка)
    //   e^T  * x                     = 1   (1 строка)
    // где h1 = lambda1 / 2, h2 = lambda2 / 2 - половинки множителей Лагранжа (как у Мартина).
    const std::size_t k = active.size();
    const std::size_t m = k + 2;

    // A и b созданы в solve один раз с запасом по памяти.
    // assign заново заполняет нулями нужное число элементов, новая память не выделяется.
    A.assign(m * m, 0.0);
    b.assign(m, 0.0);

    for (std::size_t i = 0; i < k; ++i) {
        for (std::size_t j = 0; j < k; ++j) A[i * m + j] = S(active[i], active[j]);
        A[i * m + k] = mu[active[i]];
        A[i * m + k + 1] = 1.0;

        A[k * m + i] = mu[active[i]];
        A[(k + 1) * m + i] = 1.0;
    }
    b[k] = E;
    b[k + 1] = 1.0;

    if (!gaussSolve(A, b, m)) return false;

    xa.assign(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(k));
    return true;
}

PortfolioResult PortfolioOptimizer::solve(double E) const {
    // E должна лежать между минимальной и максимальной доходностью активов
    double lo = mu[active0[0]];
    double hi = lo;
    for (std::size_t i : active0) {
        if (mu[i] < lo) lo = mu[i];
        if (mu[i] > hi) hi = mu[i];
    }
    if (E < lo - eps || E > hi + eps)
        throw std::invalid_argument("E вне диапазона доходностей активов");

    std::vector<std::size_t> active = active0;

    // Буферы создаём один раз, а не на каждом шаге цикла.
    // Максимальный размер системы: все активы + 2 уравнения.
    std::vector<double> A;
    std::vector<double> b;
    std::vector<double> xa;
    A.reserve((active.size() + 2) * (active.size() + 2));
    b.reserve(active.size() + 2);
    xa.reserve(active.size());

    // На каждом шаге выбрасывается ровно один актив, поэтому шагов не больше, чем активов (n-1).
    while (!active.empty()) {
        if (active.size() == 1) {
            // остался один актив: подходит только если его доходность равна E
            if (std::fabs(E - mu[active[0]]) > eps)
                throw std::runtime_error("решение не найдено");
            xa = {1.0};
        } else if (!solveSystem(active, E, A, b, xa)) {
            // Без этого throw xa остался бы пустым и функция вернула бы нулевой портфель.
            throw std::runtime_error("система вырождена, решение не найдено");
        }

        // ищем самую отрицательную долю
        std::size_t worst = xa.size();
        double minv = -eps;
        for (std::size_t i = 0; i < xa.size(); ++i) {
            if (xa[i] < minv) {
                minv = xa[i];
                worst = i;
            }
        }

        if (worst == xa.size()) {
            // все доли неотрицательны - собираем ответ длины n
            PortfolioResult res;
            res.x.assign(n, 0.0);
            for (std::size_t i = 0; i < active.size(); ++i)
                res.x[active[i]] = xa[i] > 0.0 ? xa[i] : 0.0;

            res.V = 0.0;
            res.E = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                res.E += mu[i] * res.x[i];
                for (std::size_t j = 0; j < n; ++j) res.V += res.x[i] * S(i, j) * res.x[j];
            }
            return res;
        }

        // выбрасываем актив с отрицательной долей и решаем заново
        active.erase(active.begin() + static_cast<std::ptrdiff_t>(worst));
    }

    // Сюда не попадаем: при одном активе мы либо вернули ответ, либо бросили исключение.
    throw std::runtime_error("решение не найдено");
}



// геометрическая интерпритация
// гаусс хотспот?
// профилировщик
// while true
// можно ли создавать вектор не в цикле
// тесты, рассмотреть пример из хадачник
// crjkmrj d while
// итерационные методы вместо гаусса
