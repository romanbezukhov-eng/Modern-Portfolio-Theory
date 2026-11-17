#include "portfolio.hpp"

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace {

    const double kDetEps = 1e-14;  // определитель системы 2 x 2 меньше этого считаем нулём

    // скалярное произведение двух векторов одинаковой длины для Шура
    double dot(const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
        return s;
    }
    // Метод сопряжённых градиентов для системы sigma_active * y = b.
    // sigma_active - подматрица sigma из строк и столбцов с номерами из active (саму матрицу не копируем).
    // sigma - плоская матрица n x n по строкам, b и y имеют длину active.size().
    // r, p, ap - рабочие буферы, их создаёт solve один раз.
    // Метод работает, потому что sigma_active симметрична и положительно определена.
    // Идея: стартуем с y = 0 и на каждом шаге сдвигаемся вдоль направления p, по которому невязка r = b - A*y убывает быстрее всего,
    // при этом новое направление строим так, чтобы оно не портило то, что уже улучшили (сопряжённость).
    // Останавливаемся, когда длина невязки стала меньше eps от длины b.
    // Возвращает false, если за 10 * k шагов не сошлись или матрица оказалась не положительно определённой.
    bool conjugateGradient(const std::vector<double>& sigma, std::size_t n,
                           const std::vector<std::size_t>& active, const std::vector<double>& b,
                           std::vector<double>& y, std::vector<double>& r, std::vector<double>& p,
                           std::vector<double>& ap, double eps) {
        const std::size_t k = active.size();

        y.assign(k, 0.0);
        r = b;  // при y = 0 невязка b - A*y равна просто b
        p = r;  // первое направление - сама невязка
        ap.assign(k, 0.0);

        double rs = dot(r, r);
        const double tol2 = eps * eps * dot(b, b);  // сравниваем квадраты, чтобы не брать корень
        const std::size_t maxIter = 10 * k;         // в точной арифметике хватило бы k шагов

        for (std::size_t it = 0; it < maxIter && rs > tol2; ++it) {
            // ap = sigma_active * p
            for (std::size_t i = 0; i < k; ++i) {
                const double* row = &sigma[active[i] * n];
                double s = 0.0;
                for (std::size_t j = 0; j < k; ++j) s += row[active[j]] * p[j];
                ap[i] = s;
            }

            const double pap = dot(p, ap);
            if (pap <= 0.0) return false;  // для положительно определённой матрицы так быть не должно

            const double alpha = rs / pap;  // длина шага вдоль p
            for (std::size_t i = 0; i < k; ++i) {
                y[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }

            const double rsNew = dot(r, r);
            const double beta = rsNew / rs;
            for (std::size_t i = 0; i < k; ++i) p[i] = r[i] + beta * p[i];
            rs = rsNew;
        }
        return rs <= tol2;

        // Метод Шура: по готовым y1 = sigma^-1 * mu и y2 = sigma^-1 * 1 (по активным активам) собираем портфель.
        // Считаем три числа A = mu^T y1, B = mu^T y2, C = 1^T y2 и решаем систему 2 x 2
        //   A*l1 + B*l2 = E
        //   B*l1 + C*l2 = 1
        // по формулам Крамера. Потом w = l1 * y1 + l2 * y2 и риск R = l1 * E + l2.
        // w записывается в w (длина active.size()), риск в R.
        // Возвращает false, если определитель системы 2 x 2 равен нулю (например, у всех активных активов одинаковое mu).
    }

    bool schur(const std::vector<double>& mu, const std::vector<std::size_t>& active,
           const std::vector<double>& y1, const std::vector<double>& y2, double E,
           std::vector<double>& w, double& R) {
        const std::size_t k = active.size();

        double A = 0.0;
        double B = 0.0;
        double C = 0.0;
        for (std::size_t i = 0; i < k; ++i) {
            const double m = mu[active[i]];
            A += m * y1[i];
            B += m * y2[i];
            C += y2[i];
        }

        const double D = A * C - B * B;
        if (std::fabs(D) < kDetEps) return false;

        const double l1 = (E * C - B) / D;
        const double l2 = (A - E * B) / D;

        w.resize(k);
        for (std::size_t i = 0; i < k; ++i) w[i] = l1 * y1[i] + l2 * y2[i];
        R = l1 * E + l2;  // w^T * sigma * w, но без двойного цикла
        return true;
    }


}  // конец namespace


// проверяет, что mu не пустой, размер sigma равен n * n и sigma симметрична
// ничего не считает, данные считаем невырожденными (sigma положительно определена)
PortfolioOptimizer::PortfolioOptimizer(std::vector<double>& mu_, std::vector<double>& sigma_,
                                       double eps_)

// n (размер вектора mu_) - общее количество активов
// mu - вектор ожидаемых доходностей для каждого актива (длина n)
// sigma - плоская матрица ковариаций (размер n * n), описывает риски и взаимосвязи активов
// eps - погрешность для сравнения вещественных чисел (например, 1e-7)

    : n(mu_.size()), mu(std::move(mu_)), sigma(std::move(sigma_)), eps(eps_) {
    if (n == 0) throw std::invalid_argument("mu пустой");
    if (sigma.size() != n * n)
        throw std::invalid_argument("размер sigma должен быть n * n");

    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = i + 1; j < n; ++j)
            if (std::fabs(S(i, j) - S(j, i)) > 1e-12)
                throw std::invalid_argument("sigma должна быть симметричной");
}


// передаем желаемую доходность E, ищем портфель с мин. риском
PortfolioResult PortfolioOptimizer::solve(double E) const {
    // E должна лежать между минимальной и максимальной доходностью активов
    double lo = mu[0];
    double hi = lo;
    for (std::size_t i = 1; i < n; ++i) {
        if (mu[i] < lo) lo = mu[i];
        if (mu[i] > hi) hi = mu[i];
    }
    if (E < lo - eps || E > hi + eps)
        throw std::invalid_argument("E вне диапазона доходностей активов");

    // active - активы, которые ещё участвуют в расчёте, сначала все
    std::vector<std::size_t> active(n);
    std::iota(active.begin(), active.end(), std::size_t{0});

    // Буферы создаём один раз, а не на каждом шаге цикла.
    // Максимальный размер - все активы, дальше они только уменьшаются.
    std::vector<double> rhs1, rhs2;  // правые части: mu и вектор из единиц (по активным активам)
    std::vector<double> y1, y2;      // решения sigma * y1 = mu и sigma * y2 = 1
    std::vector<double> r, p, ap;    // рабочие буферы метода сопряжённых градиентов
    std::vector<double> w;           // доли активных активов
    rhs1.reserve(n);
    rhs2.reserve(n);
    y1.reserve(n);
    y2.reserve(n);
    r.reserve(n);
    p.reserve(n);
    ap.reserve(n);
    w.reserve(n);

    // На каждом шаге выбрасывается ровно один актив, поэтому шагов не больше, чем активов (n-1).
    // Когда остался один актив, мы либо возвращаем ответ, либо бросаем исключение, так что цикл конечный.
    while (!active.empty()) {
        const std::size_t k = active.size();
        double R = 0.0;

        if (k == 1) {
            // остался один актив: подходит только если его доходность равна E
            // (метод Шура тут не работает, определитель D равен нулю)
            if (std::fabs(E - mu[active[0]]) > eps)
                throw std::runtime_error("решение не найдено");
            w.assign(1, 1.0);
            R = S(active[0], active[0]);
        } else {
            rhs1.resize(k);
            for (std::size_t i = 0; i < k; ++i) rhs1[i] = mu[active[i]];
            rhs2.assign(k, 1.0);

            if (!conjugateGradient(sigma, n, active, rhs1, y1, r, p, ap, eps) ||
                !conjugateGradient(sigma, n, active, rhs2, y2, r, p, ap, eps))
                throw std::runtime_error("итерационный метод не сошёлся, решение не найдено");

            if (!schur(mu, active, y1, y2, E, w, R))
                throw std::runtime_error("система 2 x 2 вырождена, решение не найдено");
        }

        // ищем самую отрицательную долю
        std::size_t worst = k;
        double minv = -eps;
        for (std::size_t i = 0; i < k; ++i) {
            if (w[i] < minv) {
                minv = w[i];
                worst = i;
            }
        }

        if (worst == k) {
            // все доли неотрицательны - собираем ответ длины n
            PortfolioResult res;
            res.x.assign(n, 0.0);
            for (std::size_t i = 0; i < k; ++i)
                res.x[active[i]] = w[i] > 0.0 ? w[i] : 0.0;
            res.E = dot(mu, res.x);
            res.R = R;
            return res;
        }

        // выбрасываем актив с отрицательной долей и считаем заново
        active.erase(active.begin() + static_cast<std::ptrdiff_t>(worst));
    }
}



// геометрическая интерпритация
// гаусс хотспот?
// профилировщик
// while true
// можно ли создавать вектор не в цикле
// тесты, рассмотреть пример из хадачник
// crjkmrj d while
// итерационные методы вместо гаусса
