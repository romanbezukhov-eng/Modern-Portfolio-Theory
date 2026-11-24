#pragma once

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Один пример из папки data/.
// Формат файла: n, потом mu (n чисел), sigma (n*n чисел), E, потом (необязательно) ожидаемые x.
// Строки, которые начинаются с '#', - комментарии.
struct DataCase {
    std::vector<double> mu;
    std::vector<double> sigma;
    double E = 0.0;
    std::vector<double> expected;  // может быть пустым
};

inline DataCase loadCase(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("не могу открыть файл: " + path);

    std::stringstream text;  // сюда собираем всё, кроме комментариев
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') continue;
        text << line << '\n';
    }

    std::size_t n = 0;
    if (!(text >> n) || n == 0) throw std::runtime_error("неверное n в файле: " + path);

    DataCase c;
    c.mu.resize(n);
    c.sigma.resize(n * n);
    for (double& v : c.mu) text >> v;
    for (double& v : c.sigma) text >> v;
    text >> c.E;
    if (!text) throw std::runtime_error("файл закончился раньше времени: " + path);

    double v;
    while (text >> v) c.expected.push_back(v);
    return c;
}

// Создаёт файл со случайными, но правдоподобными данными (тот же формат, что и в data/).
// Каждый вызов с новым seed даёт новую таблицу, файл просто перезаписывается.
//
// Идея: все акции зависят от одного общего рынка (поэтому они коррелированы),
// и у каждой есть ещё своя собственная случайность.
//   beta[i]  - насколько актив чувствителен к рынку (0.5 .. 1.5)
//   mu[i]    = 3% + 6% * beta[i] + небольшой шум   (чем рискованнее, тем доходнее, около 3% - 13%)
//   sigma[i][j] = beta[i] * beta[j] * 0.15^2        (общий рынок, i != j)
//   sigma[i][i] = beta[i]^2 * 0.15^2 + own[i]^2      (плюс своя случайность, 10% - 30%)
// Такая sigma всегда симметричная и положительно определённая, разброс по активу 17% - 33%.
// E берётся как среднее из mu (значит лежит между минимумом и максимумом).
inline void writeRandomCase(const std::string& path, std::size_t n, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> betaDist(0.5, 1.5);
    std::uniform_real_distribution<double> noiseDist(-0.02, 0.02);
    std::uniform_real_distribution<double> ownDist(0.10, 0.30);

    const double marketVol = 0.15;
    std::vector<double> beta(n), own(n), mu(n);
    double E = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        beta[i] = betaDist(gen);
        own[i] = ownDist(gen);
        mu[i] = 0.03 + 0.06 * beta[i] + noiseDist(gen);
        E += mu[i] / static_cast<double>(n);
    }

    std::ofstream out(path);
    if (!out) throw std::runtime_error("не могу записать файл: " + path);
    out << std::setprecision(10);
    out << "# случайные данные, seed = " << seed << "\n" << n << "\n";
    for (std::size_t i = 0; i < n; ++i) out << mu[i] << (i + 1 < n ? ' ' : '\n');
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double s = beta[i] * beta[j] * marketVol * marketVol;
            if (i == j) s += own[i] * own[i];
            out << s << (j + 1 < n ? ' ' : '\n');
        }
    }
    out << E << "\n";
}
