#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> & operator *= (vector<double> & lhs, double rhs)
{
    for (auto & i : lhs)
        i *= rhs;
    return lhs;
}

vector<double> operator * (vector<double> lhs, double rhs){ return lhs *= rhs; }

vector<double> & operator /= (vector<double> & lhs, double rhs)
{
    for (auto & i : lhs)
        i /= rhs;
    return lhs;
}

vector<double> operator / (vector<double> lhs, double rhs){ return lhs /= rhs; }

double operator * (vector<double> lhs, vector<double> rhs)
{
    double res = 0;
    for (size_t i = 0; i < lhs.size(); ++i)
        res += lhs[i] * rhs[i];
    return res;
}

vector<double> operator * (vector<vector<double>> lhs, vector<double> rhs)
{
    vector<double> res(rhs.size());
    for (size_t i = 0; i < rhs.size(); ++i)
        res[i] = lhs[i] * rhs;
    return res;
}

vector<double> & operator += (vector<double> & lhs, const vector<double> rhs)
{
    for (size_t i = 0; i < lhs.size(); ++i)
        lhs[i] += rhs[i];
    return lhs;
}

vector<double> & operator -= (vector<double> & lhs, const vector<double> rhs)
{
    for (size_t i = 0; i < lhs.size(); ++i)
        lhs[i] -= rhs[i];
    return lhs;
}

vector<double> operator + (vector<double> lhs, vector<double> rhs)
{
    return lhs += rhs;
}

vector<double> operator - (vector<double> lhs, vector<double> rhs)
{
    return lhs -= rhs;
}

void init(vector<double> & c)
{
    for (uint i = 0; i < c.size(); ++i)
        c[i] = sin(M_PI * i / (c.size() - 1));
}

void borders(vector<double> & c)
{
    c[0] = 0;
    c[c.size()-1] = 0;
}

int main()
{
    int n = 20;
    double h = 1. / n;
    double dt = 0.0001;
    vector<double> c(n+1);
    init(c);
    borders(c);

    // матрицы
    vector<double> ad(n+1, 2 * h / 3); ad[0] = ad[n] = h / 3;
    vector<double> al(n, h / 6);
    vector<double> au(n, h / 6);

    vector<vector<double>> b(n+1);
    for (int i = 0; i <= n; ++i)
    {
        b[i].assign(n+1,0);
        for (int j = max(i-1, 0); j <= min(i+1, n); ++j)
        {
            if (i == j && i != 0 && i != n)
                b[i][j] = -2 / h;
            if (i == j+1 || i == j - 1)
                b[i][j] = 1 / h;
        }
    }
    // прогонка для формирования матрицы системы
    // прямая прогонка
    for (int i = 0; i < n; ++i)
    {
        ad[i+1] -= au[i] / ad[i] * al[i];
        b[i+1] -= b[i] * (al[i] / ad[i]);
    }
    // последняя строка
    b[n] /= ad[n];
    ad[n] = 1;
    // обратная прогонка
    for (int i = n; i > 0; --i)
    {
        for (int j = 0; j <= n; ++j)
        {
            b[i-1][j] -= b[i][j] * au[i-1];
            b[i-1][j] /= ad[i-1];
        }
    }

    // b -- матрица системы ОДУ
    FILE *fd = fopen("res.dat", "w");
    for (int k = 0; k < 1000; ++k)
    {
        // метод Рунге-Кутты 4 порядка точности
        auto k1 = b * c * dt,
             k2 = b * (c + k1 * .5) * dt,
             k3 = b * (c + k2 * .5) * dt,
             k4 = b * (c + k3) * dt;
        c += (k1 + k2 * 2. + k3 * 2. + k4) * (1./6);
        borders(c);
        if (k % 10 == 0)
        {
            for (int i = 0; i <= n; ++i)
                fprintf(fd, "%f %f %f\n", k * dt, (double)i / n, c[i]);
            fprintf(fd, "\n");
        }
    }
    fclose(fd);
    cout << "hello" << endl;
    return 0;
}
