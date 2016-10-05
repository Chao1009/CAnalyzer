#ifndef C_MATRIX_H
#define C_MATRIX_H

#include <cstddef>
#include <vector>
#include <iostream>
#include <initializer_list>

class CMatrix
{
public:
    CMatrix(const size_t &n);
    CMatrix(const size_t &n, const size_t &m);
    CMatrix(const size_t &n, const size_t &m, double **ele);
    CMatrix(const CMatrix &rhs);
    CMatrix(CMatrix &&rhs);
    virtual ~CMatrix();

    void FillElements(const std::initializer_list<double> &ele);
    double At(const size_t &n, const size_t &m) const;
    size_t DimN() const {return dim_n;};
    size_t DimM() const {return dim_m;};
    std::vector<double> Row(const size_t &i) const;
    std::vector<double> Column(const size_t &j) const;
    CMatrix Add(const CMatrix &m, const double &scale = 1.) const;
    CMatrix Multiply(const CMatrix &m) const;
    CMatrix Scale(const double &factor) const;
    CMatrix UpperLeft(const size_t &k) const;
    double Det(const size_t &k) const;

    // available only for square matrix
    bool SquareCheck(const std::string &func_name = "") const;
    bool IsSymmetric() const;
    bool IsPositiveDefinite() const;
    CMatrix Identity() const;
    CMatrix Power(const int &n) const;
    double Det() const;
    double Trace() const;
    CMatrix Cofactor() const;
    CMatrix Transpose() const;
    void TransposeSelf();
    CMatrix Inverse() const;
    CMatrix Cholesky() const;
    CMatrix LowerTri() const;
    CMatrix UpperTri() const;

    double &operator ()(const size_t &n, const size_t &m)
    {
        return elements[n][m];
    }

    CMatrix operator +(const CMatrix &rhs) const
    {
        return Add(rhs);
    }

    CMatrix operator -(const CMatrix &rhs) const
    {
        return Add(rhs, -1.);
    }

    CMatrix operator *(const CMatrix &rhs) const
    {
        return Multiply(rhs);
    }

    CMatrix operator *(const double &factor) const
    {
        return Scale(factor);
    }

    CMatrix operator /(const double &factor) const
    {
        return Scale(1./factor);
    }

    CMatrix operator ^(const int &n) const
    {
        return Power(n);
    }

    void operator =(const std::initializer_list<double> &ele)
    {
        FillElements(ele);
    }

    double **GetPtr() {return elements;};
    void NullPtr(){elements = nullptr;};

private:
    void initialize(const bool &zero = false);

private:
    size_t dim_n, dim_m;
    double **elements;
};

std::ostream &operator <<(std::ostream &os, const CMatrix &m);
CMatrix operator *(const double &factor, const CMatrix &m);

namespace cana
{
    CMatrix det(const CMatrix &m);
    CMatrix inv(const CMatrix &m);
    CMatrix transpose(const CMatrix &m);
    CMatrix trace(const CMatrix &m);
    CMatrix cholesky(const CMatrix &m);
    CMatrix tril(const CMatrix &m);
    CMatrix triu(const CMatrix &m);
    CMatrix power(const CMatrix &m, const int &n);
    CMatrix identity(const size_t &n);
};

#endif

