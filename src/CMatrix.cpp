#include "CMatrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>


std::ostream &operator <<(std::ostream &os, const CMatrix &m)
{
    for(size_t i = 0; i < m.DimN(); ++i)
    {
        for(size_t j = 0; j < m.DimM(); ++j)
        {
            os << std::setw(12) << m.At(i, j);
        }
        os << std::endl;
    }
    return os;
}

CMatrix operator *(const double &factor, const CMatrix &m)
{
    return m.Scale(factor);
}

// constructors
// n x n matrix
CMatrix::CMatrix(const size_t &n)
: dim_n(n), dim_m(n), elements(nullptr)
{
    initialize(true);
}

// n x m matrix
CMatrix::CMatrix(const size_t &n, const size_t &m)
: dim_n(n), dim_m(m), elements(nullptr)
{
    initialize(true);
}

// n x m matrix with initial values
CMatrix::CMatrix(const size_t &n, const size_t &m, double **ele)
: dim_n(n), dim_m(m), elements(nullptr)
{
    elements = ele;
}

CMatrix::CMatrix(const CMatrix &rhs)
: dim_n(rhs.DimN()), dim_m(rhs.DimM())
{
    initialize(false);
    for(size_t i = 0; i < dim_n; ++i)
    {
        for(size_t j = 0; j < dim_m ; ++j)
        {
            elements[i][j] = rhs.At(i, j);
        }
    }
}

CMatrix::CMatrix(CMatrix && rhs)
: dim_n(rhs.DimN()), dim_m(rhs.DimM())
{
    elements = rhs.GetPtr();
    rhs.NullPtr();
}

CMatrix::~CMatrix()
{
    if(!elements) return;
    for(size_t i = 0; i < dim_n; ++i)
    {
        delete [] elements[i];
    }
    delete [] elements;
}

void CMatrix::initialize(const bool &zero)
{
    elements = new double*[dim_n];
    for(size_t i = 0; i < dim_n; ++i)
    {
        elements[i] = new double[dim_m];
    }

    if(zero) {
        for(size_t i = 0; i < dim_n; ++i)
        {
            for(size_t j = 0; j < dim_m; ++j)
            {
                elements[i][j] = 0;
            }
        }
    }
}

void CMatrix::FillElements(const std::initializer_list<double> &ele)
{

    size_t i = 0, j = 0;
    for(double x : ele)
    {
        if(j >= dim_m) {
            j -= dim_m;
            i ++;
        }
        if(i >= dim_n)
            break;
        elements[i][j++] = x;
    }
}

double CMatrix::At(const size_t &n, const size_t &m) const
{
    // its user's responsibility to not exceed the matrix range
    return elements[n][m];
}

std::vector<double> CMatrix::Row(const size_t &i) const
{
    std::vector<double> res;
    if(i < dim_n)
    {
        for(size_t j = 0; j < dim_m; ++j)
        {
            res.push_back(elements[i][j]);
        }
    }
    return res;
}

std::vector<double> CMatrix::Column(const size_t &j) const
{
    std::vector<double> res;
    if(j < dim_m)
    {
        for(size_t i = 0; i < dim_n; ++i)
        {
            res.push_back(elements[i][j]);
        }
    }
    return res;
}

CMatrix CMatrix::Add(const CMatrix &rhs, const double &scale) const
{
    if(rhs.DimN() != dim_n ||
       rhs.DimM() != dim_m )
    {
        std::cerr << "CMatrix Error: a "
                  << dim_n << "x" << dim_m
                  << " matrix cannot add by a "
                  << rhs.DimN() << "x" << rhs.DimM()
                  << " matrix." << std::endl;
        return CMatrix(0);
    }

    CMatrix res(dim_n, dim_m);

    for(size_t i = 0; i < dim_n; ++i)
    {
        for(size_t j = 0; j < dim_m; ++j)
        {
            res(i, j) = At(i, j) + rhs.At(i, j)*scale;
        }
    }
    return res;
}

CMatrix CMatrix::Multiply(const CMatrix &rhs) const
{
    if(rhs.DimN() != dim_m)
    {
        std::cerr << "CMatrix Error: a "
                  << dim_n << "x" << dim_m
                  << " matrix cannot multiply by a "
                  << rhs.DimN() << "x" << rhs.DimM()
                  << " matrix." << std::endl;
        return CMatrix(0);
    }

    CMatrix res(dim_n, rhs.DimM());

    for(size_t i = 0; i < dim_n; ++i)
    {
        for(size_t j = 0; j < rhs.DimM(); ++j)
        {
            for(size_t m = 0; m < dim_m; ++m)
            {
                res(i, j) += At(i, m)*rhs.At(m, j);
            }
        }
    }

    return res;
}


CMatrix CMatrix::Identity() const
{
    if(!SquareCheck())
        return CMatrix(0);

    CMatrix res(dim_n);
    for(size_t i = 0; i < dim_n; ++i)
        res(i, i) = 1;

    return res;
}

CMatrix CMatrix::Power(const int &n) const
{
    if(n < 0)
        return Inverse().Power(-n);

    if(n == 0)
        return Identity();

    if(n == 1)
        return CMatrix(*this);

    CMatrix temp = Power(n/2);

    if(n%2 == 0)
        return temp*temp;
    else
        return temp*temp*(*this);
}

CMatrix CMatrix::Scale(const double &factor) const
{
    CMatrix res(dim_n, dim_m);
    for(size_t i = 0; i < dim_n; ++i)
    {
        for(size_t j = 0; j < dim_m; ++j)
        {
            res(i, j) = At(i, j)*factor;
        }
    }
    return res;
}

CMatrix CMatrix::UpperLeft(const size_t &k) const
{
    if(k > dim_n || k > dim_m)
    {
        std::cerr << "CMatrix Error: trying to get the "
                  << k << " upper left determinant, "
                  << "while the matrix is " << dim_n << " x " << dim_n
                  << std::endl;
        return 0.;
    }

    CMatrix res(k);

    for(size_t i = 0; i < k; ++i)
        for(size_t j = 0; j < k; ++j)
            res(i, j) = At(i, j);

    return res;
}

bool CMatrix::SquareCheck() const
{
    if(dim_n != dim_m)
    {
        std::cerr << "CMatrix Error: This is not a square matrix."
                  << std::endl;
        return false;
    }

    if(dim_n < 1)
    {
        std::cerr << "CMatrix Error: This is a 0x0 matrix, it is not correctly initialized."
                  << std::endl;
        return false;
    }

    return true;
}

double CMatrix::Det() const
{
    if(!SquareCheck())
        return 0.;

    if(dim_n == 1) { // special case
        return At(0, 0);
    } else if(dim_n == 2) { // special case
        return At(0, 0)*At(1, 1) - At(1, 0)*At(0, 1);
    } else { // general case, det of minors
        double det = 0.;
        for(size_t i = 0; i < dim_n; ++i)
        {
            CMatrix minors(dim_n - 1);
            for(size_t j = 1; j < dim_n; ++j)
            {
                size_t k2 = 0;
                for(size_t k = 0; k < dim_n; ++k)
                {
                    if(k == i) continue;
                    minors(j-1, k2) = At(j, k);
                    ++k2;
                }
            }
            det += pow(-1.0, i) * At(0, i) * minors.Det();
        }
        return det;
    }
}

double CMatrix::Det(const size_t &k) const
{
    return UpperLeft(k).Det();
}

double CMatrix::Trace() const
{
    if(!SquareCheck())
        return 0.;

    double res = 0.;
    for(size_t i = 0; i < dim_n; ++i)
        res += At(i, i);
    return res;
}

CMatrix CMatrix::Cofactor() const
{
    if(!SquareCheck())
        return CMatrix(0);

    CMatrix cofm(dim_n - 1), res(dim_n);

    for(size_t i = 0; i < dim_n; ++i)
    {
        for(size_t j = 0; j < dim_n; ++j)
        {
            // pick up an element (i, j) in the matrix
            size_t i2 = 0;
            for(size_t ii = 0; ii < dim_n; ++ii)
            {
                if(ii == i) continue; // on i row
                size_t j2 = 0;
                for(size_t jj = 0; jj < dim_n; ++jj)
                {
                    if(jj == j) continue; // on j column
                    cofm(i2, j2) = At(ii, jj); // fill minor matrix
                    ++j2;
                }
                ++i2;
            }
            // Cofactor matrix element
            res(i, j) = pow(-1.0, i+j) * cofm.Det();
        }
    }
    return res;
}

CMatrix CMatrix::Transpose() const
{
    CMatrix res(dim_m, dim_n);

    for(size_t i = 0; i < dim_n; ++i)
        for(size_t j = 0; j < dim_m; ++j)
            res(j, i) = elements[i][j];

    return res;
}

void CMatrix::TransposeSelf()
{
    if(!SquareCheck())
        return;

    for(size_t i = 1; i < dim_n; ++i)
    {
        for(size_t j = 0; j < i; ++j)
        {
            double tmp = elements[i][j];
            elements[i][j] = elements[j][i];
            elements[j][i] = tmp;
        }
    }
}

CMatrix CMatrix::Inverse() const
{
    double det = Det();

    if(det == 0.) {
        std::cerr << "CMatrix Error: This matrix's determinant is 0, no inverse matrix exists."
                  << std::endl;
        return CMatrix(0);
    }

    CMatrix res = Cofactor();
    res.TransposeSelf();

    return res/det;
}

CMatrix CMatrix::LowerTri() const
{
    if(!SquareCheck())
        return CMatrix(0);

    CMatrix res(dim_n);
    for(size_t i = 0; i < dim_n; ++i)
    {
        for(size_t j = 0; j <= i; ++j)
        {
            res(i, j) = At(i, j);
        }
    }
    return res;
}

CMatrix CMatrix::UpperTri() const
{
    if(!SquareCheck())
        return CMatrix(0);

    CMatrix res(dim_n);
    for(size_t i = 0; i < dim_n; ++i)
    {
        for(size_t j = i; j < dim_n; ++j)
        {
            res(i, j) = At(i, j);
        }
    }
    return res;
}

// return the Cholesky Decomposition Matrix L that follows M = LL^T
// TODO, add complex number support
CMatrix CMatrix::Cholesky() const
{
    if(!SquareCheck())
        return CMatrix(0);

    CMatrix res(dim_n);

    for(size_t i = 0; i < dim_n; ++i)
    {
        for(size_t j = 0; j < i + 1; ++j)
        {
            double s = 0;
            for(size_t k = 0; k < j; ++k)
            {
                s += res.At(i, k) * res.At(j, k);
            }
            if(i == j)
                res(i, j) = sqrt(At(i, i) - s);
            else
                res(i, j) = 1./res(j, j) * (At(i, j) - s);
        }
    }
    return res;
}

CMatrix cana::det(const CMatrix &m) {return m.Det();};
CMatrix cana::inv(const CMatrix &m) {return m.Inverse();};
CMatrix cana::transpose(const CMatrix &m) {return m.Transpose();};
CMatrix cana::trace(const CMatrix &m) {return m.Trace();};
CMatrix cana::cholesky(const CMatrix &m) {return m.Cholesky();};
CMatrix cana::tril(const CMatrix &m) {return m.LowerTri();};
CMatrix cana::triu(const CMatrix &m) {return m.UpperTri();};
CMatrix cana::power(const CMatrix &m, const int &n) {return m.Power(n);};
CMatrix cana::identity(const size_t &n) {CMatrix temp(n); return temp.Identity();};
