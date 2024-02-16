#ifndef CONJUGATEGRADIENT_HPP
#define CONJUGATEGRADIENT_HPP

#include <memory>
#include "MelBLAS/include/MelBLAS.hpp"

using namespace melblas;

template<typename FloatingType>
class ConjugateGradient{
    public:
        ConjugateGradient(std::unique_ptr<MelBLAS_B<FloatingType>>&& melblas): _melblas(std::move(melblas))
        {
            static_assert(std::is_floating_point<FloatingType>::value, "DataType must be floating point");
        };
        
        void solve(const FloatingType * A, const FloatingType * b, FloatingType * x, size_t size, int max_iters, FloatingType rel_error) const;

    private:    
        std::unique_ptr<MelBLAS_B<FloatingType>> _melblas;
};

template<typename FloatingType>
void ConjugateGradient<FloatingType>::solve(const FloatingType * A, const FloatingType * b, FloatingType * x, size_t size, int max_iters, FloatingType rel_error) const
{
    FloatingType alpha, beta, bb, rr, rr_new;
    FloatingType * r = new FloatingType[size];
    FloatingType * p = new FloatingType[size];
    FloatingType * Ap = new FloatingType[size];
    int num_iters;

    for(size_t i = 0; i < size; i++)
    {
        x[i] = 0.0;
        r[i] = b[i];
        p[i] = b[i];
    }

    bb = _melblas->dot(b, b, size);
    rr = bb;
    for(num_iters = 1; num_iters <= max_iters; num_iters++)
    {
        _melblas->gemv(1.0, A, p, 0.0, Ap, size, size);
        alpha = rr / _melblas->dot(p, Ap, size);
        _melblas->axpby(alpha, p, 1.0, x, size);
        _melblas->axpby(-alpha, Ap, 1.0, r, size);
        rr_new = _melblas->dot(r, r, size);
        beta = rr_new / rr;
        rr = rr_new;
        if(std::sqrt(rr / bb) < rel_error) { break; }
        _melblas->axpby(1.0, r, beta, p, size);
    }

    delete[] r;
    delete[] p;
    delete[] Ap;

    if(num_iters <= max_iters)
    {
        printf("Converged in %d iterations, relative error is %e\n", num_iters, std::sqrt(rr / bb));
    }
    else
    {
        printf("Did not converge in %d iterations, relative error is %e\n", max_iters, std::sqrt(rr / bb));
    }
}

#endif