#ifndef MELBLAS_OMP_HPP
#define MELBLAS_OMP_HPP

#include <omp.h>
#include "MelBLAS.hpp"

namespace melblas
{

    template<typename FloatingType>
    class MelBLAS_OMP :
    public MelBLAS_B<FloatingType>{

        public:
            FloatingType dot(const FloatingType * x, const FloatingType * y, size_t size) const;

            void axpby(FloatingType alpha, const FloatingType * x, FloatingType beta, FloatingType * y, size_t size) const;

            void gemv(FloatingType alpha, const FloatingType * A, const FloatingType * x, FloatingType beta, FloatingType * y, size_t num_rows, size_t num_cols) const;

            ~MelBLAS_OMP() = default;
    };


    template<typename FloatingType>
    FloatingType dot(const FloatingType * x, const FloatingType * y, size_t size) const
    {
        FloatingType result = 0.0;

        #pragma omp parallel for
        for(size_t i = 0; i < size; i++)
        {
            result += x[i] * y[i];
        }
        return result;
    }

    template<typename FloatingType>
    void axpby(FloatingType alpha, const FloatingType * x, FloatingType beta, FloatingType * y, size_t size) const
    {
        // y = alpha * x + beta * y

        #pragma omp parallel for
        for(size_t i = 0; i < size; i++)
        {
            y[i] = alpha * x[i] + beta * y[i];
        }
    }

    template<typename FloatingType>
    void gemv(FloatingType alpha, const FloatingType * A, const FloatingType * x, FloatingType beta, FloatingType * y, size_t num_rows, size_t num_cols) const;
    {

        // y = alpha * A * x + beta * y;
        #pragma omp parallel for
        for(size_t r = 0; r < num_rows; r++)
        {
            FloatingType y_val = 0.0;
            for(size_t c = 0; c < num_cols; c++)
            {
                y_val += alpha * A[r * num_cols + c] * x[c];
            }
            y[r] = beta * y[r] + y_val;
        }
    }



}


#endif