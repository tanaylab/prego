#ifndef PREGO_BLAS_PIN_H
#define PREGO_BLAS_PIN_H

#ifndef _WIN32
#include <dlfcn.h>
#endif

// Keep the BLAS serial INSIDE a RcppParallel worker thread.
//
// OpenBLAS built against OpenMP - which is what conda ships, despite the "p"
// in libopenblasp-*.so; it links libgomp and imports omp_get_max_threads -
// sizes its thread team from omp_get_max_threads() ON THE CALLING THREAD. The
// dgemm calls below run inside a RcppParallel worker, so without this each
// worker opens its own full-size OpenMP team: n_workers * omp_max OS threads,
// and the run collapses into scheduler churn.
//
// This has to happen HERE rather than from R. libgomp's thread-count ICV is
// per-thread, and a thread that has never entered a parallel region inherits
// the *global* ICV, which omp_set_num_threads() does not modify - so calling it
// on R's main thread leaves these workers untouched. Measured standalone
// against that exact OpenBLAS, 16 pthreads each running dgemm: no pin 54.91s
// at 2034 threads, pinned per worker 0.08s at 1 thread, identical to setting
// OMP_NUM_THREADS=1 before the process starts (0.09s).
//
// Resolved with dlsym rather than by linking OpenMP: prego does not compile
// with -fopenmp, and where no OpenMP runtime is loaded there is no team to pin,
// so a missing symbol is simply a no-op. Left pinned for the thread's lifetime;
// 1 is the correct value for a worker thread, and restoring it per call would
// only re-widen the next dgemm.
inline void pin_worker_blas_to_one_thread() {
#ifndef _WIN32
    typedef void (*omp_set_num_threads_fn)(int);
    static omp_set_num_threads_fn set_omp_threads =
        (omp_set_num_threads_fn) dlsym(RTLD_DEFAULT, "omp_set_num_threads");
    if (set_omp_threads != NULL) {
        set_omp_threads(1);
    }
#endif
}

// BLAS routine declaration
extern "C" {
    void dgemm_(const char* TRANSA, const char* TRANSB,
                const int* M, const int* N, const int* K,
                const double* ALPHA, const double* A, const int* LDA,
                const double* B, const int* LDB,
                const double* BETA, double* C, const int* LDC);
}

#endif // PREGO_BLAS_PIN_H
