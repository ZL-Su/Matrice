#include <private/math/_config.h>
#include <private/nonfree/inl/blas_kernel_wrapper.inl>
#include <private/nonfree/inl/lapack_kernel_wrapper.inl>

#if MATRICE_MATH_KERNEL==MATRICE_USE_MKL
DGE_MATRICE_BEGIN
using f32_t = float;
using f64_t = double;
using fptr = f32_t * ;
using dptr = f64_t * ;
_INTERNAL_BEGIN
template<>
decltype(auto) _blas_asum(const fptr x, size_t n, int inc) {
	return cblas_sasum(MKL_INT(n), x, MKL_INT(inc));
}
template<>
decltype(auto) _blas_asum(const dptr x, size_t n, int inc) {
	return cblas_dasum(MKL_INT(n), x, MKL_INT(inc));
}

template<>
decltype(auto) _blas_axpy(const f32_t a, const fptr x, fptr y, size_t n, int incx, int incy) {
	return cblas_saxpy(MKL_INT(n), a, x, MKL_INT(incx), y, MKL_INT(incy));
}
template<>
decltype(auto) _blas_axpy(const f64_t a, const dptr x, dptr y, size_t n, int incx, int incy) {
	return cblas_daxpy(MKL_INT(n), a, x, MKL_INT(incx), y, MKL_INT(incy));
}

template<>
decltype(auto) _blas_dot(size_t n, const fptr x, int incx, const fptr y, int incy) {
	return cblas_sdot(MKL_INT(n), x, MKL_INT(incx), y, MKL_INT(incy));
}
template<>
decltype(auto) _blas_dot(size_t n, const dptr x, int incx, const dptr y, int incy) {
	return cblas_ddot(MKL_INT(n), x, MKL_INT(incx), y, MKL_INT(incy));
}

template<>
decltype(auto) _blas_nrm2(size_t n, const fptr x, int inc) {
	return cblas_snrm2(MKL_INT(n), x, MKL_INT(inc));
}
template<>
decltype(auto) _blas_nrm2(size_t n, const dptr x, int inc) {
	return cblas_dnrm2(MKL_INT(n), x, MKL_INT(inc));
}

template<>
decltype(auto) _blas_scal(size_t n, const f32_t a, fptr x, int inc) {
	return cblas_sscal(MKL_INT(n), a, x, MKL_INT(inc));
}
template<>
decltype(auto) _blas_scal(size_t n, const f64_t a, dptr x, int inc) {
	return cblas_dscal(MKL_INT(n), a, x, MKL_INT(inc));
}

template<>
void _blas_gemv(size_t m, size_t n, f32_t a, const fptr A, const fptr x, int incx, f32_t b, fptr y, size_t incy) {
	return cblas_sgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, 
		MKL_INT(m), MKL_INT(n), a, A, MKL_INT(n > 1 ? n : 1), 
		x, MKL_INT(incx), b, y, MKL_INT(incy));
}
template<>
void _blas_gemv(size_t m, size_t n, f64_t a, const dptr A, const dptr x, int incx, f64_t b, dptr y, size_t incy) {
	return cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, 
		MKL_INT(m), MKL_INT(n), a, A, MKL_INT(n > 1 ? n : 1), 
		x, MKL_INT(incx), b, y, MKL_INT(incy));
}
template<>
void _blas_gemtv(size_t m, size_t n, f32_t a, const fptr A, const fptr x, int incx, f32_t b, fptr y, size_t incy) {
	return cblas_sgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans, MKL_INT(m), MKL_INT(n), a, A, MKL_INT(m > 1 ? m : 1), x, MKL_INT(incx), b, y, MKL_INT(incy));
}
template<>
void _blas_gemtv(size_t m, size_t n, f64_t a, const dptr A, const dptr x, int incx, f64_t b, dptr y, size_t incy) {
	return cblas_dgemv(CBLAS_LAYOUT::CblasRowMajor, CBLAS_TRANSPOSE::CblasTrans, MKL_INT(m), MKL_INT(n), a, A, MKL_INT(m > 1 ? m : 1), x, MKL_INT(incx), b, y, MKL_INT(incy));
}
template<>
void _blas_gemm(int lyt, int trpa, int trpb, size_t m, size_t n, size_t k, f32_t a, const fptr A, size_t lda, const fptr B, size_t ldb, f32_t b, fptr C, size_t ldc) {
	return cblas_sgemm(CBLAS_LAYOUT(lyt), CBLAS_TRANSPOSE(trpa), CBLAS_TRANSPOSE(trpb), MKL_INT(m), MKL_INT(n), MKL_INT(k), a, A, MKL_INT(lda), B, MKL_INT(ldb), b, C, MKL_INT(ldc));
}
template<>
void _blas_gemm(int lyt, int trpa, int trpb, size_t m, size_t n, size_t k, f64_t a, const dptr A, size_t lda, const dptr B, size_t ldb, f64_t b, dptr C, size_t ldc) {
	return cblas_dgemm(CBLAS_LAYOUT(lyt), CBLAS_TRANSPOSE(trpa), CBLAS_TRANSPOSE(trpb), MKL_INT(m), MKL_INT(n), MKL_INT(k), a, A, MKL_INT(lda), B, MKL_INT(ldb), b, C, MKL_INT(ldc));
}
_INTERNAL_END

_INTERNAL_BEGIN
template<>
lapack_int _lapack_potrf(int lyt, char ul, lapack_int n, float* a, lapack_int lda) {
	return LAPACKE_spotrf(lyt, ul, n, a, lda);
}
template<>
lapack_int _lapack_potrf(int lyt, char ul, lapack_int n, double* a, lapack_int lda) {
	return LAPACKE_dpotrf(lyt, ul, n, a, lda);
}
template<>
lapack_int _lapack_gesvd(int lyt, char jobu, char jobvt, lapack_int m, lapack_int n, fptr a, lapack_int lda, fptr s, fptr u, lapack_int ldu, fptr vt, lapack_int ldvt, fptr superb) {
	return LAPACKE_sgesvd(lyt, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}
template<>
lapack_int _lapack_gesvd(int lyt, char jobu, char jobvt, lapack_int m, lapack_int n, dptr a, lapack_int lda, dptr s, dptr u, lapack_int ldu, dptr vt, lapack_int ldvt, dptr superb) {
	return LAPACKE_dgesvd(lyt, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, superb);
}
template<>
lapack_int _lapack_syev(int lyt, char job, char ul, lapack_int n, fptr a, lapack_int lda, fptr w) {
	return LAPACKE_ssyev(lyt, job, ul, n, a, lda, w);
}
template<>
lapack_int _lapack_syev(int lyt, char job, char ul, lapack_int n, dptr a, lapack_int lda, dptr w) {
	return LAPACKE_dsyev(lyt, job, ul, n, a, lda, w);
}
_INTERNAL_END

DGE_MATRICE_END
#endif