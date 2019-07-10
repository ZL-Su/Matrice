
#include "../../include/Matrice/core/solver.h"
#include "../../include/Matrice/core/matrix.h"
#include "../../include/Matrice/private/nonfree/_lnalge.h"
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif

DGE_MATRICE_BEGIN _DETAIL_BEGIN

/**
 *\Perform decomposition for coeff. matrix A.
 */
template<typename _T> LinearOp::info_t LinearOp::OpBase<_T>::_Impl(view_t& A)
{
	typename view_t::pointer pCoef = A.data();
	size_t _M = A.rows(), _N = A.cols();
	if (_M != _N) throw std::runtime_error("Support only for square matrix.");

	int layout = layout_traits<view_t>::is_rmajor(A.format) ? rmaj : cmaj;
	info_t info;

	if (layout_traits<view_t>::is_symmetric(A.format)) {
		info.alg = solver_type::CHD;
		if constexpr (type_bytes<value_t>::value == 4)
#ifdef __use_mkl__
			info.status = LAPACKE_spotrf(layout, 'L', _N, (float*)pCoef, _N);
#else
			info.status = flapk::_scholy(fkl::sptr(pCoef), _N);

#endif
		if constexpr (type_bytes<value_t>::value == 8)
#ifdef __use_mkl__
			info.status = LAPACKE_dpotrf(layout, 'L', _N, (double*)pCoef, _N);
#else
			info.status = flapk::_dcholy(fkl::dptr(pCoef), _N);
#endif
		return info;
	}

	{ //general dense matrix
		info.alg = solver_type::LUF;
		Matrix_<int, view_t::CompileTimeCols, min(view_t::CompileTimeCols, 1)> iwp(_N, 1);
		if constexpr (type_bytes<value_t>::value == 4)
#ifdef __use_mkl__
			info.status = LAPACKE_sgetrf(layout, _M, _N, (float*)pCoef, _N, iwp.data());
#else
			info.status = flapk::_sLU(fkl::sptr(pCoef), iwp.data(), _N);
#endif
		if constexpr (type_bytes<value_t>::value == 8)
#ifdef __use_mkl__
			info.status = LAPACKE_dgetrf(layout, _M, _N, (double*)pCoef, _N, iwp.data());
#else
			info.status = flapk::_dLU(fkl::dptr(pCoef), iwp.data(), _N);
#endif
		for (int i = 1; i <= _N; ++i) if (i != iwp(i)) info.sign *= -1;
		return info;
	}
}
template LinearOp::info_t LinearOp::OpBase<float>::_Impl(Matrix_<value_t, __, __>&);
template LinearOp::info_t LinearOp::OpBase<double>::_Impl(Matrix_<value_t, __, __>&);

/**
 *\Solve the solution X.
 */
template<typename _T> void LinearOp::OpBase<_T>::_Impl(view_t& A, view_t& X)
{
	const auto n = A.cols() - 1;
	const int_t inc = X.cols();

	if (this->_Info.alg == solver_type::GLS) return;

	auto _fn = [&](typename view_t::iterator b)->typename view_t::iterator {
		for (int_t i = 1; i <= n; ++i) {
			value_t sum = value_t(0);
			for (int_t j = 0; j < i; ++j) sum += A[i][j] * b[j*inc];
			b[i*inc] -= sum;
		}
		b[n*inc] = b[n*inc] / A[n][n];
		return b;
	};
	
	if (this->_Info.alg == solver_type::CHD) {
		auto L = A.t();
		for (int_t c = 0; c < inc; ++c) {
			auto b = _fn(X.begin() + c);
			for (int_t i = n - 1; i >= 0; --i) {
				int_t k = i * (n + 1);
				value_t sum = value_t(0);
				for (int_t j = i + 1; j <= n; ++j)
					sum += L(j + k)*b[j*inc];
				b[i*inc] = (b[i*inc] - sum) / L(i + k);
			}
		}
		return;
	}
	if (this->_Info.alg == solver_type::LUF) {
		for (int_t c = 0; c < inc; ++c) {
			auto b = _fn(X.begin() + c);
			for (int_t i = n - 1; i >= 0; --i) {
				value_t sum = value_t(0);
				for (int_t j = i + 1; j <= n; ++j)
					sum += A[i][j] * b[j*inc];
				b[i*inc] = (b[i*inc] - sum) / A[i][i];
			}
		}
		return;
	}
}
template void LinearOp::OpBase<float>::_Impl(Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&);
template void LinearOp::OpBase<double>::_Impl(Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&);

/**
 *\Perform SVD decomposition for U.
 */
template<typename _T> LinearOp::info_t LinearOp::OpBase<_T>::_Impl(view_t& U, view_t& S, view_t& Vt)
{
	using value_t = typename view_t::value_t;
	if (S.empty) S.create(U.cols(), 1);
	if (Vt.empty) Vt.create(U.cols(), U.cols());

	info_t info;
	info.status = detail::_Lapack_kernel_impl<value_t>::svd(
		U.data(), S.data(), Vt.data(), U.shape().tiled());
	info.alg = solver_type::SVD;

	return info;
}
template LinearOp::info_t LinearOp::OpBase<float>::_Impl(Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&);
template LinearOp::info_t LinearOp::OpBase<double>::_Impl(Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&);

_DETAIL_END DGE_MATRICE_END
