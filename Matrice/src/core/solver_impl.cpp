
#include "../../include/Matrice/core/solver.h"
#include "../../include/Matrice/core/matrix.h"
#ifdef __use_mkl__
#include <mkl.h>
#else
#include <fkl.h>
#endif

namespace dgelom {namespace types {

template<typename _Ty> template<int _M, int _N, solver_type _alg>
void Solver_<_Ty>::Linear<_M, _N, _alg>::_solver_kernel(iterator b)
{
	const int_t n = A.cols() - 1;

	if constexpr (SolverType == solver_type::SVD)
	{
		if constexpr (type_bytes<value_t>::value == 4)
			flapk::_sgsvdsv(fkl::sptr(A.data()), fkl::sptr(b), A.cols(), A.rows());
		if constexpr (type_bytes<value_t>::value == 8)
			flapk::_dgsvdsv(fkl::dptr(A.data()), fkl::dptr(b), A.cols(), A.rows());
		return;
	}
}
template void Solver_<float>::Linear<0, 0, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<2, 2, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<3, 2, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<3, 3, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<3, 4, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<3, 5, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<3, 6, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<4, 2, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<4, 3, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<4, 4, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<4, 5, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<4, 6, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<5, 2, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<5, 3, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<5, 4, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<5, 5, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<5, 6, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<6, 2, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<6, 3, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<6, 4, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<6, 5, SVD>::_solver_kernel(iterator);
template void Solver_<float>::Linear<6, 6, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<0, 0, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<2, 2, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<3, 2, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<3, 3, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<3, 4, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<3, 5, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<3, 6, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<4, 2, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<4, 3, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<4, 4, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<4, 5, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<4, 6, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<5, 2, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<5, 3, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<5, 4, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<5, 5, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<5, 6, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<6, 2, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<6, 3, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<6, 4, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<6, 5, SVD>::_solver_kernel(iterator);
template void Solver_<double>::Linear<6, 6, SVD>::_solver_kernel(iterator);

template<typename _Ty> template<int _M, int _N, solver_type _alg>
void Solver_<_Ty>::Linear<_M, _N, _alg>::_solver_kernel(iterator b, int inc)
{
	const int_t n = A.cols() - 1;

	if constexpr (SolverType == solver_type::GLS)
	{
		Matrix_<value_t, __, __> A_view = A, temp_b(n + 1, 1);
		for (int_t i = 0; i < n; ++i) temp_b(i) = b[i*inc];
		temp_b = typename dgelom::solve::LinearOp::Gls(A_view, temp_b);
		for (int_t i = 0; i < n; ++i) b[i*inc] = temp_b(i);
		return;
	}
	if constexpr (SolverType == solver_type::SVD)
	{
		if constexpr (type_bytes<value_t>::value == 4)
			flapk::_sgsvdsv(fkl::sptr(A.data()), fkl::sptr(b), A.cols(), A.rows());
		if constexpr (type_bytes<value_t>::value == 8)
			flapk::_dgsvdsv(fkl::dptr(A.data()), fkl::dptr(b), A.cols(), A.rows());
		return;
	}

	for (int_t i = 1; i <= n; ++i) {
		value_t sum = value_t(0);
		for (int_t j = 0; j < i; ++j)
			sum += A[i][j] * b[j*inc];
		b[i*inc] -= sum;
	}

	b[n*inc] = b[n*inc] / A[n][n];
	if (this->options.used_alg == solver_type::CHD) {
		auto L = A.transpose();
		for (int_t i = n - 1; i >= 0; --i) {
			int_t k = i * (n + 1);
			value_t sum = value_t(0);
			for (int_t j = i + 1; j <= n; ++j)
				sum += L(j + k)*b[j*inc];
			b[i*inc] = (b[i*inc] - sum) / L(i + k);
		}
		return;
	}
	if (this->options.used_alg == solver_type::LUF) {
		for (int_t i = n - 1; i >= 0; --i) {
			value_t sum = value_t(0);
			for (int_t j = i + 1; j <= n; ++j)
				sum += A[i][j]*b[j*inc];
			b[i*inc] = (b[i*inc] - sum) / A[i][i];
		}
		return;
	}
}

template void Solver_<float>::Linear<0, 0, AUTO>::_solver_kernel(iterator, int);
template void Solver_<float>::Linear<2, 2, AUTO>::_solver_kernel(iterator, int);
template void Solver_<float>::Linear<3, 3, AUTO>::_solver_kernel(float*, int);
template void Solver_<float>::Linear<4, 4, AUTO>::_solver_kernel(float*, int);
template void Solver_<float>::Linear<5, 5, AUTO>::_solver_kernel(float*, int);
template void Solver_<float>::Linear<6, 6, AUTO>::_solver_kernel(float*, int);
template void Solver_<double>::Linear<0, 0, AUTO>::_solver_kernel(iterator, int);
template void Solver_<double>::Linear<2, 2, AUTO>::_solver_kernel(iterator, int);
template void Solver_<double>::Linear<3, 3, AUTO>::_solver_kernel(double*, int);
template void Solver_<double>::Linear<4, 4, AUTO>::_solver_kernel(double*, int);
template void Solver_<double>::Linear<5, 5, AUTO>::_solver_kernel(double*, int);
template void Solver_<double>::Linear<6, 6, AUTO>::_solver_kernel(double*, int);

template<typename _Ty> template<int _M, int _N, solver_type _alg>
void Solver_<_Ty>::Linear<_M, _N, _alg>::_Pre_solve()
{
	using Matrix = Matrix_<_M, _N>;
	typename Matrix::pointer pCoef = A.data(); 
	int layout = layout_traits<decltype(A)>::is_rmajor(A.format) ? rmaj : cmaj;

	if constexpr(SolverType == solver_type::AUTO && CompileTimeRows == CompileTimeCols)
	{
		if (layout_traits<decltype(A)>::is_symmetric(A.format)) {
			if constexpr (type_bytes<value_t>::value == 4)
#ifdef __use_mkl__
				options.status = LAPACKE_spotrf(layout, 'L', _N, (float*)pCoef, _N);
#else
				options.status = flapk::_scholy(fkl::sptr(pCoef), _N);
			
#endif
			if constexpr (type_bytes<value_t>::value == 8)
#ifdef __use_mkl__
				options.status = LAPACKE_dpotrf(layout, 'L', _N, (double*)pCoef, _N);
#else
				options.status = flapk::_dcholy(fkl::dptr(pCoef), _N);
#endif
			options.used_alg = solver_type::CHD;
		}
		else {
			if constexpr (type_bytes<value_t>::value == 4)
#ifdef __use_mkl__
				options.status = LAPACKE_sgetrf(layout, _M, _N, (float*)pCoef, _N, options.iwp.data());
#else
				options.status = flapk::_sLU(fkl::sptr(pCoef), options.iwp, _N);
#endif
			if constexpr (type_bytes<value_t>::value == 8)
#ifdef __use_mkl__
				options.status = LAPACKE_dgetrf(layout, _M, _N, (double*)pCoef, _N, options.iwp.data());
#else
				options.status = flapk::_dLU(fkl::dptr(pCoef), options.iwp, _N);
#endif
			options.used_alg = solver_type::LUF;
		}
	}
	if constexpr (CompileTimeRows != CompileTimeCols) {
		options.used_alg = solver_type::SVD;
	}
}

template void Solver_<float>::Linear<0, 0, solver_type::AUTO>::_Pre_solve();
template void Solver_<float>::Linear<2, 2, solver_type::AUTO>::_Pre_solve();
template void Solver_<float>::Linear<3, 3, solver_type::AUTO>::_Pre_solve();
template void Solver_<float>::Linear<4, 4, solver_type::AUTO>::_Pre_solve();
template void Solver_<float>::Linear<5, 5, solver_type::AUTO>::_Pre_solve();
template void Solver_<float>::Linear<6, 6, solver_type::AUTO>::_Pre_solve();
template void Solver_<double>::Linear<0, 0, solver_type::AUTO>::_Pre_solve();
template void Solver_<double>::Linear<2, 2, solver_type::AUTO>::_Pre_solve();
template void Solver_<double>::Linear<3, 3, solver_type::AUTO>::_Pre_solve();
template void Solver_<double>::Linear<4, 4, solver_type::AUTO>::_Pre_solve();
template void Solver_<double>::Linear<5, 5, solver_type::AUTO>::_Pre_solve();
template void Solver_<double>::Linear<6, 6, solver_type::AUTO>::_Pre_solve();

}}

MATRICE_NAMESPACE_BEGIN_TYPES

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

template<typename _T> void LinearOp::OpBase<_T>::_Impl(view_t& A, view_t& X)
{
	const int_t n = A.cols() - 1;
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
		auto L = A.transpose();
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

template<typename _T> LinearOp::info_t LinearOp::OpBase<_T>::_Impl(view_t& U, view_t& S, view_t& V)
{
	using fkl::sptr; using fkl::dptr;
	using value_t = typename view_t::value_t;
	if (S.empty) S.create(U.cols(), 1);
	if (V.empty) V.create(U.cols(), U.cols());

	info_t info;
	if constexpr (type_bytes<value_t>::value == 4)
		info.status = flapk::_sgesvd(sptr(U.data()), sptr(S.data()), sptr(V.data()), U.rows(), U.cols());
	if constexpr (type_bytes<value_t>::value == 8)
		info.status = flapk::_dgesvd(dptr(U.data()), dptr(S.data()), dptr(V.data()), U.rows(), U.cols());
	info.alg = solver_type::SVD;

	return info;
}
template LinearOp::info_t LinearOp::OpBase<float>::_Impl(Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&);
template LinearOp::info_t LinearOp::OpBase<double>::_Impl(Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&, Matrix_<value_t, __, __>&);

MATRICE_NAMESPACE_END_TYPES
