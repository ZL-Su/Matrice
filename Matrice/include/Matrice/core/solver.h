#pragma once
#include <algorithm>
#include <functional>
#include <numeric>
#include <future>
#include <type_traits>
#include "../private/_expr_type_traits.h"
#include "../util/_type_defs.h"
#include "../util/util.h"
#include "../util/genalgs.h"

MATRICE_NAMESPACE_BEGIN_TYPES

template<typename _Ty, int _M, int _N> class Matrix_;

template<class _Derived> class SolverBase
{
	using derived_t = _Derived;
public:
	struct Options{ int iwp[max(1, derived_t::CompileTimeCols)]; int status; solver_type used_alg = solver_type(derived_t::SolverType);};
	constexpr SolverBase() = default;

protected:
	template<typename... _Args>
	void _Prep(_Args... args) { static_cast<derived_t*>(this)->_Pre_solve(args...); };
	template<typename... _Args>
	void _Perf(_Args... args) { static_cast<derived_t*>(this)->_solver_kernel(args...); }
	template<typename... _Args>
	auto _Impl(_Args... args) { return static_cast<derived_t*>(this)->m_op(args...); }
};

struct LinearOp final
{
	struct info_t { solver_type alg = AUTO; int status = 1; int sign = 1; };
	template<typename _T> class OpBase 
	{
	public:
		using value_t = typename conditional <std::is_scalar<_T>::value, _T, default_type>::type;
	private:
		using view_t = Matrix_<value_t, __, __>;
	protected:
		info_t _Info;
		view_t _Aview, _Bview, _Cview, _Dview;
		std::future<info_t> _Future;
		std::function<void()> _Launch = [&]{ if (_Future.valid()) _Info = _Future.get(); };
		info_t _Impl(view_t& A);
		void   _Impl(view_t& A, view_t& x); 
		info_t _Impl(view_t& U, view_t& S, view_t& V);
	};
	template<typename _T> struct Auto : public OpBase<typename _T::value_t>
	{
		using base_t = OpBase<typename _T::value_t>;
		using typename base_t::value_t;
		using _Mty = typename std::remove_reference<_T>::type;
		enum { N = _Mty::CompileTimeCols };
		enum { option = solver_type::AUTO };

		constexpr Auto(const _Mty& arg) : A(arg) {
			OpBase<value_t>::_Future = std::async(std::launch::async, [&]()->info_t { 
				return OpBase<value_t>::_Impl(OpBase<value_t>::_Aview = A); });
		}

		template<typename _Ret = Matrix_<value_t, N, min(N, 1)>> constexpr _Ret& operator() (_Ret& X)
		{
			OpBase<value_t>::_Launch();
			OpBase<value_t>::_Impl(OpBase<value_t>::_Aview, OpBase<value_t>::_Bview = X);
			return (X);
		}
	private:
		const _Mty& A;
	};
	template<typename _T> struct Gls : public OpBase<typename _T::value_t>
	{
		using typename OpBase<typename _T::value_t>::value_t;
		using _Mty = typename std::remove_reference<_T>::type;
		enum { N = _Mty::CompileTimeCols };
		enum { option = solver_type::GLS };

		constexpr Gls(const _Mty& arg) : A(arg) { OpBase<value_t>::_Info.alg = option; };

		template<typename _Rhs> constexpr _Rhs& operator() (_Rhs& b) const { 
			return OpBase<value_t>::_Impl(OpBase<value_t>::_Aview = A, OpBase<value_t>::_Bview = b);
		};

	private:
		const _Mty& A;
	};
	template<typename _T> struct Svd : public OpBase<typename _T::value_t>
	{
		using typename OpBase<typename _T::value_t>::value_t;
		using _Mtx = typename std::remove_reference<_T>::type;
		enum { N = _Mtx::CompileTimeCols };
		enum { option = solver_type::SVD };

		constexpr Svd(const _Mtx& args) : U(args) {
			S.create(U.cols(), 1), V.create(U.cols(), U.cols());
			OpBase<value_t>::_Future = std::async(std::launch::async, [&] {
				using Op = OpBase<value_t>;
				return OpBase<value_t>::_Impl(Op::_Aview = U, Op::_Bview = S, Op::_Cview = V);
			});
		};
		constexpr auto operator() (std::_Ph<0> _ph = {}) {
			OpBase<value_t>::_Launch();
			Matrix_<value_t, _Mtx::CompileTimeCols, min(_Mtx::CompileTimeCols, 1)> X(U.cols(), 1);
			transform([&](value_t val)->value_t{return (val);}, V.begin(), V.end(), V.cols(), X.begin());
			return (X);
		}
		template<typename _Ret = Matrix_<value_t, N, min(N, 1)>> constexpr _Ret& operator() (_Ret& X){
			OpBase<value_t>::_Launch();

			Matrix_<value_t, _Mtx::CompileTimeCols, min(_Mtx::CompileTimeCols, 1)> Z(X.rows(), 1); //inverse of singular values

			return (X);
		}
		//\return singular values
		MATRICE_HOST_FINL auto& sv() { OpBase<value_t>::_Launch(); return (S); }
		//\return V^{T} expression
		MATRICE_HOST_FINL auto vt() { OpBase<value_t>::_Launch(); return (V.transpose()); }
	private:
		const _Mtx& U;
		Matrix_<value_t, _Mtx::CompileTimeCols, min(_Mtx::CompileTimeCols,1)> S;
		Matrix_<value_t, _Mtx::CompileTimeCols, _Mtx::CompileTimeCols> V; //V not V^{T}
	};
	template<typename _T> struct Evv : public OpBase<typename _T::value_t>
	{
		using typename OpBase<typename _T::value_t>::value_t;
		using _Mty = typename std::remove_reference<_T>::type;
		enum { N = _Mty::CompileTimeCols };
		enum { option = solver_type::SVD };

		constexpr Evv(const _Mty& args) : A(args) {}

	private:
		const _Mty& A;
	};
};

struct Solver final
{
	template<typename _Op = LinearOp::Auto<Matrix<default_type>>> class Linear_ : public SolverBase<Linear_<_Op>>
	{
		using Base = SolverBase<Linear_<_Op>>;
		using value_t = typename _Op::value_t;
		using typename Base::Options;
	public:
		_Op m_op;
		template<typename... _Args> constexpr Linear_(const _Args&... args) : m_op(args...) {};
		template<typename... _Args> constexpr auto solve(const _Args&... args) { return Base::_Impl(args...); }
	};
};

template<typename _Ty> class Solver_ final
{
public:
	const int dynamic = 0;
	using default_type = double;
	using value_t = typename conditional<std::is_scalar<_Ty>::value, _Ty, default_type>::type;
	template<int _M, int _N> using Matrix_ = Matrix_<value_t, _M, _N>;

	template<int _M, int _N, solver_type _alg = solver_type::AUTO>
	class Linear : public SolverBase<Linear<_M, _N, _alg>>
	{
		using Base = SolverBase<Linear<_M, _N, _alg>>;
		using typename Base::Options;
		using iterator = typename Matrix_<_M, _N>::iterator;
	public:
		using op_t = std::plus<value_t>;
		enum { SolverType = _alg };
		enum { CompileTimeRows = _M, CompileTimeCols = _N };
		constexpr Linear() = default;
		constexpr Linear(const Matrix_<_M, _N>& coef) : A(coef) { _Pre_solve(); }

		auto operator()(const Matrix_<CompileTimeRows, CompileTimeCols>& coef)
		{
			Matrix_<CompileTimeCols, min(CompileTimeCols, 1)> x(CompileTimeCols, 1);
			A = coef; Base::_Perf(x.begin());
			return (x);
		}
		
		template<int ComileTimeNrhs> auto& solve(Matrix_<CompileTimeCols, ComileTimeNrhs>& b)
		{
			for (std::ptrdiff_t offset = 0; offset < ComileTimeNrhs; ++offset) 
				Base::_Perf(b.begin() + offset, ComileTimeNrhs);
			return (b);
		}
		void _solver_kernel(iterator x);
		void _solver_kernel(iterator x, int inc);
	private:
		void _Pre_solve();
		Options options;
		Matrix_<CompileTimeRows, CompileTimeCols> A;
	};
};
MATRICE_NAMESPACE_END_TYPES