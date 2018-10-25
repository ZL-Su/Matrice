#pragma once
#include <algorithm>
#include <functional>
#include <numeric>
#include <future>
#include <type_traits>
#include "../private/_type_traits.h"
#include "../private/math/_svd.h"
#include "../util/_type_defs.h"
#include "../util/utils.h"
#include "../util/genalgs.h"

DGE_MATRICE_BEGIN
_TYPES_BEGIN
template<typename _Ty, int _M, int _N> class Matrix_;

template<class _Derived> class SolverBase
{
	using derived_t = _Derived;
public:
	struct Options{ 
		int status;
		int iwp[max(1, derived_t::CompileTimeCols)]; 
		solver_type used_alg = solver_type(derived_t::SolverType);
	};
	constexpr SolverBase() = default;

protected:
	template<typename... _Args> MATRICE_GLOBAL_FINL
	auto _Impl(_Args... args) { return static_cast<derived_t*>(this)->m_op(args...); }
};

struct LinearOp MATRICE_NONHERITABLE
{
	struct info_t { solver_type alg = AUTO; int status = 1; int sign = 1; };
	template<typename _T> class OpBase 
	{
	public:
		using value_t = typename conditional <std::is_scalar<_T>::value, _T, default_type>::type;
		using view_t = Matrix_<value_t, __, __>;

	protected:
		info_t _Info;
		view_t _Aview, _Bview, _Cview, _Dview;
		std::future<info_t> _Future;
		std::function<void()> _Launch = [&]{ 
			if (_Future.valid()) _Info = _Future.get(); 
			if (!_Info.status) throw std::runtime_error("Solver kernel " + std::to_string(_Info.alg) + " error: " + std::to_string(_Info.status));
		};
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

		template<typename _Ret = Matrix_<value_t, N, min(N, 1)>> 
		MATRICE_HOST_INL constexpr _Ret& operator() (_Ret& X) {
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

		MATRICE_GLOBAL_FINL constexpr Gls(const _Mty& arg) : A(arg) {
			OpBase<value_t>::_Info.alg = option; 
		};

		template<typename _Rhs> 
		MATRICE_HOST_INL constexpr _Rhs& operator() (_Rhs& b) const {
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

		MATRICE_GLOBAL_FINL constexpr Svd(const _Mtx& _Coeff) : U(_Coeff) {
			S.create(U.cols(), 1), Vt.create(U.cols(), U.cols());
			OpBase<value_t>::_Future = std::async(std::launch::async, [&] {
				using Op = OpBase<value_t>;
				return OpBase<value_t>::_Impl(Op::_Aview = U, Op::_Bview = S, Op::_Cview = Vt);
			});
		};
		MATRICE_GLOBAL_FINL constexpr auto operator() (std::_Ph<0> _ph = {}) {
			OpBase<value_t>::_Launch();
			auto X = Vt.rview(Vt.rows() - 1).eval<_Mtx::CompileTimeCols>().transpose().eval();
			return std::forward<decltype(X)>(X);
		}
		template<typename _Ret = Matrix_<value_t, N, min(N, 1)>> 
		MATRICE_GLOBAL_FINL constexpr _Ret& operator() (_Ret& X){
			OpBase<value_t>::_Launch();

			Matrix_<value_t, _Mtx::CompileTimeCols, min(_Mtx::CompileTimeCols, 1)> Z(X.rows(), 1); //inverse of singular values

			return (X);
		}
		//\return singular values
		MATRICE_HOST_FINL auto& sv() { OpBase<value_t>::_Launch(); return (S); }
		//\return V^{T} expression
		MATRICE_HOST_FINL auto vt() { OpBase<value_t>::_Launch(); return (Vt); }
	private:
		const _Mtx& U;
		Matrix_<value_t, _Mtx::CompileTimeCols, min(_Mtx::CompileTimeCols,1)> S;
		Matrix_<value_t, _Mtx::CompileTimeCols, _Mtx::CompileTimeCols> Vt; //V^{T}
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

struct Solver MATRICE_NONHERITABLE
{
	template<typename _Op = LinearOp::Auto<Matrix_<default_type, __, __>>> 
	class Linear_ : public SolverBase<Linear_<_Op>>
	{
		using Base = SolverBase<Linear_<_Op>>;
		using value_t = typename _Op::value_t;
		using typename Base::Options;
	public:
		_Op m_op;
		template<typename... _Args> MATRICE_GLOBAL_FINL
		constexpr Linear_(const _Args&... args) : m_op(args...) {};
		template<typename... _Args> MATRICE_GLOBAL_FINL
		constexpr auto solve(const _Args&... args) { return Base::_Impl(args...); }
	};
};
_TYPES_END


/**
 *\Linear solver, default _Op is auto-solver-kernel.
 */
template<typename _Op = types::LinearOp::Auto<types::Matrix_<default_type, __, __>>> 
using linear_solver = types::Solver::Linear_<_Op>;

DGE_MATRICE_END