/*  *************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
*	*************************************************************************/
#pragma once
#include <algorithm>
#include <functional>
#include <numeric>
#include <future>
#include <type_traits>
#include "_svd.h"
#include "../_type_traits.h"
#include "util/_type_defs.h"
#include "util/utils.h"
#include "util/genalgs.h"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty, int _M, int _N> class Matrix_;

template<class _Derived> class SolverBase
{
	using derived_t = _Derived;
public:
	struct Options{ 
		int status;
		int iwp[max(1, derived_t::cols_at_compiletime)]; 
		solver_type used_alg = solver_type(derived_t::SolverType);
	};
	constexpr SolverBase() = default;

protected:
	template<typename... _Args> MATRICE_GLOBAL_FINL
	auto _Impl(_Args... args) { return static_cast<derived_t*>(this)->m_op(args...); }
};

struct LinearOp MATRICE_NONHERITABLE
{
	template<typename _Ty, int _M, int _N=_M>
	using Matrix_ = detail::Matrix_<_Ty, _M, _N>;
	struct info_t { solver_type alg = AUTO; int status = 1; int sign = 1; };
	template<typename _T> class OpBase 
	{
		MATRICE_GLOBAL_FINL std::string _Op_name(solver_type _type) const {
			switch (_type)
			{
			case dgelom::AUTO: return "Auto"; break; case dgelom::LUF: return "LU factorisation"; break;
			case dgelom::CHD: return "Cholesky decomposition"; break; case dgelom::QRD: return "QR decomposition"; break;
			case dgelom::SVD: return "SVD decomposition"; break; case dgelom::ITER: return "Iterator"; break;
			case dgelom::SPAR: return "Sparse"; break; case dgelom::GLS: return "Gaussian elimination"; break;
			default: return "Unknown"; break;
			}
		}
	public:
		using value_t = typename conditional <is_scalar_v<_T>, _T, default_type>::type;
		using view_t = Matrix_<value_t, ::dynamic>;

	protected:
		info_t _Info;
		view_t _Aview, _Bview, _Cview, _Dview;
		std::future<info_t> _Future;
		std::function<void()> _Launch = [&]{ 
			if (_Future.valid()) _Info = _Future.get(); 
			if (_Info.status) throw std::runtime_error("Solver kernel " + _Op_name(_Info.alg) + " error: " + std::to_string(_Info.status));
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
		enum { N = _Mty::cols_at_compiletime };
		enum { option = solver_type::AUTO };

		MATRICE_HOST_INL constexpr Auto(const _Mty& arg) : A(arg) {
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
		enum { N = _Mty::cols_at_compiletime };
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
		using _Mty = typename std::remove_reference<_T>::type;
		enum { N = _Mty::cols_at_compiletime };
		enum { option = solver_type::SVD };

		MATRICE_GLOBAL_FINL constexpr Svd(const _Mty& _Coeff) : U(_Coeff) {
			S.create(U.cols(), 1), Vt.create(U.cols(), U.cols());
			using Op = OpBase<value_t>;
			Op::_Aview.ref(U), Op::_Bview.ref(S), Op::_Cview.ref(Vt);
			Op::_Impl(Op::_Aview, Op::_Bview, Op::_Cview);
		};
		/**
		 * Return nullspace $X$ of the src matrix, s.t. $AX=0$.
		 */
		MATRICE_GLOBAL_FINL constexpr auto operator()()const noexcept {
			Matrix_<value_t, N, 1> X(Vt.cols(), 1);
			const size_t last = Vt.rows() - 1;
			const value_t* head = vt()[last];
			X.from(head);
			return std::forward<decltype(X)>(X);
		}
		template<typename _Ret = Matrix_<value_t, N, min(N, 1)>> 
		MATRICE_GLOBAL_FINL constexpr _Ret& operator() (_Ret& X){
			Matrix_<value_t, _Mty::cols_at_compiletime, min(_Mty::cols_at_compiletime, 1)> Z(X.rows(), 1); //inverse of singular values

			return (X);
		}
		//\return singular values
		MATRICE_HOST_FINL const auto& sv()const noexcept { return (S); }
		MATRICE_HOST_FINL auto& sv()noexcept { return (S); }
		//\return $V^T$ expression
		MATRICE_HOST_FINL const auto& vt()const noexcept { return (Vt); }
		MATRICE_HOST_FINL auto& vt()noexcept { return (Vt); }
	private:
		const _Mty& U;
		Matrix_<value_t, _Mty::cols_at_compiletime, min(_Mty::cols_at_compiletime,1)> S;
		Matrix_<value_t, _Mty::cols_at_compiletime, _Mty::cols_at_compiletime> Vt; //V^{T}
	};
	template<typename _T> struct Evv : public OpBase<typename _T::value_t>
	{
		using typename OpBase<typename _T::value_t>::value_t;
		using _Mty = typename std::remove_reference<_T>::type;
		enum { N = _Mty::cols_at_compiletime };
		enum { option = solver_type::SVD };

		constexpr Evv(const _Mty& args) : A(args) {}

	private:
		const _Mty& A;
	};
};
_DETAIL_END
DGE_MATRICE_END