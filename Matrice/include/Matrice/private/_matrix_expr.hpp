/**************************************************************************
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
**************************************************************************/
#pragma once
#include <functional>
#include <cassert>
#ifdef __AVX__
#include <immintrin.h>
#endif
#include "../util/_macros.h"
#include "_expr_type_traits.h"
#include "_storage.hpp"

#pragma warning(disable: 4715)

MATRICE_NAMESPACE_EXPR_BEGIN

template<typename _Ty> class Matrix;
template<typename _Ty, int _Rows, int _cols> class Matrix_;
template<class _Lhs, class _Rhs, typename _BinaryOp> class MatBinaryExpr;
template<class _Oprd, typename _UnaryOp> class MatUnaryExpr;

struct Expr {
	using size_t = std::size_t;
	using default_type = double;
	enum OpFlag { ewise = 0, mmul = 1, inv = 2, trp = 3, undef = -1 };

	/*factorial_t<N> = N!*/
	template<int N> struct factorial_t { enum { value = factorial_t<N - 1>::value*N }; };
	template<> struct factorial_t<0> { enum { value = 1 }; };

	template<int _Option> struct option {enum{value = _Option | expr}; };

	struct Op
	{
		template<typename _Ty> using EwiseSum = std::plus<_Ty>;
		template<typename _Ty> using EwiseMin = std::minus<_Ty>;
		template<typename _Ty> using EwiseMul = std::multiplies<_Ty>;
		template<typename _Ty> using EwiseDiv = std::divides<_Ty>;
		template<typename _Ty> struct MatInv { enum { flag = inv }; MATRICE_GLOBAL _Ty* operator()(int M, _Ty* Out, _Ty* In = nullptr) const; };
		template<typename _Ty> struct MatTrp { enum { flag = trp }; MATRICE_HOST_FINL _Ty operator()(int M, _Ty* Out, _Ty* i) const { return (Out[int(i[1])*M + int(i[0])]); } };
		template<typename _Ty> struct SpreadMul
		{
			enum { flag = undef };
			template<typename _Lhs, typename _Rhs> MATRICE_GLOBAL_FINL
				_Ty operator() (const _Lhs& lhs, const _Rhs& rhs, int r, int c) const { return (lhs(r) * rhs(c)); }
		};
		template<typename _Ty> struct MatMul
		{
			enum { flag = mmul };
			template<typename _Rhs> MATRICE_GLOBAL_FINL
				_Ty operator() (const _Ty* lhs, const _Rhs& rhs, int c, int _plh = 0) const
			{
				_Ty val = _Ty(0);
				const int K = rhs.rows(), N = rhs.cols();
#ifdef __disable_simd__
				for (int k = 0; k < K; ++k) val += lhs[k] * rhs(k*N + c);
#else
#ifdef __AVX__
				for (int k = 0; k < K; ++k) val += lhs[k] * rhs(k*N + c);
#endif
#endif
				return (val);
			}
			template<typename _Lhs, typename _Rhs> MATRICE_GLOBAL_FINL
				_Ty operator() (const _Lhs& lhs, const _Rhs& rhs, int r, int c) const
			{
				_Ty _Ret = _Ty(0), val = 0;
				const int K = rhs.rows(), N = rhs.cols(), _Idx = r * lhs.cols();
#ifdef __disable_simd__
				for (int k = 0; k < K; ++k) {
					_Ret += lhs(_Idx + k) * rhs(k*N + c);
				}
#else
#ifdef __AVX__
				for (int k = 0; k < K; ++k) _Ret += lhs(_Idx + k) * rhs(k*N + c);
#endif
#endif
				return (_Ret);
			}
		};
	};

	template<typename _Op> struct Base_
	{
		using Derived = _Op;
		//using value_type = typename dgelom::conditonal<std::is_class<_Op>::value, typename _Op::value_t, default_type>::type;
		template<typename _Ty> MATRICE_GLOBAL_FINL operator Matrix<_Ty>()
		{
			Matrix<_Ty> ret(static_cast<Derived*>(this)->rows(), static_cast<Derived*>(this)->cols());
			return (ret = *static_cast<Derived*>(this));
		}
		template<typename _Ty, int _M, int _N> MATRICE_GLOBAL_FINL operator Matrix_<_Ty, _M, _N>()
		{
			Matrix<_Ty> ret; return (ret = *static_cast<Derived*>(this));
		}
		template<typename _OutType> MATRICE_GLOBAL_INL _OutType eval() const
		{
			_OutType ret(static_cast<const Derived*>(this)->rows(), static_cast<const Derived*>(this)->cols());
			return (ret = *static_cast<const Derived*>(this));
		}
		template<typename _OutType> MATRICE_GLOBAL_INL _OutType& assign(_OutType& res) const
		{
#ifdef _DEBUG
			if (res.size() == 0) res.create(M, N);
#endif
			static_cast<const Derived*>(this)->assign_to(res);
			return (res);
		}
		template<typename _OutType> MATRICE_GLOBAL_INL void synch(_OutType& res) const
		{
#ifdef _DEBUG
			if (res.size() == 0) res.create(M, N);
#endif
			static_cast<const Derived*>(this)->assign(res);
		}
		template<typename _Rhs> MATRICE_GLOBAL_INL
		MatBinaryExpr<Derived, _Rhs, Op::MatMul<typename _Rhs::value_t>> mul(const _Rhs& _rhs)
		{
			return MatBinaryExpr<Derived, _Rhs, Op::MatMul<typename _Rhs::value_t>>(*static_cast<Derived*>(this), _rhs);
		}
		MATRICE_GLOBAL_INL double sum() const
		{
			double val = 0.0;
			int n = static_cast<const Derived*>(this)->size();
			for (int i = 0; i < n; ++i) val += static_cast<const Derived*>(this)->operator()(i);
			return (val);
		}
		MATRICE_GLOBAL_FINL size_t size() const { return M*N; }
		MATRICE_GLOBAL_FINL size_t rows() const { return M; }
		MATRICE_GLOBAL_FINL size_t cols() const { return N; }
	protected:
		size_t M, K, N;
	};

	template<typename T, typename U, typename _BinaryOp> class EwiseBinaryExpr 
		: public Base_<EwiseBinaryExpr<T, U, _BinaryOp>>
	{
	public:
		enum { options = option<ewise>::value };
		using Lopd_t = T;
		using Ropd_t = U;
		using value_t = conditional_t<std::is_scalar_v<T>, T, typename T::value_t>;
		static constexpr value_t inf = std::numeric_limits<value_t>::infinity();

		MATRICE_GLOBAL_INL EwiseBinaryExpr(const value_t _scalar, const Ropd_t& _rhs) noexcept
			: _Scalar(_scalar), _LHS(_rhs), _RHS(_rhs) { M = _LHS.rows(), N = _RHS.cols(); }
		MATRICE_GLOBAL_INL EwiseBinaryExpr(const Lopd_t& _lhs, const value_t _scalar) noexcept
			: _Scalar(_scalar), _LHS(_lhs), _RHS(_lhs) { M = _LHS.rows(), N = _RHS.cols(); }
		MATRICE_GLOBAL_INL EwiseBinaryExpr(const Lopd_t& _lhs, const Ropd_t& _rhs) noexcept
			: _LHS(_lhs), _RHS(_rhs) { M = _LHS.rows(), N = _RHS.cols(); }

		MATRICE_GLOBAL_INL value_t operator() (size_t _idx)
		{
			if (_Scalar == inf) return _Op(_LHS(_idx), _RHS(_idx));
			if (_Scalar != inf) return _Op(_Scalar, _RHS(_idx));
		}
		MATRICE_GLOBAL_INL const value_t operator() (size_t _idx) const
		{
			if (_Scalar == inf) return _Op(_LHS(_idx), _RHS(_idx));
			if (_Scalar != inf) return _Op(_Scalar, _RHS(_idx));
		}
		template<typename _Mty> MATRICE_GLOBAL_INL void assign_to(_Mty& res) const
		{
			for (size_t i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}

	private:
		const value_t _Scalar = std::numeric_limits<value_t>::infinity();
		const Lopd_t& _LHS; const Ropd_t& _RHS;
		//std::function<value_t(value_t, value_t)> _Op = _BinaryOp();
		_BinaryOp _Op;
		using Base_<EwiseBinaryExpr<T, U, _BinaryOp>>::M;
		using Base_<EwiseBinaryExpr<T, U, _BinaryOp>>::N;
	};

	template<class _T1, class _T2, typename _BinaryOp> class MatBinaryExpr 
		: public Base_<MatBinaryExpr<_T1, _T2, _BinaryOp>>
	{
	public:
		enum{options = option<_BinaryOp::flag>::value};
		using Lopd_t = _T1; using Ropd_t = _T2;
		using value_t = typename dgelom::conditonal<std::is_class<_T1>::value, typename _T1::value_t, default_type>::type;

		MATRICE_GLOBAL_INL MatBinaryExpr(const Lopd_t& _lhs, const Ropd_t& _rhs) noexcept 
			: _LHS(_lhs), _RHS(_rhs) {
			M = _LHS.rows(), K = _LHS.cols(), N = _RHS.cols();
			if (K != _RHS.rows() && K == N) N = _RHS.rows();
		}

		MATRICE_HOST_INL value_t* operator() (value_t* res) const
		{
			return _Op(_LHS.data(), _RHS.data(), res, M, K, N);
		}
		MATRICE_GLOBAL_INL value_t operator() (int r, int c) const
		{
			return _Op(_LHS, _RHS, r, c);
		}
		MATRICE_GLOBAL_INL value_t operator() (int _idx) const
		{
			int r = _idx / N, c = _idx - r * N;
			return _Op(_LHS, _RHS, r, c);
		}
		template<typename _Mty> MATRICE_GLOBAL_INL void assign_to(_Mty& res) const
		{
			for (int i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}
	private:
		const Lopd_t& _LHS; 
		const Ropd_t& _RHS;
		_BinaryOp _Op;
		using Base_<MatBinaryExpr<_T1, _T2, _BinaryOp>>::M;
		using Base_<MatBinaryExpr<_T1, _T2, _BinaryOp>>::K;
		using Base_<MatBinaryExpr<_T1, _T2, _BinaryOp>>::N;
	};

	template<class _Oprd, typename _UnaryOp> class MatUnaryExpr 
		: public Base_<MatUnaryExpr<_Oprd, _UnaryOp>>
	{
	public:
		enum {options = option<_UnaryOp::flag>::value};
		using value_t = typename dgelom::conditonal<std::is_class<_Oprd>::value, typename _Oprd::value_t, typename dgelom::conditonal<std::is_scalar<_Oprd>::value, _Oprd, default_type>::type>::type;
		MATRICE_GLOBAL_FINL MatUnaryExpr(const _Oprd& inout) noexcept
			: _RHS(inout), _ANS(inout) {
			M = _RHS.rows(), N = _RHS.cols();
			if constexpr (options & trp == trp) std::swap(M, N);
		}
		MATRICE_GLOBAL_FINL MatUnaryExpr(const _Oprd& _rhs, _Oprd& _ans)
			: _RHS(_rhs), _ANS(_ans) {
			M = _RHS.rows(), N = _RHS.cols();
			if constexpr (options & trp == trp) std::swap(M, N);
		}
		MATRICE_GLOBAL_FINL MatUnaryExpr(MatUnaryExpr&& _other)
			: options(_other.options) {
			*this = std::move(_other);
		}

		MATRICE_GLOBAL_FINL auto operator()() const 
		{ 
			return _Op(M, _ANS.data(), _RHS.data());
		}
		MATRICE_GLOBAL_FINL value_t operator() (int r, int c) const 
		{ 
			if constexpr (options & inv == inv) return _ANS(c + r * N);
			if constexpr (options & trp == trp) return _ANS(r + c * N);
		}
		MATRICE_GLOBAL_FINL value_t operator() (int _idx) const 
		{ 
			if constexpr (options & trp == trp)
				_idx = _idx * M + _idx / N * (1 - N * M);
			return _ANS(_idx);
		}
		MATRICE_GLOBAL_FINL _Oprd& ans() { return _ANS; }
		MATRICE_GLOBAL_FINL const _Oprd& ans() const { return _ANS; }
		template<typename _Rhs, typename _Ret = MatBinaryExpr<MatUnaryExpr, _Rhs, Op::MatMul<value_t>>> MATRICE_GLOBAL_INL _Ret mul(const _Rhs& _rhs)
		{
			if constexpr (options & inv == inv) this->operator()();
			if constexpr (options & trp == trp)                   ;
			return _Ret(*this, _rhs);
		}
		template<typename _Mty> MATRICE_GLOBAL_INL void assign_to(_Mty& res) const
		{
			if constexpr (options & inv == inv) {
				_Op(M, res.data(), _RHS.data());
			}
			if constexpr (options & trp == trp) {
				for (int i = 0; i < size(); ++i) res(i) = (*this)(i);
			}
		} /*i*N+(1-N*M)*(i/M) is replaced by i at left hand*/
	private:
		const _Oprd& _RHS;
		const _Oprd& _ANS;
		_UnaryOp _Op;
		using Base_<MatUnaryExpr<_Oprd, _UnaryOp>>::M;
		using Base_<MatUnaryExpr<_Oprd, _UnaryOp>>::N;
		using Base_<MatUnaryExpr<_Oprd, _UnaryOp>>::size;
	};
};

MATRICE_NAMESPACE_EXPR_END