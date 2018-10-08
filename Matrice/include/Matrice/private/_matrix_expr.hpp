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
#include "_type_traits.h"
//#include "_storage.hpp"

#pragma warning(disable: 4715)

MATRICE_NAMESPACE_BEGIN_
// \forward declaration 
_TYPES_BEGIN
template<typename _Ty> class Matrix;
template<typename _Ty, int _Rows, int _cols> class Matrix_;
template<typename _Derived, typename _Traits, typename _Ty> class Base_;
_TYPES_END

template<typename _Ty, int _M, int _N> 
struct is_matrix <types::Matrix_<_Ty, _M, _N>> : std::true_type {};

template<typename _Derived, typename _Traits, typename _Ty> 
struct is_matrix<types::Base_<_Derived, _Traits, _Ty>> : std::true_type {};

template<typename _Ty, int _Rows, int _Cols> 
struct matrix_traits<types::Matrix_<_Ty, _Rows, _Cols>> {
	using type = _Ty;
	enum { _M = _Rows, _N = _Cols };
	struct size { struct rows { enum { value = _M }; }; struct cols { enum { value = _N }; }; };
	static constexpr bool Is_base = std::false_type::value;
};

template<typename _Derived, typename _Traits, typename _Ty>
struct matrix_traits<types::Base_<_Derived, _Traits, _Ty>> {
	using type = _Ty;
	enum {_M = _Traits::_M, _N = _Traits::_N};
	struct size { struct rows { enum { value = _M }; }; struct cols { enum { value = _N }; }; };
	static constexpr bool Is_base = std::true_type::value;
};
_MATRICE_NAMESPACE_END

MATRICE_NAMESPACE_EXPR_BEGIN
template<class _Lhs, class _Rhs, typename _BinaryOp> class MatBinaryExpr;
template<class _Lhs, typename _UnaryOp> class MatUnaryExpr;

struct Expr {
	using size_t = std::size_t;
	using default_type = double;
	enum OpFlag { ewise = 0, mmul = 1, inv = 2, trp = 3, undef = -1 };
	template<int _Option> struct option { enum { value = _Option | expr }; };

	/*factorial_t<N> = N!*/
	template<int N> struct factorial_t { enum { value = factorial_t<N - 1>::value*N }; };
	template<> struct factorial_t<0> { enum { value = 1 }; };

	struct Op
	{
		template<typename _Ty> using EwiseSum = std::plus<_Ty>;
		template<typename _Ty> using EwiseMin = std::minus<_Ty>;
		template<typename _Ty> using EwiseMul = std::multiplies<_Ty>;
		template<typename _Ty> using EwiseDiv = std::divides<_Ty>;
		template<typename _Ty> struct EwiseSqrt
		{
			enum { flag = ewise };
			MATRICE_GLOBAL_FINL auto operator() (const _Ty& _Val) const {
				return (std::sqrt(_Val));
			}
		};
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

	template<typename T, typename U> struct _Base_v1 {
		using matrix_type = conditional_t<is_matrix_v<T>, T, U>;
		using value_type = typename matrix_type::value_t;
		enum {
			_M = matrix_traits<matrix_type>::_M,
			_N = conditional_size_v<is_matrix_v<U>, matrix_traits<U>::_N, matrix_traits<matrix_type>::_N>
		};
	};
	template<typename _Op> struct Base_
	{
#define _CDTHIS static_cast<const_derived_pointer>(this)

		using Derived = _Op;
		using derived_pointer = std::add_pointer_t<Derived>;
		using const_derived = std::add_const_t<Derived>;
		using const_derived_pointer = std::add_pointer_t<const_derived>;
		//enum {_M = expression_traits<Derived>::_M, _N = expression_traits<Derived>::_N };

		/*template<typename _Ty> MATRICE_GLOBAL_FINL operator Matrix<_Ty>() {
			Matrix<_Ty> ret(static_cast<Derived*>(this)->rows(), static_cast<Derived*>(this)->cols());
			return (ret = *static_cast<Derived*>(this));
		}
		template<typename _Ty, int _M, int _N> MATRICE_GLOBAL_INL operator Matrix_<_Ty, _M, _N>() {
			Matrix_<_Ty, _M, _N> ret; 
			return (ret = *static_cast<Derived*>(this));
		}*/

		// \evaluate the matrix expression
		MATRICE_GLOBAL_INL auto eval() const { return _CDTHIS->eval_impl(); }
		template<typename _OutType> MATRICE_GLOBAL_INL _OutType eval() const {
			_OutType ret(_CDTHIS->rows(), _CDTHIS->cols());
			return std::forward<_OutType>(ret = *_CDTHIS);
		}
		// \retrieve the evaluated result
		template<typename _OutType> MATRICE_GLOBAL_INL _OutType& assign(_OutType& res) const {
#ifdef _DEBUG
			if (res.size() == 0) res.create(M, N);
#endif
			_CDTHIS->assign_to(res);
			return (res);
		}
		// \retrieve the evaluated result
		template<typename _OutType> MATRICE_GLOBAL_INL void synch(_OutType& res) const {
#ifdef _DEBUG
			if (res.size() == 0) res.create(M, N);
#endif
			_CDTHIS->assign(res);
		}
		// \formulate matmul expression
		template<typename _Rhs> MATRICE_GLOBAL_INL auto mul(const _Rhs& _rhs) const {
			return MatBinaryExpr<Derived, _Rhs, Op::MatMul<typename _Rhs::value_t>>(*_CDTHIS, _rhs);
		}
		// \summation over all entries
		MATRICE_GLOBAL_INL auto sum() const {
			default_type _Ret = 0.0; int n = _CDTHIS->size();
#pragma omp parallel if(n>100)
		{
#pragma omp for reduction(+:_Ret)
			for (int i = 0; i < n; ++i) _Ret += _CDTHIS->operator()(i);
		}
			return (_Ret);
		}

		MATRICE_GLOBAL_FINL size_t size() const { return M*N; }
		MATRICE_GLOBAL_FINL size_t rows() const { return M; }
		MATRICE_GLOBAL_FINL size_t cols() const { return N; }

	protected:
		size_t M, K, N;

#undef _CDTHIS
	};

	// \matrix element-wise binary operation expression: +, -, *, /, ...
	template<typename T, typename U, typename _BinaryOp> class EwiseBinaryExpr 
		: public Base_<EwiseBinaryExpr<T, U, _BinaryOp>>//, _Base_v1<T, U>
	{
	public:
		enum { options = option<ewise>::value };
		enum {_M = 0, _N = 0};
		using value_t = conditional_t<std::is_scalar_v<T>, T, typename T::value_t>;
		static constexpr value_t inf = std::numeric_limits<value_t>::infinity();
		static constexpr auto _Is_T_exp = is_expression_v<T>;

		MATRICE_GLOBAL_INL EwiseBinaryExpr(const value_t _scalar, const U& _rhs) noexcept
			: _Scalar(_scalar), _LHS(_rhs), _RHS(_rhs) { M = _LHS.rows(), N = _RHS.cols(); }
		MATRICE_GLOBAL_INL EwiseBinaryExpr(const T& _lhs, const value_t _scalar) noexcept
			: _Scalar(_scalar), _LHS(_lhs), _RHS(_lhs) { M = _LHS.rows(), N = _RHS.cols(); }
		MATRICE_GLOBAL_INL EwiseBinaryExpr(const T& _lhs, const U& _rhs) noexcept
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
		MATRICE_GLOBAL_INL auto eval_impl() const {
			types::Matrix_<value_t, _M, _N> _Ret(M, N);
			this->assign_to(_Ret);
			return std::forward<decltype(_Ret)>(_Ret);
		}
		template<typename _Mty> MATRICE_GLOBAL_INL void assign_to(_Mty& res) const
		{
#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}

	private:
		const value_t _Scalar = std::numeric_limits<value_t>::infinity();
		const T& _LHS; const U& _RHS;
		_BinaryOp _Op;
		using Base_<EwiseBinaryExpr<T, U, _BinaryOp>>::M;
		using Base_<EwiseBinaryExpr<T, U, _BinaryOp>>::N;
	};

	// \matrix element-wise unary operation expression: abs(), log(), sqrt() ....
	template<typename T, typename _UnaryOp> class EwiseUnaryExpr 
		: public Base_<EwiseUnaryExpr<T, _UnaryOp>>//, _Base_v1<T, T>
	{
		using my_type = EwiseUnaryExpr;
		using const_reference_t = add_const_reference_t<T>;
	public:
		enum { options = expression_options<_UnaryOp>::value,
			_M = conditional_size_v<is_matrix_v<T>, matrix_traits<T>::_M, 0>, 
			_N = conditional_size_v<is_matrix_v<T>, matrix_traits<T>::_N, 0> };
		using value_t = conditional_t<std::is_scalar_v<T>, T, typename T::value_t>;
		EwiseUnaryExpr(const_reference_t _rhs) noexcept
			:_RHS(_rhs) { M = _rhs.rows(), N = _rhs.cols(); }

		MATRICE_GLOBAL_INL value_t operator() (size_t _idx) {
			return _Op(_RHS(_idx));
		}
		MATRICE_GLOBAL_INL const value_t operator() (size_t _idx) const {
			return _Op(_RHS(_idx));
		}
		MATRICE_GLOBAL_INL auto eval_impl() const {
			types::Matrix_<value_t, _M, _N> _Ret(M, N);
			this->assign_to(_Ret);
			return std::forward<decltype(_Ret)>(_Ret);
		}
		template<typename _Mty> MATRICE_GLOBAL_INL void assign_to(_Mty& res) const {
#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}

	private:
		_UnaryOp _Op;
		const_reference_t _RHS;
		using Base_<my_type>::M;
		using Base_<my_type>::N;
	};

	// \matrix binary operation expression: matmul(), ...
	template<class T, class U, typename _BinaryOp> class MatBinaryExpr 
		: public Base_<MatBinaryExpr<T, U, _BinaryOp>>//, _Base_v1<T, U>
	{
	public:
		enum{options = option<_BinaryOp::flag>::value};
		enum {_M = 0, _N = 0};
		using value_t = conditional_t<std::is_class<T>::value, typename T::value_t, default_type>;

		MATRICE_GLOBAL_INL MatBinaryExpr(const T& _lhs, const U& _rhs) noexcept 
			: _LHS(_lhs), _RHS(_rhs) {
			M = _LHS.rows(), K = _LHS.cols(), N = _RHS.cols();
			if (K != _RHS.rows() && K == N) N = _RHS.rows();
		}

		MATRICE_HOST_INL value_t* operator() (value_t* res) const {
			return _Op(_LHS.data(), _RHS.data(), res, M, K, N);
		}
		MATRICE_GLOBAL_INL value_t operator() (int r, int c) const {
			return _Op(_LHS, _RHS, r, c);
		}
		MATRICE_GLOBAL_INL value_t operator() (int _idx) const {
			int r = _idx / N, c = _idx - r * N;
			return _Op(_LHS, _RHS, r, c);
		}
		MATRICE_GLOBAL_INL auto eval_impl() const {
			types::Matrix_<value_t, _M, _N> _Ret(M, N);
			this->assign_to(_Ret);
			return std::forward<decltype(_Ret)>(_Ret);
		}

		template<typename _Mty> MATRICE_GLOBAL_INL void assign_to(_Mty& res) const
		{
#pragma omp parallel for
			for (int i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}
	private:
		const T& _LHS; 
		const U& _RHS;
		_BinaryOp _Op;
		using Base_<MatBinaryExpr<T, U, _BinaryOp>>::M;
		using Base_<MatBinaryExpr<T, U, _BinaryOp>>::K;
		using Base_<MatBinaryExpr<T, U, _BinaryOp>>::N;
	};

	// \matrix unary operation expression: inverse(), transpose(), ...
	template<class T, typename _UnaryOp> class MatUnaryExpr 
		: public Base_<MatUnaryExpr<T, _UnaryOp>>//, _Base_v1<T, T>
	{
	public:
		enum {options = option<_UnaryOp::flag>::value};
		enum {_M = 0, _N = 0};
		using value_t = conditional_t<std::is_class_v<T>, typename T::value_t, conditional_t<std::is_scalar_v<T>, T, default_type>>;
		MATRICE_GLOBAL_FINL MatUnaryExpr(const T& inout) noexcept
			: _RHS(inout), _ANS(inout) {
			M = _RHS.rows(), N = _RHS.cols();
			if constexpr (options & trp == trp) std::swap(M, N);
		}
		MATRICE_GLOBAL_FINL MatUnaryExpr(const T& _rhs, T& _ans)
			: _RHS(_rhs), _ANS(_ans) {
			M = _RHS.rows(), N = _RHS.cols();
			if constexpr (options & trp == trp) std::swap(M, N);
		}
		MATRICE_GLOBAL_FINL MatUnaryExpr(MatUnaryExpr&& _other)
			: options(_other.options) {
			*this = std::move(_other);
		}

		MATRICE_GLOBAL_FINL auto operator()() const { 
			return _Op(M, _ANS.data(), _RHS.data());
		}
		MATRICE_GLOBAL_FINL value_t operator() (int r, int c) const { 
			if constexpr (options & inv == inv) return _ANS(c + r * N);
			if constexpr (options & trp == trp) return _ANS(r + c * N);
		}
		MATRICE_GLOBAL_FINL value_t operator() (int _idx) const { 
			if constexpr (options & trp == trp)
				_idx = _idx * M + _idx / N * (1 - N * M);
			return _ANS(_idx);
		}
		MATRICE_GLOBAL_FINL T& ans() { return _ANS; }
		MATRICE_GLOBAL_FINL const T& ans() const { return _ANS; }
		template<typename _Rhs, typename _Ret = MatBinaryExpr<MatUnaryExpr, _Rhs, Op::MatMul<value_t>>> MATRICE_GLOBAL_INL _Ret mul(const _Rhs& _rhs) {
			if constexpr (options & inv == inv) this->operator()();
			if constexpr (options & trp == trp)                   ;
			return _Ret(*this, _rhs);
		}

		MATRICE_GLOBAL_INL auto eval_impl() const {
			types::Matrix_<value_t, options & trp == trp ? _N : _M, options & trp == trp ? _M : _N> _Ret(M, N);
			this->assign_to(_Ret);
			return std::forward<decltype(_Ret)>(_Ret);
		}

		template<typename _Mty>MATRICE_GLOBAL_INL void assign_to(_Mty& res)const {
			if constexpr (options & inv == inv) {
				_Op(M, res.data(), _RHS.data());
			}
			if constexpr (options & trp == trp) {
#pragma omp parallel for
				for (int i = 0; i < size(); ++i) res(i) = (*this)(i);
			}
		} /*i*N+(1-N*M)*(i/M) is replaced by i at left hand*/
	private:
		const T& _RHS;
		const T& _ANS;
		_UnaryOp _Op;
		using Base_<MatUnaryExpr<T, _UnaryOp>>::M;
		using Base_<MatUnaryExpr<T, _UnaryOp>>::N;
		using Base_<MatUnaryExpr<T, _UnaryOp>>::size;
	};
};
template<typename _Derived>
struct is_expression<Expr::Base_<_Derived>> :std::true_type {};
template<typename T, typename U, typename _Op>
struct is_expression<Expr::EwiseBinaryExpr<T, U, _Op>> :std::true_type {};
template<typename T, typename _Op> 
struct is_expression<Expr::EwiseUnaryExpr<T, _Op>> :std::true_type {};
template<typename T, typename U, typename _Op>
struct is_expression<Expr::MatBinaryExpr<T, U, _Op>> :std::true_type {};
template<typename T, typename _Op>
struct is_expression<Expr::MatUnaryExpr<T, _Op>> :std::true_type {};
template<typename T, typename U, typename _Op>
struct expression_traits<Expr::EwiseBinaryExpr<T, U, _Op>> {
	using abstract_type = conditional_t<is_matrix_v<T>, T, U>;
	using matrix_type = conditional_t<is_matrix_v<abstract_type>, abstract_type, conditional_t<is_expression_v<abstract_type>, typename abstract_type::matrix_type, types::Matrix<typename abstract_type::value_t>>>;
	enum {
		_M = matrix_traits<matrix_type>::size::rows::value,
		_N = matrix_traits<matrix_type>::size::cols::value
	};
};
template<typename T, typename _Op>
struct expression_traits<Expr::EwiseUnaryExpr<T, _Op>> {
	using matrix_type = conditional_t<is_matrix_v<T>, T, conditional_t<is_expression_v<T>, typename T::matrix_type, types::Matrix<typename T::value_t>>>;
	enum {
		_M = matrix_traits<matrix_type>::size::rows::value,
		_N = matrix_traits<matrix_type>::size::cols::value
	};
};
template<typename T, typename U, typename _Op>
struct expression_traits<Expr::MatBinaryExpr<T, U, _Op>> {
	using matrix_type_1 = conditional_t<is_matrix_v<T>, T, conditional_t<is_expression_v<T>, typename T::matrix_type, types::Matrix<typename T::value_t>>>;
	using matrix_type_2 = conditional_t<is_matrix_v<U>, U, conditional_t<is_expression_v<U>, typename U::matrix_type, types::Matrix<typename U::value_t>>>;
	enum {
		_M = matrix_traits<matrix_type_1>::size::rows::value,
		_N = matrix_traits<matrix_type_2>::size::cols::value
	};
	using matrix_type = types::Matrix_<typename matrix_type_1::value_t, _M, _N>;
};
template<typename T, typename _Op>
struct expression_traits<Expr::MatUnaryExpr<T, _Op>> {
	using matrix_type = conditional_t<is_matrix_v<T>, T, conditional_t<is_expression_v<T>, typename T::matrix_type, types::Matrix<typename T::value_t>>>;
	enum {
		_M = matrix_traits<matrix_type>::size::rows::value,
		_N = matrix_traits<matrix_type>::size::cols::value
	};
};
template<
	typename _Exp, 
	typename = std::enable_if_t<is_expression_v<_Exp>>
	>
using _Matrix_exp = typename Expr::Base_<_Exp>;

MATRICE_NAMESPACE_EXPR_END