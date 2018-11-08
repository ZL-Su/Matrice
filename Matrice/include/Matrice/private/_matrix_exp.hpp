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
#include "_size_traits.h"

#pragma warning(disable: 4715)

DGE_MATRICE_BEGIN
#pragma region <!-- Forward declarations and matrix traits supplements -->
// \forward declarations 
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
#pragma endregion
DGE_MATRICE_END

MATRICE_NAMESPACE_EXPR_BEGIN
template<typename _T, typename _U, typename _Op> class EwiseBinaryExpr;
template<typename _T, typename _U, typename _Op> class MatBinaryExpr;
template<typename _T, typename _Op> class MatUnaryExpr;

struct Expr {
	using size_t = std::size_t;
	using default_type = double;
	enum OpFlag { ewise = 0, mmul = 1, inv = 2, trp = 3, undef = -1 };
	template<int _Option> struct option { enum { value = _Option | expr }; };

	/*factorial_t<N> = N!*/
	template<int N> struct factorial_t { enum { value = factorial_t<N - 1>::value*N }; };
	template<> struct factorial_t<0> { enum { value = 1 }; };

#define _MATRICE_DEFEXP_EWISEUOP(_Name) \
template<typename _Ty> struct _Ewise_##_Name { \
	enum {flag = ewise}; \
	MATRICE_GLOBAL_FINL constexpr auto operator() (const _Ty& _Val) const {\
		return (std::_Name(_Val)); \
	}\
};
#define _MATRICE_DEFEXP_EWISEBOP(_Name) \
template<typename _Ty> struct _Ewise_##_Name { \
	enum {flag = ewise}; \
	template<typename _Uy = _Ty> \
	MATRICE_GLOBAL_FINL constexpr auto operator() (const _Ty& _Left, const _Uy& _Right) const {\
		return (std::_Name(_Left, _Right)); \
	}\
};
#define _MATRICE_DEFEXP_ARITHOP(OP, NAME) \
template<typename _Rhs, typename  = std::enable_if_t<std::true_type::value>> \
MATRICE_GLOBAL_FINL auto operator##OP(const _Rhs& _Right) { \
	return EwiseBinaryExpr<derived_t, _Rhs, Op::_Ewise_##NAME<value_t>>(*_CDTHIS, _Right); \
} \
template<typename _Lhs, typename = std::enable_if_t<std::true_type::value>> friend \
	MATRICE_GLOBAL_FINL auto operator##OP(const _Lhs& _Left, const_derived& _Right) { \
	return EwiseBinaryExpr<_Lhs, derived_t, Op::_Ewise_##NAME<value_t>>(_Left, _Right); \
}
	struct Op {
		template<typename _Ty> using _Ewise_add = std::plus<_Ty>;
		template<typename _Ty> using _Ewise_sub = std::minus<_Ty>;
		template<typename _Ty> using _Ewise_mul = std::multiplies<_Ty>;
		template<typename _Ty> using _Ewise_div = std::divides<_Ty>;

		_MATRICE_DEFEXP_EWISEBOP(max)
		_MATRICE_DEFEXP_EWISEBOP(min)
		_MATRICE_DEFEXP_EWISEUOP(sqrt)
		_MATRICE_DEFEXP_EWISEUOP(exp)
		_MATRICE_DEFEXP_EWISEUOP(abs)
		_MATRICE_DEFEXP_EWISEUOP(log)
		_MATRICE_DEFEXP_EWISEUOP(log2)
		_MATRICE_DEFEXP_EWISEUOP(log10)

		template<typename _Ty> struct _Mat_inv { 
			enum { flag = inv }; 
			MATRICE_GLOBAL _Ty* operator()(int M, _Ty* Out, _Ty* In = nullptr) const; 
		};
		template<typename _Ty> struct _Mat_trp { 
			enum { flag = trp }; 
			MATRICE_HOST_FINL _Ty operator()(int M, _Ty* Out, _Ty* i) const { return (Out[int(i[1])*M + int(i[0])]); } 
		};
		// *\matrix spread multiply
		template<typename _Ty> struct _Mat_sprmul {
			enum { flag = undef };
			template<typename _Lhs, typename _Rhs> MATRICE_GLOBAL_FINL
				_Ty operator() (const _Lhs& lhs, const _Rhs& rhs, int r, int c) const { return (lhs(r) * rhs(c)); }
		};
		template<typename _Ty> struct _Mat_mul {
			enum { flag = mmul };
			template<typename _Rhs> MATRICE_GLOBAL_FINL _Ty operator() (const _Ty* lhs, const _Rhs& rhs, int c, int _plh = 0) const
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
			template<typename _Lhs, typename _Rhs> MATRICE_GLOBAL_FINL _Ty operator() (const _Lhs& lhs, const _Rhs& rhs, int r, int c) const
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
#define _CDTHIS static_cast<const_derived_pointer>(this)
		using myt_traits = expression_traits<_Op>;
		using value_type = typename myt_traits::value_type;
		using value_t = value_type;
		using matrix_type = typename myt_traits::auto_matrix_type;
		using derived_t = typename myt_traits::type;
		using derived_pointer = std::add_pointer_t<derived_t>;
		using const_derived = std::add_const_t<derived_t>;
		using const_derived_pointer = std::add_pointer_t<const_derived>;
		enum {options = myt_traits::options,
			CompileTimeRows = myt_traits::rows, 
			CompileTimeCols = myt_traits::cols};

		// \evaluate the matrix expression
		MATRICE_GLOBAL_INL auto eval() const { 
			matrix_type _Ret(M, N);
			_CDTHIS->assign_to(_Ret);
			return std::forward<matrix_type>(_Ret);
		}
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
			return MatBinaryExpr<derived_t, _Rhs, Op::_Mat_mul<value_t>>(*_CDTHIS, _rhs);
		}
		// \summation over all entries
		MATRICE_GLOBAL_INL auto sum() const {
			auto _Ret = value_t(0); int _Size = _CDTHIS->size();
#pragma omp parallel if(_Size > 100)
		{
#pragma omp for reduction(+:_Ret)
			for (int i = 0; i < _Size; ++i) _Ret += _CDTHIS->operator()(i);
		}
			return (_Ret);
		}

		MATRICE_GLOBAL_INL auto operator() (std::size_t x, std::size_t y) {
			return _CDTHIS->operator()(x + y * N);
		}
		MATRICE_GLOBAL_INL const auto operator() (std::size_t x, std::size_t y) const {
			return _CDTHIS->operator()(x + y * N);
		}

		MATRICE_GLOBAL_FINL constexpr std::size_t size() const { return M*N; }
		MATRICE_GLOBAL_FINL constexpr std::size_t rows() const { return M; }
		MATRICE_GLOBAL_FINL constexpr std::size_t cols() const { return N; }
		MATRICE_GLOBAL_FINL constexpr auto shape() const { return std::tie(M, N); }

		_MATRICE_DEFEXP_ARITHOP(+, add)
		_MATRICE_DEFEXP_ARITHOP(-, sub)
		_MATRICE_DEFEXP_ARITHOP(*, mul)
		_MATRICE_DEFEXP_ARITHOP(/, div)

	protected:
		std::size_t M, K, N;

#undef _CDTHIS
	};

	// \matrix element-wise binary operation expression: +, -, *, /, ...
	template<typename T, typename U, typename _BinaryOp, bool _T = std::is_scalar_v<T>, bool _U = std::is_scalar_v<U>> class EwiseBinaryExpr {};
	template<typename T, typename U, typename _BinaryOp>
	class EwiseBinaryExpr<T, U, _BinaryOp, false, false> : public Base_<EwiseBinaryExpr<T, U, _BinaryOp, false, false>>
	{
		using _Base_type = Base_<EwiseBinaryExpr<T, U, _BinaryOp, false, false>>;
	public:
		using _Base_type::CompileTimeRows;
		using _Base_type::CompileTimeCols;
		using typename _Base_type::value_t;
		using typename _Base_type::matrix_type;
		using _Base_type::operator();
		enum { options = option<ewise>::value };

		MATRICE_GLOBAL_INL EwiseBinaryExpr(const T& _lhs, const U& _rhs) noexcept
			: _LHS(_lhs), _RHS(_rhs) { M = _LHS.rows(), N = _RHS.cols(); }

		MATRICE_GLOBAL_INL value_t operator() (std::size_t _idx) { 
			return _Op(_LHS(_idx), _RHS(_idx));
		}
		MATRICE_GLOBAL_INL const value_t operator() (std::size_t _idx) const {
			return _Op(_LHS(_idx), _RHS(_idx));
		}

		template<typename _Mty> 
		MATRICE_GLOBAL_INL void assign_to(_Mty& res) const {
#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}

	private:
		const T& _LHS; const U& _RHS;
		_BinaryOp _Op;
		using _Base_type::M;
		using _Base_type::N;
	};
	template<typename T, typename U, typename _BinaryOp>
	class EwiseBinaryExpr<T, U, _BinaryOp, true, false> : public Base_<EwiseBinaryExpr<T, U, _BinaryOp, true, false>>
	{
		using _Base_type = Base_<EwiseBinaryExpr<T, U, _BinaryOp, true, false>>;
	public:
		using _Base_type::CompileTimeRows;
		using _Base_type::CompileTimeCols;
		using typename _Base_type::value_t;
		using typename _Base_type::matrix_type;
		using _Base_type::operator();
		enum { options = option<ewise>::value };

		MATRICE_GLOBAL_INL EwiseBinaryExpr(const T _scalar, const U& _rhs) noexcept
			: _Scalar(_scalar), _RHS(_rhs) {
			M = _RHS.rows(), N = _RHS.cols();
		}

		MATRICE_GLOBAL_INL value_t operator() (size_t _idx) {
			return _Op(_Scalar, _RHS(_idx));
		}
		MATRICE_GLOBAL_INL const value_t operator() (size_t _idx) const {
			return _Op(_Scalar, _RHS(_idx));
		}

		template<typename _Mty> 
		MATRICE_GLOBAL_INL void assign_to(_Mty& res) const {
#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}

	private:
		const value_t _Scalar = std::numeric_limits<value_t>::infinity();
		const U& _RHS;
		_BinaryOp _Op;
		using _Base_type::M;
		using _Base_type::N;
	};
	template<typename T, typename U, typename _BinaryOp>
	class EwiseBinaryExpr<T, U, _BinaryOp, false, true> : public Base_<EwiseBinaryExpr<T, U, _BinaryOp, false, true>>
	{
		using _Base_type = Base_<EwiseBinaryExpr<T, U, _BinaryOp, false, true>>;
	public:
		using _Base_type::CompileTimeRows;
		using _Base_type::CompileTimeCols;
		using typename _Base_type::value_t;
		using typename _Base_type::matrix_type;
		using _Base_type::operator();
		enum { options = option<ewise>::value };

		MATRICE_GLOBAL_INL EwiseBinaryExpr(const T& _lhs, const U _scalar) noexcept
			: _Scalar(_scalar), _LHS(_lhs) {
			M = _LHS.rows(), N = _LHS.cols();
		}

		MATRICE_GLOBAL_INL value_t operator() (size_t _idx) {
			return _Op(_LHS(_idx), _Scalar);
		}
		MATRICE_GLOBAL_INL const value_t operator() (size_t _idx) const {
			return _Op(_LHS(_idx), _Scalar);
		}

		template<typename _Mty> 
		MATRICE_GLOBAL_INL void assign_to(_Mty& res) const {
#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}

	private:
		const value_t _Scalar = std::numeric_limits<value_t>::infinity();
		const T& _LHS;
		_BinaryOp _Op;
		using _Base_type::M;
		using _Base_type::N;
	};

	// \matrix element-wise unary operation expression: abs(), log(), sqrt() ....
	template<typename T, typename _UnaryOp> 
	class EwiseUnaryExpr : public Base_<EwiseUnaryExpr<T, _UnaryOp>>
	{
		using const_reference_t = add_const_reference_t<T>;
		using _Base_type = Base_<EwiseUnaryExpr<T, _UnaryOp>>;
	public:
		using _Base_type::CompileTimeRows;
		using _Base_type::CompileTimeCols;
		using typename _Base_type::value_t;
		using typename _Base_type::matrix_type;
		using _Base_type::operator();
		enum { options = expression_options<_UnaryOp>::value };

		EwiseUnaryExpr(const_reference_t _rhs) noexcept
			:_RHS(_rhs) { M = _rhs.rows(), N = _rhs.cols(); }

		MATRICE_GLOBAL_INL value_t operator() (size_t _idx) {
			return _Op(_RHS(_idx));
		}
		MATRICE_GLOBAL_INL const value_t operator() (size_t _idx) const {
			return _Op(_RHS(_idx));
		}

		template<typename _Mty> MATRICE_GLOBAL_INL void assign_to(_Mty& res) const {
#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}

	private:
		_UnaryOp _Op;
		const_reference_t _RHS;
		using _Base_type::M;
		using _Base_type::N;
	};

	// \matrix binary operation expression: matmul(), ...
	template<class T, class U, typename _BinaryOp> 
	class MatBinaryExpr : public Base_<MatBinaryExpr<T, U, _BinaryOp>>
	{
		using _Base_type = Base_<MatBinaryExpr<T, U, _BinaryOp>>;
	public:
		using _Base_type::CompileTimeRows;
		using _Base_type::CompileTimeCols;
		using typename _Base_type::value_t;
		using typename _Base_type::matrix_type;
		enum{options = option<_BinaryOp::flag>::value};

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

		template<typename _Mty> MATRICE_GLOBAL_INL void assign_to(_Mty& res) const
		{
#pragma omp parallel for
			for (int i = 0; i < res.size(); ++i) res(i) = this->operator()(i);
		}
	private:
		const T& _LHS; 
		const U& _RHS;
		_BinaryOp _Op;
		using _Base_type::M;
		using _Base_type::K;
		using _Base_type::N;
	};

	// \matrix unary operation expression: inverse(), transpose(), ...
	template<class T, typename _UnaryOp> 
	class MatUnaryExpr : public Base_<MatUnaryExpr<T, _UnaryOp>>
	{
		using _Base_type = Base_<MatUnaryExpr<T, _UnaryOp>>;
	public:
		using _Base_type::CompileTimeRows;
		using _Base_type::CompileTimeCols;
		using typename _Base_type::value_t;
		using typename _Base_type::matrix_type;
		enum {options = option<_UnaryOp::flag>::value};

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
		template<typename _Rhs, typename _Ret = MatBinaryExpr<MatUnaryExpr, _Rhs, Op::_Mat_mul<value_t>>> 
		MATRICE_GLOBAL_INL _Ret mul(const _Rhs& _rhs) {
			if constexpr (options & inv == inv) this->operator()();
			if constexpr (options & trp == trp)                   ;
			return _Ret(*this, _rhs);
		}

		template<typename _Mty>MATRICE_GLOBAL_INL void assign_to(_Mty& res)const {
			if constexpr (options & inv == inv) {
				_Op(M, res.data(), _RHS.data());
			}
			if constexpr (options & trp == trp) {
				const int _Size = this->size();
				for (int i = 0; i < _Size; ++i) res(i) = (*this)(i);
			}
		} /*i*N+(1-N*M)*(i/M) is replaced by i at left hand*/
	private:
		const T& _RHS;
		const T& _ANS;
		_UnaryOp _Op;
		using _Base_type::M;
		using _Base_type::N;
		using _Base_type::size;
	};
#undef _MATRICE_DEFEXP_ARITHOP
#undef _MATRICE_DEFEXP_EWISEUOP
#undef _MATRICE_DEFEXP_EWISEBOP
};

template<typename _Exp, typename = std::enable_if_t<is_expression_v<_Exp>>>
using _Matrix_exp = typename Expr::Base_<_Exp>;

MATRICE_NAMESPACE_EXPR_END

DGE_MATRICE_BEGIN

using _Exp_op = exprs::Expr::Op;

template<typename _Derived>
struct is_expression<exprs::Expr::Base_<_Derived>> :std::true_type {};
template<typename T, typename U, typename _Op, bool _T, bool _U>
struct is_expression<exprs::Expr::EwiseBinaryExpr<T, U, _Op, _T, _U>> :std::true_type {};
template<typename T, typename _Op>
struct is_expression<exprs::Expr::EwiseUnaryExpr<T, _Op>> :std::true_type {};
template<typename T, typename U, typename _Op>
struct is_expression<exprs::Expr::MatBinaryExpr<T, U, _Op>> :std::true_type {};
template<typename T, typename _Op>
struct is_expression<exprs::Expr::MatUnaryExpr<T, _Op>> :std::true_type {};

template<typename T, typename U, typename _Op>
struct expression_traits<exprs::Expr::EwiseBinaryExpr<T, U, _Op, false, false>> {
	using type = exprs::Expr::EwiseBinaryExpr<T, U, _Op, false, false>;
	using value_type = std::common_type_t<typename T::value_t, typename U::value_t>;
	enum {
		options = exprs::Expr::option<exprs::Expr::OpFlag::ewise>::value,
		rows = T::CompileTimeRows,
		cols = U::CompileTimeCols,
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};
template<typename T, typename U, typename _Op>
struct expression_traits<exprs::Expr::EwiseBinaryExpr<T, U, _Op, true, false>> {
	using type = exprs::Expr::EwiseBinaryExpr<T, U, _Op, true, false>;
	using value_type = typename U::value_t;
	enum {
		options = exprs::Expr::option<exprs::Expr::OpFlag::ewise>::value,
		rows = U::CompileTimeRows,
		cols = U::CompileTimeCols,
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};
template<typename T, typename U, typename _Op>
struct expression_traits<exprs::Expr::EwiseBinaryExpr<T, U, _Op, false, true>> {
	using type = exprs::Expr::EwiseBinaryExpr<T, U, _Op, false, true>;
	using value_type = typename T::value_t;
	enum {
		options = exprs::Expr::option<exprs::Expr::OpFlag::ewise>::value,
		rows = T::CompileTimeRows,
		cols = T::CompileTimeCols,
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};

template<typename T, typename _Op>
struct expression_traits<exprs::Expr::EwiseUnaryExpr<T, _Op>> {
	using type = exprs::Expr::EwiseUnaryExpr<T, _Op>;
	using value_type = typename T::value_t;
	enum {
		options = expression_options<_Op>::value,
		rows = T::CompileTimeRows,
		cols = T::CompileTimeCols
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};
template<typename T, typename U, typename _Op>
struct expression_traits<exprs::Expr::MatBinaryExpr<T, U, _Op>> {
	using type = exprs::Expr::MatBinaryExpr<T, U, _Op>;
	using value_type = std::common_type_t<typename T::value_t, typename U::value_t>;
	enum {
		options = expression_options<_Op>::value,
		rows = conditional_size_v<T::CompileTimeRows <= 0 || U::CompileTimeRows <= 0, 0, max_integer_v<T::CompileTimeRows, 0>>,
		cols = conditional_size_v<T::CompileTimeCols <= 0 || U::CompileTimeCols <= 0, 0, max_integer_v<U::CompileTimeCols, 0>>
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};
template<typename T, typename _Op>
struct expression_traits<exprs::Expr::MatUnaryExpr<T, _Op>> {
	using type = exprs::Expr::MatUnaryExpr<T, _Op>;
	using value_type = typename T::value_t;
	enum {
		options = expression_options<_Op>::value,
		rows = conditional_size_v<options&exprs::Expr::OpFlag::trp == exprs::Expr::OpFlag::trp, T::CompileTimeCols, T::CompileTimeRows>,
		cols = conditional_size_v<options&exprs::Expr::OpFlag::trp == exprs::Expr::OpFlag::trp, T::CompileTimeRows, T::CompileTimeCols>
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};

using exprs::Expr;
// *\element-wise addition
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<std::is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_add<value_t>>>
	MATRICE_GLOBAL_FINL auto operator+ (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise subtraction
template<
	typename _Lhs, class _Rhs,
	typename value_t = conditional_t<std::is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_sub<value_t>>>
	MATRICE_GLOBAL_FINL auto operator- (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise multiplication
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<std::is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_mul<value_t>>>
	MATRICE_GLOBAL_FINL auto operator* (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise division
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<std::is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_div<value_t>>>
	MATRICE_GLOBAL_FINL auto operator/ (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise maximum
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<std::is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_max<value_t>>>
	MATRICE_GLOBAL_FINL auto max(const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise minimum
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<std::is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_min<value_t>>>
	MATRICE_GLOBAL_FINL auto min(const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }


// *\element-wise sqrt()
template<
	typename _Rhs,
	typename value_t = std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_sqrt<value_t>>>
	MATRICE_GLOBAL_FINL auto sqrt(const _Rhs& _right) { return _Op(_right); }

// *\element-wise exp()
template<
	typename _Rhs,
	typename value_t = std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_exp<value_t>>>
	MATRICE_GLOBAL_FINL auto exp(const _Rhs& _right) { return _Op(_right); }

// *\element-wise log()
template<
	typename _Rhs,
	typename value_t = std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_log<value_t>>>
	MATRICE_GLOBAL_FINL auto log(const _Rhs& _right) { return _Op(_right); }

// *\element-wise abs()
template<
	typename _Rhs,
	typename value_t = std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename     _Op = Expr::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_abs<value_t>>>
	MATRICE_GLOBAL_FINL auto abs(const _Rhs& _right) { return _Op(_right); }

// *\transpose expression
template<
	typename _Rhs,
	typename value_t = std::enable_if_t<std::is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename _Op = Expr::MatUnaryExpr<_Rhs, _Exp_op::_Mat_trp<value_t>>>
	MATRICE_GLOBAL_FINL auto transpose(const _Rhs& _right) { return _Op(_right); }

// *\outer product expression : xy^T
template<
	typename _Lhs, typename _Rhs,
	typename value_t = std::common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = Expr::MatBinaryExpr<_Lhs, _Rhs, _Exp_op::_Mat_sprmul<value_t>>>
	MATRICE_GLOBAL_FINL auto outer_product(const _Lhs& _left, const _Rhs& _right) {
	return _Op(_left, _right);
}

DGE_MATRICE_END