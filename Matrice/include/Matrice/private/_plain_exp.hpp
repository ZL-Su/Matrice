/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include "_plain_shape.hpp"
#include "_type_traits.h"
#include "_size_traits.h"
#include "_tag_defs.h"
#include "../private/math/_primitive_funcs.hpp"
#if defined(MATRICE_SIMD_ARCH)
#include "../arch/ixpacket.h"
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)
#endif

DGE_MATRICE_BEGIN
#pragma region <!-- Forward declarations and matrix traits supplements -->
// \forward declarations 
_TYPES_BEGIN
template<typename _Ty> class Matrix;
template<typename _Ty, int _Rows, int _cols> class Matrix_;
template<typename _Derived, typename _Traits, typename _Ty> class Base_;
_TYPES_END
_DETAIL_BEGIN
//template<typename _Ty> class _Tensor;
template<typename _Ty, size_t _Depth> class _Tensor;
_DETAIL_END

template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
constexpr static _Ty zero = static_cast<_Ty>(0);
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
constexpr static _Ty one = static_cast<_Ty>(1);
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
constexpr static _Ty two = static_cast<_Ty>(2);

template<typename _Ty, int _M, int _N> 
struct is_matrix<types::Matrix_<_Ty, _M, _N>> : std::true_type {};

template<typename _Ty, int _M>
struct is_fxdvector<types::Matrix_<_Ty, _M, 1>> : std::true_type {};

template<typename _Derived, typename _Traits, typename _Ty> 
struct is_matrix<types::Base_<_Derived, _Traits, _Ty>> : std::true_type {};

template<typename _Ty, size_t _Depth>
struct is_tensor<detail::_Tensor<_Ty, _Depth>> : std::true_type {};

template<typename _Derived, typename _Traits, typename _Ty>
struct matrix_traits<types::Base_<_Derived, _Traits, _Ty>> {
	using type = _Ty;
	enum { _M = _Traits::_M, _N = _Traits::_N };
	static constexpr bool Is_base = std::true_type::value;
};
template<typename _Ty, int _Rows, int _Cols> 
struct matrix_traits<types::Matrix_<_Ty, _Rows, _Cols>> {
	using type = _Ty;
	using category = tag::_Matrix_tag;
	enum { _M = _Rows, _N = _Cols };
	static constexpr bool Is_base = std::false_type::value;
};
template<typename _Ty>
struct matrix_traits<detail::_Tensor<_Ty, 0>> {
	using type = _Ty;
	using category = tag::_Tensor_tag;
	static constexpr auto _M = 0, _N = 0;
	static constexpr auto Is_base = std::false_type::value;
};

#pragma endregion
DGE_MATRICE_END

MATRICE_NAMESPACE_EXPR_BEGIN
template<typename _T, typename _U, typename _Op> class EwiseBinaryExpr;
template<typename _T, typename _U, typename _Op> class MatBinaryExpr;
template<typename _T, typename _Op> class MatUnaryExpr;

struct Expr {
#define _MATRICE_DEFEXP_EWISEUOP(NAME) \
template<typename _Ty> struct _Ewise_##NAME { \
	enum {flag = ewise}; \
	using category = tag::_Ewise_##NAME##_tag; \
	MATRICE_GLOBAL_FINL constexpr auto operator() (const _Ty& _Val) const {\
		return (NAME(_Val)); \
	}\
};
#define _MATRICE_DEFEXP_EWISEBOP(NAME) \
template<typename _Ty> struct _Ewise_##NAME { \
	enum {flag = ewise}; \
	using category = tag::_Ewise_##NAME##_tag; \
	template<typename _Uy = _Ty> \
	MATRICE_GLOBAL_FINL constexpr auto operator() (const _Ty& _Left, const _Uy& _Right) const {\
		return (NAME(_Left, _Right)); \
	}\
};
#define _MATRICE_DEFEXP_ARITHOP(OP, NAME) \
template<typename _Rhs, MATRICE_ENABLE_IF(std::true_type::value)> \
MATRICE_GLOBAL_FINL auto operator##OP(const _Rhs& _Right) { \
	return EwiseBinaryExpr<derived_t, _Rhs, Op::_Ewise_##NAME<value_t>>(*_CDTHIS, _Right); \
} \
template<typename _Lhs, MATRICE_ENABLE_IF(std::true_type::value)> friend \
MATRICE_GLOBAL_FINL auto operator##OP(const _Lhs& _Left, const_derived& _Right) { \
	return EwiseBinaryExpr<_Lhs, derived_t, Op::_Ewise_##NAME<value_t>>(_Left, _Right); \
} 
	using default_type = double;
	/**
	 * \expression tags
	 */
	enum OpFlag {
		ewise = 0, mmul = 1, inv = 2, trp = 3, sum = 4, undef = -1
	};
	/**
	 * \expression option
	 */
	template<int _Option> struct option {
		enum { value = _Option };
	};
	/**
	 * \factorial_t<N> = N!
	 */
	template<int N> struct factorial_t {
		enum { value = factorial_t<N - 1>::value*N };
	};
	template<> struct factorial_t<0> {
		enum { value = 1 };
	};

	/**
	 * \expression operators
	 */
	struct Op {
		_MATRICE_DEFEXP_EWISEBOP(add);
		_MATRICE_DEFEXP_EWISEBOP(sub);
		_MATRICE_DEFEXP_EWISEBOP(mul);
		_MATRICE_DEFEXP_EWISEBOP(div);
		_MATRICE_DEFEXP_EWISEBOP(max);
		_MATRICE_DEFEXP_EWISEBOP(min);
		_MATRICE_DEFEXP_EWISEUOP(sqrt);
		_MATRICE_DEFEXP_EWISEUOP(exp);
		_MATRICE_DEFEXP_EWISEUOP(abs);
		_MATRICE_DEFEXP_EWISEUOP(log);
		_MATRICE_DEFEXP_EWISEUOP(log2);
		_MATRICE_DEFEXP_EWISEUOP(log10);
		_MATRICE_DEFEXP_EWISEUOP(floor);

		/**
		 * \accumulates all elements of input expression
		 */
		template<typename _Ty> struct _Accum_exp {
			enum { flag = sum };
			using value_type = _Ty;
			MATRICE_GLOBAL_INL auto operator()(const value_type& _Exp) {
				return _Exp.sum();
			}
		};

		/**
		 * \matrix inverse operator
		 */
		template<typename _Ty> struct _Mat_inv { 
			enum { flag = inv };
			using value_type = _Ty;
			using category = tag::_Matrix_inv_tag;
			MATRICE_GLOBAL _Ty* operator()(int M, value_type* Out, value_type* In = nullptr) const;
		};
		
		/**
		 * \matrix transpose operator
		 */
		template<typename _Ty> struct _Mat_trp { 
			enum { flag = trp };
			using value_type = _Ty;
			using category = tag::_Matrix_trp_tag;
			MATRICE_HOST_FINL value_type operator()(int M, value_type* Out, value_type* i) const noexcept { 
				return (Out[int(i[1])*M + int(i[0])]); 
			}
		};
		
		/**
		 * \matrix spread multiply
		 */
		template<typename _Ty> struct _Mat_sprmul {
			enum { flag = undef };
			using value_type = _Ty;
			template<typename _Lhs, typename _Rhs> MATRICE_GLOBAL_FINL
			value_type operator() (const _Lhs& lhs, const _Rhs& rhs, int r, int c) const noexcept { 
				return (lhs(r) * rhs(c)); 
			}
		};

		/**
		 * \matrix multiplication
		 */
		template<typename _Ty> struct _Mat_mul {
			enum { flag = mmul };
			using category = tag::_Matrix_mul_tag;
			using value_type = _Ty;
			using packet_type = simd::Packet_<value_type, packet_size_v>;

			template<typename _Rhs> MATRICE_GLOBAL_FINL
			value_type operator() (const _Ty* lhs, const _Rhs& rhs, int c, int _plh = 0) const noexcept
			{
				value_type val = value_type(0);
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
			value_type operator() (const _Lhs& lhs, const _Rhs& rhs, int r, int c) const noexcept {
				value_type _Ret = value_type(0);
				const int K = rhs.rows(), N = rhs.cols(), _Idx = r * lhs.cols();
#ifdef __disable_simd__
				for (auto k = 0; k < K; ++k) _Ret += lhs(_Idx + k) * rhs(k*N + c);
#else
#ifdef __AVX__
				for (auto k = 0; k < K; ++k) _Ret += lhs(_Idx + k) * rhs(k*N + c);
#endif
#endif
				return (_Ret);
			}
			template<typename _Lhs, typename _Rhs> MATRICE_GLOBAL_FINL
			value_type operator() (const _Lhs& _L, const _Rhs& _R, int r, int c, tag::_Matrix_tag, tag::_Matrix_tag) {
				const int K = _R.rows(), N = _R.cols(), _Idx = r * _L.cols();
				return detail::template _Blas_kernel_impl<value_type>::dot(_L(_Idx), _R(c), K, 1, N);
			}
			template<typename _Lhs, typename _Rhs> MATRICE_GLOBAL_FINL
			value_type operator() (const _Lhs& _L, const _Rhs& _R, int r, int c, tag::_Matrix_tag, tag::_Matrix_view_tag) {
				const int K = _R.rows(), N = _R.cols(), _Idx = r * _L.cols();
				value_type _Ret = value_type(0);
				for (auto k = 0; k < K; ++k) _Ret += _L(_Idx + k) * _R(k*N + c);
				return (_Ret);
			}
		};

	};

	/**
	 * \expression base class
	 */
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

		MATRICE_GLOBAL_FINL Base_() {}
		MATRICE_GLOBAL_FINL Base_(const basic_shape_t& _Shape) 
			: Shape(_Shape), M(_Shape.rows()), N(_Shape.cols()) {}

		// \evaluate the matrix expression
		MATRICE_GLOBAL_INL auto eval() const {
			matrix_type _Ret(M, N);
			_CDTHIS->assign_to(_Ret);
			return std::forward<matrix_type>(_Ret);
		}
		template<typename _Ret> 
		MATRICE_GLOBAL_INL _Ret eval() const {
			_Ret ret(_CDTHIS->rows(), _CDTHIS->cols());
			return std::forward<_Ret>(ret = *_CDTHIS);
		}
		// \retrieve the evaluated result
		template<typename _Ret>
		MATRICE_GLOBAL_INL _Ret& assign(_Ret& res) const noexcept {
#ifdef _DEBUG
			if (res.size() == 0) res.create(M, N);
#endif
			_CDTHIS->assign_to(res);
			return (res);
		}
		// \formulate matmul expression
		template<typename _Rhs>
		MATRICE_GLOBAL_INL auto mul(const _Rhs& _rhs) const noexcept {
			return MatBinaryExpr<derived_t, _Rhs, Op::_Mat_mul<value_t>>(*_CDTHIS, _rhs);
		}

		/**
		 *\brief spread this to multiply with _rhs element-wisely; if _rhs can be evaluated to a square matrix, this should be unrolled to a 1D array if it is not.
		 *\param [_rhs] must be a matrix with a type of Matrix_<>, or an expression can be evaluated to a matrix.
		 */
		template<typename _Rhs>
		MATRICE_GLOBAL_INL auto spreadmul(const _Rhs& _rhs) const {
			conditional_t<is_matrix_v<_Rhs>, 
				_Rhs, typename _Rhs::matrix_type> _Res(_rhs.shape());
			//spread this along each row to mul. with _rhs
			if (N == 1 || size() == _rhs.rows()) {
#pragma omp parallel for if(_Res.rows() > 100)
				{
					for (index_t r = 0; r < _Res.rows(); ++r) {
						const auto _Val = _CDTHIS->operator()(r);
						auto _Ptr = _Res[r];
						for (index_t c = 0; c < _Res.cols(); ++c) {
							_Ptr[c] = _Val * this->operator()(c, r);
						}
					}
				}
			}
			// spread each entry along column to mul. with _rhs
			else if (M == 1 || size() == _Res.cols()) {
#pragma omp parallel for if(_Res.cols() > 100)
				{
					for (index_t r = 0; r < _Res.rows(); ++r) {
						auto _Ptr = _Res[r];
						for (index_t c = 0; c < _Res.cols(); ++c) {
							_Ptr[c] = _CDTHIS->operator()(c)*this->operator()(c,r);
						}
					}
				}
			}
			else {
				DGELOM_ERROR("Only one-dimension array spread is supported.");
			}
			return forward<decltype(_Res)>(_Res);
		}

		/**
		 * \sum over all expression entries
		 */
		MATRICE_GLOBAL_INL auto sum() const noexcept {
			auto _Ret = value_t(0); int _Size = _CDTHIS->size();
#pragma omp parallel if(_Size > 1000)
		{
#pragma omp for reduction(+:_Ret)
			for (int i = 0; i < _Size; ++i) _Ret += _CDTHIS->operator()(i);
		}
			return (_Ret);
		}
		/**
		 * \average of all entries
		 */
		MATRICE_GLOBAL_INL auto avg() const noexcept {
			return (_CDTHIS->sum() / _CDTHIS->size());
		}
		/**
		 * \variation of all entries
		 */
		MATRICE_GLOBAL_INL auto var() const;

		/**
		 * \operator for expression evaluation
		 */
		MATRICE_GLOBAL_INL auto operator()(tag::_Var_tag = tag::_Expression_eval_tag()) const {
			matrix_type _Ret(M, N);
			_CDTHIS->assign_to(_Ret);
			return std::forward<matrix_type>(_Ret);
		}
		/**
		 * \2-d ewise evaluation operator
		 */
		MATRICE_GLOBAL_FINL auto operator()(size_t x, size_t y) noexcept {
			return _CDTHIS->operator()(x + y * N);
		}
		MATRICE_GLOBAL_FINL const auto operator()(size_t x, size_t y)const noexcept {
			return _CDTHIS->operator()(x + y * N);
		}

		MATRICE_GLOBAL_FINL constexpr size_t size() const noexcept { return M*N; }
		MATRICE_GLOBAL_FINL constexpr size_t rows() const noexcept { return M; }
		MATRICE_GLOBAL_FINL constexpr size_t cols() const noexcept { return N; }

		/**
		 *\brief Get full shape of dst. {N, {D, {H, W}}}
		 */
		MATRICE_GLOBAL_FINL const basic_shape_t& shape() const noexcept {
			return (Shape); 
		}
		MATRICE_GLOBAL_FINL basic_shape_t& shape() noexcept {
			return (Shape);
		}

		/**
		 *\brief Get full dims {N,{C,{H,W}}}
		 */
		MATRICE_GLOBAL_FINL constexpr auto& dims() const noexcept {
			return (Shape);
		}

		_MATRICE_DEFEXP_ARITHOP(+, add)
		_MATRICE_DEFEXP_ARITHOP(-, sub)
		_MATRICE_DEFEXP_ARITHOP(*, mul)
		_MATRICE_DEFEXP_ARITHOP(/, div)

	protected:
		size_t M, K, N;
		basic_shape_t Shape;
#undef _CDTHIS
	};

	// \matrix element-wise binary operation expression: +, -, *, /, ...
	template<typename T, typename U, typename _BinaryOp, 
		bool _T = is_scalar_v<T>, bool _U = is_scalar_v<U>> 
		class EwiseBinaryExpr {};

	template<typename T, typename U, typename _BinaryOp>
	class EwiseBinaryExpr<T, U, _BinaryOp, false, false> 
		: public Base_<EwiseBinaryExpr<T, U, _BinaryOp, false, false>>
	{
		using _Mybase = Base_<EwiseBinaryExpr<T,U,_BinaryOp,false,false>>;
	public:
		using _Mybase::CompileTimeRows;
		using _Mybase::CompileTimeCols;
		using typename _Mybase::value_t;
		using typename _Mybase::matrix_type;
		using _Mybase::operator();
		using category = category_type_t<_BinaryOp>;
		enum { options = option<ewise>::value };

		MATRICE_GLOBAL_INL EwiseBinaryExpr(const T& _lhs, const U& _rhs)
		 :_Mybase(union_shape(_lhs.shape(),_rhs.shape())),
			_LHS(_lhs), _RHS(_rhs) {}

		/*MATRICE_GLOBAL_FINL value_t operator()(size_t _idx) noexcept { 
			return _Op(_LHS(_idx), _RHS(_idx));
		}*/
		MATRICE_GLOBAL_FINL value_t operator()(size_t _idx) const noexcept{
			return _Op(_LHS(_idx), _RHS(_idx));
		}

		template<typename _Mty> 
		MATRICE_GLOBAL_INL void assign_to(_Mty& _Res) const noexcept{
			if (_LHS.size() == _RHS.size()) { //element-wise operation
//#pragma omp parallel for if(_Res.size() > 100)
				for (index_t i = 0; i < _Res.size(); ++i)
					_Res(i) = this->operator()(i);
			}
			else { //spreaded element-wise operation
//#pragma omp parallel for if(_Res.size() > 100)
				for (index_t i = 0; i < _Res.size(); ++i)
					_Res(i) = _Op(_LHS(i/N), _RHS(i));
			}
		}

	private:
		const T& _LHS; const U& _RHS;
		_BinaryOp _Op;
		using _Mybase::M;
		using _Mybase::N;
	};
	template<typename T, typename U, typename _BinaryOp>
	class EwiseBinaryExpr<T, U, _BinaryOp, true, false> 
		: public Base_<EwiseBinaryExpr<T, U, _BinaryOp, true, false>>
	{
		using _Mybase = Base_<EwiseBinaryExpr<T, U, _BinaryOp, true, false>>;
	public:
		using _Mybase::CompileTimeRows;
		using _Mybase::CompileTimeCols;
		using typename _Mybase::value_t;
		using typename _Mybase::matrix_type;
		using _Mybase::operator();
		using category = category_type_t<_BinaryOp>;
		enum { options = option<ewise>::value };

		MATRICE_GLOBAL_INL EwiseBinaryExpr(const T _scalar, const U& _rhs)
			noexcept :_Mybase(_rhs.shape()), _Scalar(_scalar), _RHS(_rhs) {}

		/*MATRICE_GLOBAL_FINL value_t operator() (size_t _idx) {
			return _Op(_Scalar, _RHS(_idx));
		}*/
		MATRICE_GLOBAL_FINL value_t operator()(size_t _idx) const noexcept {
			return _Op(_Scalar, _RHS(_idx));
		}

		template<typename _Mty> 
		MATRICE_GLOBAL_INL void assign_to(_Mty& res) const noexcept {
//#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) 
				res(i) = this->operator()(i);
		}

	private:
		const value_t _Scalar = std::numeric_limits<value_t>::infinity();
		const U& _RHS;
		_BinaryOp _Op;
		using _Mybase::M;
		using _Mybase::N;
	};
	template<typename T, typename U, typename _BinaryOp>
	class EwiseBinaryExpr<T, U, _BinaryOp, false, true> 
		: public Base_<EwiseBinaryExpr<T, U, _BinaryOp, false, true>>
	{
		using _Mybase = Base_<EwiseBinaryExpr<T, U, _BinaryOp, false, true>>;
	public:
		using _Mybase::CompileTimeRows;
		using _Mybase::CompileTimeCols;
		using typename _Mybase::value_t;
		using typename _Mybase::matrix_type;
		using _Mybase::operator();
		using category = category_type_t<_BinaryOp>;
		enum { options = option<ewise>::value };

		MATRICE_GLOBAL_INL EwiseBinaryExpr(const T& _lhs, const U _scalar)
			noexcept :_Mybase(_lhs.shape()), _Scalar(_scalar), _LHS(_lhs) {}

		/*MATRICE_GLOBAL_FINL value_t operator() (size_t _idx) {
			return _Op(_LHS(_idx), _Scalar);
		}*/
		MATRICE_GLOBAL_FINL value_t operator() (size_t _idx) const noexcept {
			return _Op(_LHS(_idx), _Scalar);
		}

		template<typename _Mty> 
		MATRICE_GLOBAL_INL void assign_to(_Mty& res) const noexcept {
//#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) 
				res(i) = this->operator()(i);
		}

	private:
		const value_t _Scalar = std::numeric_limits<value_t>::infinity();
		const T& _LHS;
		_BinaryOp _Op;
		using _Mybase::M;
		using _Mybase::N;
	};

	// \matrix element-wise unary operation expression: abs(), log(), sqrt() ....
	template<typename T, typename _UnaryOp> 
	class EwiseUnaryExpr : public Base_<EwiseUnaryExpr<T, _UnaryOp>>
	{
		using const_reference_t = add_const_reference_t<T>;
		using _Mybase = Base_<EwiseUnaryExpr<T, _UnaryOp>>;
	public:
		using _Mybase::CompileTimeRows;
		using _Mybase::CompileTimeCols;
		using typename _Mybase::value_t;
		using typename _Mybase::matrix_type;
		using _Mybase::operator();
		using category = category_type_t<_UnaryOp>;
		enum { options = expression_options<_UnaryOp>::value };

		EwiseUnaryExpr(const_reference_t _rhs) noexcept
			:_Mybase(_rhs.shape()), _RHS(_rhs) {}

		/*MATRICE_GLOBAL_FINL value_t operator() (size_t _idx) {
			return _Op(_RHS(_idx));
		}*/
		MATRICE_GLOBAL_FINL const value_t operator() (size_t _idx) const noexcept {
			return _Op(_RHS(_idx));
		}

		template<typename _Mty> 
		MATRICE_GLOBAL_INL void assign_to(_Mty& res) const noexcept {
//#pragma omp parallel for
			for (index_t i = 0; i < res.size(); ++i) 
				res(i) = this->operator()(i);
		}

	private:
		_UnaryOp _Op;
		const_reference_t _RHS;
		using _Mybase::M;
		using _Mybase::N;
	};

	// \matrix binary operation expression: matmul(), ...
	template<class T, class U, typename _BinaryOp> 
	class MatBinaryExpr : public Base_<MatBinaryExpr<T, U, _BinaryOp>>
	{
		using _Myt = MatBinaryExpr;
		using _Mybase = Base_<MatBinaryExpr<T, U, _BinaryOp>>;
	public:
		using _Mybase::CompileTimeRows;
		using _Mybase::CompileTimeCols;
		using typename _Mybase::value_t;
		using typename _Mybase::matrix_type;
		using category = category_type_t<_BinaryOp>;
		enum{options = option<_BinaryOp::flag>::value};

		MATRICE_GLOBAL_INL MatBinaryExpr(const T& _lhs, const U& _rhs) 
			noexcept :_Mybase(_lhs.shape()), _LHS(_lhs), _RHS(_rhs) {
			M = _LHS.rows(), K = _LHS.cols(), N = _RHS.cols();
			if (K != _RHS.rows() && K == N) N = _RHS.rows();
		}

		MATRICE_HOST_FINL value_t* operator() (value_t* res) const {
			return _Op(_LHS.data(), _RHS.data(), res, M, K, N);
		}
		MATRICE_GLOBAL_FINL value_t operator() (int r, int c) const {
			return _Op(_LHS, _RHS, r, c);
		}
		MATRICE_GLOBAL_FINL value_t operator() (int _idx) const noexcept {
			int r = _idx / N, c = _idx - r * N;
			return _Op(_LHS, _RHS, r, c);
		}

		template<typename _Mty>
		MATRICE_GLOBAL_INL void assign_to(_Mty& res) const noexcept {
//#pragma omp parallel for if (res.size() > 100)
			for (int i = 0; i < res.size(); ++i) 
				res(i) = this->operator()(i);
		}
	private:
		const T& _LHS; 
		const U& _RHS;
		_BinaryOp _Op;
		using _Mybase::M;
		using _Mybase::K;
		using _Mybase::N;
	};

	// \matrix unary operation expression: inverse(), transpose(), ...
	template<class T, typename _UnaryOp> 
	class MatUnaryExpr : public Base_<MatUnaryExpr<T, _UnaryOp>>
	{
		using _Mybase = Base_<MatUnaryExpr<T, _UnaryOp>>;
	public:
		using _Mybase::CompileTimeRows;
		using _Mybase::CompileTimeCols;
		using typename _Mybase::value_t;
		using typename _Mybase::matrix_type;
		using category = category_type_t<_UnaryOp>;
		enum {options = option<_UnaryOp::flag>::value};

		MATRICE_GLOBAL_INL MatUnaryExpr(const T& inout) noexcept
			:_Mybase(inout.shape()), _RHS(inout), _ANS(inout) {
			if constexpr (options == trp) { 
				std::swap(M, N); 
				_Mybase::Shape.get(2) = M;
				_Mybase::Shape.get(3) = N;
			}
		}
		MATRICE_GLOBAL_INL MatUnaryExpr(const T& _rhs, T& _ans)
			: _Mybase(_rhs.shape()), _RHS(_rhs), _ANS(_ans) {
			if constexpr (options == trp) {
				std::swap(M, N);
				_Mybase::Shape.get(2) = M;
				_Mybase::Shape.get(3) = N;
			}
		}
		MATRICE_GLOBAL_INL MatUnaryExpr(MatUnaryExpr&& _other)
			: _Mybase(_other.shape()), options(_other.options) {
			*this = move(_other);
		}

		MATRICE_GLOBAL_FINL auto operator()() const noexcept { 
			return _Op(M, _ANS.data(), _RHS.data());
		}
		MATRICE_GLOBAL_FINL value_t operator() (int r, int c) const { 
			if constexpr (options == inv) return _ANS(c + r * N);
			if constexpr (options == trp) return _ANS(r + c * N);
		}
		MATRICE_GLOBAL_FINL value_t operator() (int _idx) const { 
			if constexpr (options == trp)
				_idx = _idx * M + _idx / N * (1 - N * M);
			return _ANS(_idx);
		}
		MATRICE_GLOBAL_FINL T& ans() { return _ANS; }
		MATRICE_GLOBAL_FINL const T& ans() const { return _ANS; }
		template<typename _Rhs, 
			typename _Ret = MatBinaryExpr<MatUnaryExpr, _Rhs, Op::_Mat_mul<value_t>>> 
		MATRICE_GLOBAL_INL _Ret mul(const _Rhs& _rhs) {
			if constexpr (options == inv) this->operator()();
			if constexpr (options == trp)                   ;
			return _Ret(*this, _rhs);
		}

		template<typename _Mty>
		MATRICE_GLOBAL_INL void assign_to(_Mty& res) const noexcept {
			if constexpr (options == inv) {
				_Op(M, res.data(), _RHS.data());
			}
			if constexpr (options == trp) {
				const int _Size = this->size();
				if (&res == &_RHS) {
					_Mty _Res(res.shape());
					for (int i = 0; i < _Size; ++i)
							_Res(i) = (*this)(i);
					res = forward<_Mty>(_Res);
				}
				else {
					for (int i = 0; i < _Size; ++i)
						res(i) = (*this)(i);
				}
			}
		} /*i*N+(1-N*M)*(i/M) is replaced by i at left hand*/
	private:
		const T& _RHS;
		const T& _ANS;
		_UnaryOp _Op;
		using _Mybase::M;
		using _Mybase::N;
		using _Mybase::size;
	};

#undef _MATRICE_DEFEXP_ARITHOP
#undef _MATRICE_DEFEXP_EWISEUOP
#undef _MATRICE_DEFEXP_EWISEBOP
};

template<typename _Exp, MATRICE_ENABLE_IF(is_expression_v<_Exp>)>
using _Matrix_exp = typename Expr::Base_<_Exp>;

MATRICE_NAMESPACE_EXPR_END

///<-------------- I am the lovely seperate line -------------->

DGE_MATRICE_BEGIN
using _Exp = exprs::Expr;
using _Exp_op = _Exp::Op;
using _Exp_tag = _Exp::OpFlag;

template<typename _Derived>
struct is_expression<_Exp::Base_<_Derived>>:std::true_type {};
template<typename T, typename U, typename _Op, bool _T, bool _U>
struct is_expression<_Exp::EwiseBinaryExpr<T,U,_Op,_T,_U>>:std::true_type {};
template<typename T, typename _Op>
struct is_expression<_Exp::EwiseUnaryExpr<T,_Op>>:std::true_type {};
template<typename T, typename U, typename _Op>
struct is_expression<_Exp::MatBinaryExpr<T,U,_Op>>:std::true_type {};
template<typename T, typename _Op>
struct is_expression<_Exp::MatUnaryExpr<T,_Op>>:std::true_type {};

template<typename T, typename U, typename _Op>
struct expression_traits<_Exp::EwiseBinaryExpr<T, U, _Op, false, false>> {
	using type = _Exp::EwiseBinaryExpr<T, U, _Op, false, false>;
	using value_type = common_type_t<typename T::value_t, typename U::value_t>;
	enum {
		options = _Exp::option<_Exp_tag::ewise>::value,
		rows = max_integer_v<T::CompileTimeRows, U::CompileTimeRows>,
		cols = max_integer_v<T::CompileTimeCols, U::CompileTimeCols>,
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};
template<typename T, typename U, typename _Op>
struct expression_traits<_Exp::EwiseBinaryExpr<T, U, _Op, true, false>> {
	using type = _Exp::EwiseBinaryExpr<T, U, _Op, true, false>;
	using value_type = typename U::value_t;
	enum {
		options = _Exp::option<_Exp_tag::ewise>::value,
		rows = U::CompileTimeRows,
		cols = U::CompileTimeCols,
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};
template<typename T, typename U, typename _Op>
struct expression_traits<_Exp::EwiseBinaryExpr<T, U, _Op, false, true>> {
	using type = _Exp::EwiseBinaryExpr<T, U, _Op, false, true>;
	using value_type = typename T::value_t;
	enum {
		options = _Exp::option<_Exp_tag::ewise>::value,
		rows = T::CompileTimeRows,
		cols = T::CompileTimeCols,
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};

template<typename T, typename _Op>
struct expression_traits<_Exp::EwiseUnaryExpr<T, _Op>> {
	using type = _Exp::EwiseUnaryExpr<T, _Op>;
	using value_type = typename T::value_t;
	enum {
		options = expression_options<_Op>::value,
		rows = T::CompileTimeRows,
		cols = T::CompileTimeCols
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};
template<typename T, typename U, typename _Op>
struct expression_traits<_Exp::MatBinaryExpr<T, U, _Op>> {
	using type = _Exp::MatBinaryExpr<T, U, _Op>;
	using value_type = common_type_t<typename T::value_t, typename U::value_t>;
	enum {
		options = expression_options<_Op>::value,
		rows = conditional_size_v<T::CompileTimeRows <= 0 || U::CompileTimeRows <= 0, 0, max_integer_v<T::CompileTimeRows, 0>>,
		cols = conditional_size_v<T::CompileTimeCols <= 0 || U::CompileTimeCols <= 0, 0, max_integer_v<U::CompileTimeCols, 0>>
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};
template<typename T, typename _Op>
struct expression_traits<_Exp::MatUnaryExpr<T, _Op>> {
	using type = _Exp::MatUnaryExpr<T, _Op>;
	using value_type = typename T::value_t;
	enum {
		options = expression_options<_Op>::value,
		rows = conditional_size_v<(options&_Exp_tag::trp) == _Exp_tag::trp, T::CompileTimeCols, T::CompileTimeRows>,
		cols = conditional_size_v<(options&_Exp_tag::trp) == _Exp_tag::trp, T::CompileTimeRows, T::CompileTimeCols>
	};
	using auto_matrix_type = types::Matrix_<value_type, rows, cols>;
};

using exprs::Expr;
// *\element-wise addition
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_add<value_t>>>
	MATRICE_GLOBAL_FINL auto operator+ (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise subtraction
template<
	typename _Lhs, class _Rhs,
	typename value_t = conditional_t<is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename     _Op = _Exp::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_sub<value_t>>>
	MATRICE_GLOBAL_FINL auto operator- (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise multiplication
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_mul<value_t>>>
	MATRICE_GLOBAL_FINL auto operator* (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise division
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<is_scalar_v<_Rhs>, _Rhs, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_div<value_t>>>
	MATRICE_GLOBAL_FINL auto operator/ (const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise maximum
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<is_scalar_v<_Rhs>, _Rhs, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_max<value_t>>>
	MATRICE_GLOBAL_FINL auto max(const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }
// *\element-wise minimum
template<
	typename _Lhs, typename _Rhs,
	typename value_t = conditional_t<is_scalar_v<_Rhs>, typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseBinaryExpr<_Lhs, _Rhs, _Exp_op::_Ewise_min<value_t>>>
	MATRICE_GLOBAL_FINL auto min(const _Lhs& _left, const _Rhs& _right) { return _Op(_left, _right); }

// *\element-wise sqr()
template<typename _Rhs>
	MATRICE_GLOBAL_FINL auto sqr(const _Rhs& _right) { return (_right*_right); }

// *\element-wise sqrt()
template<
	typename _Rhs,
	typename value_t = enable_if_t<is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_sqrt<value_t>>>
	MATRICE_GLOBAL_FINL auto sqrt(const _Rhs& _right) { return _Op(_right); }

// *\element-wise exp()
template<
	typename _Rhs,
	typename value_t = enable_if_t<is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_exp<value_t>>>
	MATRICE_GLOBAL_FINL auto exp(const _Rhs& _right) { return _Op(_right); }

// *\element-wise log()
template<
	typename _Rhs,
	typename value_t = enable_if_t<is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_log<value_t>>>
	MATRICE_GLOBAL_FINL auto log(const _Rhs& _right) { return _Op(_right); }

// *\element-wise abs()
template<
	typename _Rhs,
	typename value_t = enable_if_t<is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_abs<value_t>>>
	MATRICE_GLOBAL_FINL auto abs(const _Rhs& _right) { return _Op(_right); }

// *\element-wise floor()
template<
	typename _Rhs,
	typename value_t = enable_if_t<is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename _Op = _Exp::EwiseUnaryExpr<_Rhs, _Exp_op::_Ewise_floor<value_t>>>
	MATRICE_GLOBAL_FINL auto floor(const _Rhs& _right) { return _Op(_right); }

// *\transpose expression
template<
	typename _Rhs,
	typename value_t = enable_if_t<is_scalar_v<typename _Rhs::value_t>, typename _Rhs::value_t>,
	typename _Op = _Exp::MatUnaryExpr<_Rhs, _Exp_op::_Mat_trp<value_t>>>
	MATRICE_GLOBAL_FINL auto transpose(const _Rhs& _right) { return _Op(_right); }

// *\outer product expression : xy^T
template<
	typename _Lhs, typename _Rhs,
	typename value_t = common_type_t<typename _Lhs::value_t, typename _Rhs::value_t>,
	typename _Op = _Exp::MatBinaryExpr<_Lhs, _Rhs, _Exp_op::_Mat_sprmul<value_t>>>
	MATRICE_GLOBAL_FINL auto outer_product(const _Lhs& _left, const _Rhs& _right) {
	return _Op(_left, _right);
}

// *\summation of expression
template<
	typename _Rhs, 
	typename _Op = _Exp_op::_Accum_exp<_Rhs>>
	MATRICE_GLOBAL_FINL auto accum(const _Rhs& _right, _Op&& _op = _Op()) {
	return _op(_right);
}

template<typename _Op>
MATRICE_GLOBAL_INL auto _Exp::Base_<_Op>::var() const {
	const auto _Avg = this->avg();
	auto _Diff = _Exp::EwiseBinaryExpr<
		derived_t, decltype(_Avg),
		_Exp_op::_Ewise_sub<value_t>>(
			*static_cast<const derived_t*>(this), 
			_Avg);
	return (_Diff*_Diff).avg();
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif

DGE_MATRICE_END