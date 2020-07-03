/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

#include <valarray>
#include <functional>
#include "_type_traits.h"
#include "_plain_exp.hpp"
#include "_matrix_ops.hpp"
#include "_storage.hpp"
#include "_iterator.h"
#include "_shape.hpp"
#include "_view.h"
#include "util/_type_defs.h"
#include "util/_conditional_macros.h"
#include "util/_exception.h"
#include "util/_property.hpp"
#include "core/solver.h"
#ifdef MATRICE_USE_OPENCV
#include "../../../addin/interface.h"
#endif

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4514 4786 4503 4127)
#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)
#endif

#ifndef _HAS_CXX17
#  error Matrice library must be compiled as C++ 17 or above.
#endif

DGE_MATRICE_BEGIN

template<typename _Ty> class Scalar;

// \CLASS TEMPLATE to make plane view for all matrix types
template<typename _Ty, int _Type = _View_trait<_Ty>::value>
class _Basic_plane_view_base {
	enum { MAGIC_VAL = 0x42FF0000 };
	size_t stride = type_bytes_v<_Ty>;
	int type = _Type, flags = MAGIC_VAL | _Type; 
	int nval = 1, hval = 1, wval = 1;
public:
	using plvt_type = tuple<int, int, std::add_pointer_t<_Ty>, bool>;

	MATRICE_GLOBAL_FINL constexpr _Basic_plane_view_base() = default;
	MATRICE_GLOBAL_FINL constexpr _Basic_plane_view_base(int _rows, int _cols, _Ty* _data = nullptr) noexcept
		: m_shape(_rows,_cols), m_data(_data) {
#ifdef MATRICE_DEBUG
		this->_Flush_view_buf();
#endif
	};
	MATRICE_GLOBAL_FINL constexpr _Basic_plane_view_base(const shape_t<3>& _Shape, _Ty* _data) noexcept
		: m_data(_data), m_shape(_Shape) {
#ifdef MATRICE_DEBUG
		this->_Flush_view_buf();
#endif
	};
	MATRICE_GLOBAL_FINL constexpr _Basic_plane_view_base(const shape_t<3>& _Shape) noexcept
		: m_shape(_Shape) {
#ifdef MATRICE_DEBUG
		this->_Flush_view_buf();
#endif
	};

	/**
	 * \shape of a matrix or tensor
	 * Example: 
		auto [_Rows, _Cols, _Depth(=1)] = _Matrix.shape();
		auto [_Rows, _Cols, _Depth] = _Tensor.shape();
	 */
	MATRICE_GLOBAL_FINL constexpr decltype(auto) shape() const noexcept {
		return m_shape;
	}
	template<typename _T, MATRICE_ENABLE_IF(is_scalar_v<_T>)>
	MATRICE_GLOBAL_FINL constexpr auto shape(_T _Scale) const noexcept {
		return shape_t<2>(m_rows*_Scale, m_cols*_Scale);
	}
	template<typename _T1, typename _T2 = _T1, 
		MATRICE_ENABLE_IF(is_scalar_v<_T1>&&is_scalar_v<_T2>)>
	MATRICE_GLOBAL_FINL constexpr auto shape(_T1 _Rsf, _T2 _Csf) const noexcept {
		return shape_t<2>(m_rows*_Rsf, m_cols*_Csf);
	}

	/**
	 *\brief Get full dims {N,{C,{H,W}}}
	 */
	//MATRICE_GLOBAL_FINL constexpr auto& dims() const noexcept {
	//	return (m_shape);
	//}

	/**
	 * \a tuple of plane view of the matrix for fast memory access.
	 * \param <_Axis> indicates view type: x for row vector view, y for column vector view, all for original plane view. 
	 * \param [require_t] refers to the view will be transposed or not.
	 */
	template<axis _Axis = axis::all>
	MATRICE_GLOBAL_FINL constexpr plvt_type plvt(bool require_t = 0) noexcept {
		if constexpr (_Axis == axis::x)
			return std::make_tuple(1, m_rows*m_cols, m_data, require_t);
		else if constexpr (_Axis == axis::y)
			return std::make_tuple(m_rows*m_cols, 1, m_data, require_t);
		else
			return std::make_tuple(m_rows, m_cols, m_data, require_t);
	}
	template<axis _Axis = axis::all>
	MATRICE_GLOBAL_FINL constexpr plvt_type plvt(bool require_t = 0) const noexcept {
		if constexpr (_Axis == axis::x)
			return std::make_tuple(1, m_rows*m_cols, m_data, require_t);
		else if constexpr (_Axis == axis::y)
			return std::make_tuple(m_rows*m_cols, 1, m_data, require_t);
		else
			return std::make_tuple(m_rows, m_cols, m_data, require_t);
	}
	MATRICE_GLOBAL_FINL constexpr operator plvt_type() noexcept {
		return std::make_tuple(m_rows, m_cols, m_data);
	}

	/**
	 * \raw plane view of the matrix
	 * Example: auto& _Raw = _Matrix.raw();
	 */
	MATRICE_GLOBAL_FINL constexpr auto& raw() noexcept {
		return *this;
	}
	MATRICE_GLOBAL_FINL constexpr const auto& raw() const noexcept {
		return *this;
	}

protected:
	MATRICE_GLOBAL_FINL constexpr void _Flush_view_buf() noexcept {
		stride = m_cols * type_bytes_v<_Ty>;
		nval = m_shape.d;
		hval = m_shape.h;
		wval = m_shape.w;
	}

	shape_t<3> m_shape;
	size_t m_rows = m_shape.rows();
	size_t m_cols = m_shape.cols();
	_Ty* m_data = nullptr;
};

_DETAIL_BEGIN

/*******************************************************************
	              Generic Base for Matrix Class
	    Copyright (c) : Zhilong (Dgelom) Su, since 14/Feb/2018
 ******************************************************************/
template<
	typename _Derived, 
	typename _Mytraits = matrix_traits<_Derived>, 
	typename _Valty = typename _Mytraits::type>
class Base_ : public _Basic_plane_view_base<_Valty>
{
#pragma region MACRO DEFINITIONS
#define MATRICE_LINK_PTR { \
	m_data = _Myalloc.data(); \
}

#define MATRICE_EVALEXP_TOTHIS { \
	m_data = _Myalloc.data(); \
	exp.assign(*this); \
}

#define MATRICE_EXPAND_SHAPE \
_Shape.h, _Shape.w

#define MATRICE_MAKE_EXPOP_TYPE(DESC, NAME) \
typename _Exp_op::_##DESC##_##NAME<_Valty>

#define MATRICE_MAKE_ARITHOP(OP, NAME) \
template<typename _Rhs> MATRICE_GLOBAL_INL \
auto operator##OP(const _Rhs& _Right) const noexcept { \
	return Exp::EwiseBinaryExp<_Myt, _Rhs, _Xop_ewise_##NAME>(*this, _Right); \
} \
template<typename _Lhs, MATRICE_ENABLE_IF(is_scalar_v<_Lhs>)> friend \
MATRICE_GLOBAL_FINL auto operator##OP(const _Lhs& _Left, const _Derived& _Right) noexcept { \
	return Exp::EwiseBinaryExp<_Lhs, _Derived, _Xop_ewise_##NAME>(_Left, _Right); \
}

#define MATRICE_MAKE_EXP_ASSIGNOP(NAME) \
template<typename _Lhs, typename _Rhs, typename _Op> \
MATRICE_GLOBAL_INL \
auto& operator=(const Exp##NAME##BinaryExp<_Lhs,_Rhs,_Op>& _Ex)noexcept{ \
return (*static_cast<_Derived*>(&_Ex.assign(*this))); \
} \
template<typename _Rhs, typename _Op> \
MATRICE_GLOBAL_INL \
auto& operator=(const Exp##NAME##UnaryExp<_Rhs, _Op>& _Ex) noexcept{ \
return (*static_cast<_Derived*>(&_Ex.assign(*this))); \
}
#pragma endregion

	enum { _M = _Mytraits::_M, _N = _Mytraits::_N };
	using _Mylayout = plain_layout::row_major;
	using _Myalloc_t = Allocator<_Valty, _M, _N, allocator_traits_v<_M, _N>, _Mylayout>;

	using _Myt = Base_;
	using _Mybase = _Basic_plane_view_base<_Valty>;
	using _Myreference = std::add_lvalue_reference_t<_Myt>;
	using _Myt_fwd_iterator = _Matrix_forward_iterator<_Valty>;
	using _Myt_rwise_iterator = _Matrix_rwise_iterator<_Valty>;
	using _Myt_cwise_iterator = _Matrix_cwise_iterator<_Valty>;
	using _Myt_rview_type = _Matrix_rview<_Valty, _N>;
	using _Myt_cview_type = _Matrix_cview<_Valty, _M>;
	using _Myt_blockview_type = _Matrix_block<_Valty>;
	using _Xop_ewise_add   = MATRICE_MAKE_EXPOP_TYPE(Ewise, add);
	using _Xop_ewise_sub   = MATRICE_MAKE_EXPOP_TYPE(Ewise, sub);
	using _Xop_ewise_mul   = MATRICE_MAKE_EXPOP_TYPE(Ewise, mul);
	using _Xop_ewise_div   = MATRICE_MAKE_EXPOP_TYPE(Ewise, div);
	using _Xop_ewise_sqrt  = MATRICE_MAKE_EXPOP_TYPE(Ewise, sqrt);
	using _Xop_ewise_exp   = MATRICE_MAKE_EXPOP_TYPE(Ewise, exp);
	using _Xop_ewise_log   = MATRICE_MAKE_EXPOP_TYPE(Ewise, log);
	using _Xop_ewise_log2  = MATRICE_MAKE_EXPOP_TYPE(Ewise, log2);
	using _Xop_ewise_log10 = MATRICE_MAKE_EXPOP_TYPE(Ewise, log10);
	using _Xop_mat_mul     = MATRICE_MAKE_EXPOP_TYPE(Mat,   mul);
	using _Xop_mat_inv     = MATRICE_MAKE_EXPOP_TYPE(Mat,   inv);
	using _Xop_mat_trp     = MATRICE_MAKE_EXPOP_TYPE(Mat,   trp);
public:
	using base_t = _Mybase;
	using value_t = _Valty;
	using value_type = value_t;
	using scalar_type = Scalar<value_type>;
	using pointer = value_type*;
	using reference = value_type&;
	using iterator = pointer;
	using const_iterator = _Matrix_const_iterator<value_type>;
	using const_initlist = std::add_const_t<initlist<value_t>>;
	using derived_t = _Derived;
	using loctn_t = Location;
	using category = typename _Mytraits::category;
	template<typename _Xop> using exp_base_type = Exp::Base_<_Xop>;

	/**
	 *\brief static properties
	 */
	enum { options = _Allocator_traits<_Myalloc_t>::options };
	/**
	 *\brief static property for querying the rows at compile-time
	 */
	static constexpr long long rows_at_compiletime = _M;
	/**
	 *\brief static property for querying the cols at compile-time
	 */
	static constexpr long long cols_at_compiletime = _N;
	/**
	 *\brief static property for querying the size at compile-time
	 */
	static constexpr long long Size = rows_at_compiletime * cols_at_compiletime;
	/**
	 *\brief for static querying memory location
	 */
	static constexpr auto location = options;
	/**
	 *\brief for querying infinity attribute of the value type
	 */
	static constexpr auto inf = std::numeric_limits<value_type>::infinity();
	/**
	 *\brief for querying round error attribute of the value type
	 */
	static constexpr auto eps = std::numeric_limits<value_type>::epsilon();

	MATRICE_GLOBAL_INL constexpr Base_()noexcept
		:_Mybase(max_integer_v<0, _M>, max_integer_v<0, _N>), _Myalloc() {
		m_data = _Myalloc.data();
	}
	MATRICE_GLOBAL_INL constexpr Base_(size_t _rows)noexcept
		:_Mybase(_rows < _M ? _M : _rows, _N > 1 ? _N : 1),
		_Myalloc(_rows < _M ? _M : _rows, _N > 1 ? _N : 1) {
		m_data = _Myalloc.data();
	}
	MATRICE_GLOBAL_INL constexpr Base_(size_t _rows, size_t _cols)noexcept
		:_Mybase(_rows < _M ? _M : _rows, _cols < _N ? _N : _cols),
		_Myalloc(_rows < _M ? _M : _rows, _cols < _N ? _N : _cols) {
		m_data = _Myalloc.data();
	}
	MATRICE_GLOBAL_INL explicit Base_(size_t _rows, size_t _cols, pointer data)noexcept
		:_Mybase(_rows, _cols, data) {
	}
	MATRICE_GLOBAL_INL explicit Base_(size_t _rows, size_t _cols, value_t _val)noexcept
		:Base_(_rows, _cols) {
		_Myalloc = (_val);
	}
	MATRICE_GLOBAL_INL explicit Base_(shape_t<2>&& _Shape)noexcept
		:Base_(MATRICE_EXPAND_SHAPE) {
	}
	MATRICE_GLOBAL_INL explicit Base_(shape_t<2>&& _Shape, pointer _Data)noexcept
		:_Mybase(MATRICE_EXPAND_SHAPE, _Data) {
	}
	MATRICE_GLOBAL_INL explicit Base_(shape_t<2>&& _Shape, value_t _Val)noexcept
		:Base_(MATRICE_EXPAND_SHAPE) {
		_Myalloc = (_Val);
	}
	MATRICE_GLOBAL_INL explicit Base_(shape_t<3>&& _Shape) noexcept
		:_Mybase(_Shape), _Myalloc(_Shape.rows(), _Shape.cols()) {
		m_data = _Myalloc.data();
	}
	MATRICE_GLOBAL_INL constexpr Base_(const_initlist _list) noexcept
		:Base_(_M > 0 ? _M : _N > 0 ? 1 : _list.size(), _N > 0 ? _N : 1) {
		_Myalloc = ((pointer)_list.begin());
	}
	MATRICE_GLOBAL_INL constexpr Base_(const _Myt& _other) noexcept
		:_Mybase(_other.shape()), _Myalloc(_other._Myalloc) {
		m_data = _Myalloc.data();
	}
	MATRICE_GLOBAL_INL constexpr Base_(_Myt&& _other) noexcept
		:_Mybase(_other.shape()), _Myalloc(move(_other._Myalloc)) {
		m_data = _Myalloc.data();
	}
	/**
	 *\from STD vector<value_t>, while no memory is copied.
	 */
	MATRICE_HOST_INL Base_(std::vector<value_t>&_other, int _cols=1) noexcept
		:Base_(_other.size()/_cols, _cols, _other.data()) {}
	/**
	 *\from STD valarray<...>, while no memory is copied.
	 */
	MATRICE_HOST_INL Base_(std::valarray<value_t>& _other, int _cols = 1) noexcept
		:Base_(_other.size() / _cols, _cols, &_other[0]) {}
	/**
	 *\from explicit specified matrix type
	 */
	template<int _CTR, int _CTC,
		typename _Mty = Matrix_<value_t, _CTR, _CTC>, 
		typename = enable_if_t<is_not_same_v<_Mty, _Derived>>>
	MATRICE_GLOBAL_INL constexpr Base_(const _Mty& _Oth) noexcept
		:_Mybase(_Oth.rows(), _Oth.cols()),_Myalloc(_Oth.allocator()){
		m_data = _Myalloc.data();
	}
#ifdef MATRICE_ENABLE_CUDA
	MATRICE_GLOBAL_INL constexpr Base_(const Matrix_<value_t, ::device> & _Oth) noexcept
		:_Mybase(_Oth.rows(), _Oth.cols()), _Myalloc(_Oth.allocator()) {
		m_data = _Myalloc.data();
	}
#endif
	/**
	 *\from expression
	 */
	template<typename _Exp>
	MATRICE_GLOBAL_FINL Base_(const exprs::_Matrix_exp<_Exp>& exp)
		:Base_(exp.rows(), exp.cols()) MATRICE_EVALEXP_TOTHIS
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Exp::EwiseBinaryExp<_Lhs,_Rhs,_Op>& exp)
		:Base_(exp.shape()) MATRICE_EVALEXP_TOTHIS
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Exp::MatBinaryExp<_Lhs, _Rhs, _Op>& exp)
		:Base_(exp.shape()) MATRICE_EVALEXP_TOTHIS
	template<typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Exp::MatUnaryExp<_Rhs, _Op>& exp)
		:Base_(exp.shape()) MATRICE_EVALEXP_TOTHIS

	/**
	 *\interfaces for opencv if it is enabled
	 */
#ifdef MATRICE_USE_OPENCV
	MATRICE_HOST_INL Base_(cv::Mat&& mat):Base_(mat.rows, mat.cols, mat.ptr<value_t>()){
		mat.flags = 0x42FF0000; mat.dims = mat.rows = mat.cols = 0;
		mat.data = nullptr; mat.datastart = nullptr; mat.dataend = nullptr; mat.datalimit = nullptr;
		mat.allocator = nullptr;
		mat.u = nullptr;
	}
	template<typename _Fn>
	MATRICE_HOST_INL Base_(const cv::Mat& mat, _Fn&& _Op) : Base_(mat) { each(_Op); }
	MATRICE_HOST_INL cv::Mat cvmat() { return cv::Mat(m_rows, m_cols, cv::DataType<value_t>::type, m_data); }
	MATRICE_HOST_INL const cv::Mat cvmat() const { return cv::Mat(m_rows, m_cols, cv::DataType<value_t>::type, m_data); }
#endif
public:
	/**
	 *\create a matrix with dynamic (host or device) memory allocation
	 */
	MATRICE_HOST_INL _Derived& create(diff_t _Rows, diff_t _Cols = (1)) {
		if constexpr (_M <= 0 || _N <= 0) 
			this->derived().__create_impl(_Rows, _Cols);
		return this->derived();
	};
	template<typename _Uy, MATRICE_ENABLE_IF(is_scalar_v<_Uy>)>
	MATRICE_HOST_INL _Derived& create(diff_t _Rows, diff_t _Cols, _Uy _Val) {
		return this->create(_Rows, _Cols) = value_type(_Val);
	};
	MATRICE_HOST_INL _Derived& create(const shape_t<2>& _Shape) {
		if constexpr (_M <= 0 || _N <= 0) 
			this->derived().__create_impl(MATRICE_EXPAND_SHAPE);
		return this->derived();
	};
	template<typename _Uy, MATRICE_ENABLE_IF(is_scalar_v<_Uy>)>
	MATRICE_HOST_INL _Derived& create(const shape_t<2>& _Shape, _Uy _Val) {
		return this->create(_Shape) = value_type(_Val);
	};

	/**
	 *\returns pointer to y-th row
	 *\sa ptr()
	 */
	MATRICE_GLOBAL_FINL pointer operator[](index_t y) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(y < rows(), "Matrix_ subscript out of row range.");
#endif
		return (m_data + y * m_cols); 
	}
	MATRICE_GLOBAL_FINL const pointer operator[](index_t y) const noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(y < rows(), "Matrix_ subscript out of row range.");
#endif
		return (m_data + y * m_cols); 
	}
	/**
	 *\1D index random accessor to get i-th element reference
	 */
	MATRICE_GLOBAL_FINL reference operator()(index_t i) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(i < size(), "Matrix_ subscript out of range.");
#endif
		return m_data[i]; 
	}
	MATRICE_GLOBAL_FINL const reference operator()(index_t i)const noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(i < size(), "Matrix_ subscript out of range.");
#endif
		return m_data[i]; 
	}
	/**
	 *\2D index random accessor to get element reference at r-th row and c-th col.
	 */
	MATRICE_GLOBAL_INL reference operator()(index_t r, index_t c) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(r < rows(), "Matrix_ subscript out of row range.");
		DGELOM_CHECK(c < cols(), "Matrix_ subscript out of column range.");
#endif
		return (*this)[r][c]; 
	}
	MATRICE_GLOBAL_INL const reference operator()(index_t r, index_t c) const noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(r < rows(), "Matrix_ subscript out of row range.");
		DGELOM_CHECK(c < cols(), "Matrix_ subscript out of column range.");
#endif
		return (*this)[r][c]; 
	}

	/**
	 *\brief returns pointer to the object memory
	 */
	MATRICE_GLOBAL_FINL constexpr pointer data()noexcept { return (m_data); }
	MATRICE_GLOBAL_FINL constexpr const pointer data()const noexcept { return (m_data); }

	/**
	 *\brief returns pointer to y-th row
	 *\sa operator[]
	 */
	MATRICE_GLOBAL_FINL pointer ptr(int y = 0) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(y < rows(), "Matrix_ subscript out of row range.");
#endif
		return (m_data + (m_cols) * (y)); 
	}
	MATRICE_GLOBAL_FINL const pointer ptr(int y = 0) const {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(y < rows(), "Matrix_ subscript out of row range.");
#endif
		return (m_data + (m_cols) * (y)); 
	}

	/**
	 *\brief returns a reference to the derived object
	 */
	MATRICE_GLOBAL_INL const _Derived& derived() const noexcept {
		return *static_cast<_Derived*>(this);
	}
	MATRICE_GLOBAL_INL _Derived& derived() noexcept {
		return *static_cast<_Derived*>(this);
	}

	/**
	 *\brief returns a const reference to the derived object
	 */
	MATRICE_GLOBAL_INL const _Derived& const_derived() const noexcept {
		return *static_cast<const _Derived*>(this);
	}

	/**
	 * \brief returns reference to the derived object
	 */
	MATRICE_GLOBAL_FINL _Derived& eval() noexcept {
		return (this->derived());
	}
	/**
	 * \brief returns const reference to the derived object
	 */
	MATRICE_GLOBAL_FINL constexpr const _Derived& eval()const noexcept {
		return (this->derived());
	}

#pragma region <!-- iterators -->
	/**
	 *\brief returns STL-stype element-wise iterator
	 */
	MATRICE_GLOBAL_FINL const_iterator begin()const noexcept {
		return { m_data, size() };
	}
	MATRICE_GLOBAL_FINL const_iterator end()const noexcept {
		return { m_data + size(), size() };
	}
	MATRICE_GLOBAL_FINL iterator begin()noexcept { 
		return (m_data); 
	}
	MATRICE_GLOBAL_FINL iterator end()noexcept { 
		return (m_data + size()); 
	}
	/**
	 *\brief column iterator for accessing elements in i-th column
	 *\example:
	 *		auto _A = Matrix_<float,3,3>::rand();
	 *		auto _Fwd_col = _A.cbegin(1); //get 1-th column iterator
	 *		for(auto& _It : _Fwd_col) _It = float(0); //set this column to zero
	 */
	MATRICE_GLOBAL_FINL _Myt_fwd_iterator cbegin(size_t i) noexcept {
		return _Myt_fwd_iterator(m_data + i, m_rows, m_cols);
	}
	MATRICE_GLOBAL_FINL _Myt_fwd_iterator cend(size_t i) noexcept {
		return _Myt_fwd_iterator(_End(m_data + i, m_rows, m_cols));
	}
	MATRICE_GLOBAL_FINL const _Myt_fwd_iterator cbegin(size_t i) const noexcept {
		return _Myt_fwd_iterator(m_data + i, m_rows, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_fwd_iterator cend(size_t i) const noexcept {
		return _Myt_fwd_iterator(_End(m_data + i, m_rows, m_cols));
	}

	/**
	 *\brief row iterator for accessing elements in i-th row
	 */
	MATRICE_GLOBAL_FINL _Myt_fwd_iterator rbegin(size_t i) noexcept {
		return _Myt_fwd_iterator(m_data + i * m_cols, m_cols);
	}
	MATRICE_GLOBAL_FINL _Myt_fwd_iterator rend(size_t i) noexcept {
		return _Myt_fwd_iterator(_End(m_data + i * m_cols, m_cols));
	}
	MATRICE_GLOBAL_FINL const _Myt_fwd_iterator rbegin(size_t i) const noexcept {
		return _Myt_fwd_iterator(m_data + i * m_cols, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_fwd_iterator rend(size_t i) const noexcept {
		return _Myt_fwd_iterator(_End(m_data + i * m_cols, m_cols));
	}

	/**
	 *\brief column-wise iterator to get all elements in i-th column
	 */
	MATRICE_GLOBAL_FINL _Myt_cwise_iterator cwbegin(size_t i = 0) noexcept {
		return _Myt_cwise_iterator(m_data + i * m_rows, m_cols, m_rows);
	}
	MATRICE_GLOBAL_FINL _Myt_cwise_iterator cwend() noexcept {
		return _Myt_cwise_iterator(_End(m_data, m_cols));
	}
	MATRICE_GLOBAL_FINL const _Myt_cwise_iterator cwbegin(size_t i = 0) const noexcept {
		return _Myt_cwise_iterator(m_data + i * m_rows, m_cols, m_rows);
	}
	MATRICE_GLOBAL_FINL const _Myt_cwise_iterator cwend(size_t i = 0) const noexcept {
		return _Myt_cwise_iterator(_End(m_data, m_cols));
	}

	/**
	 * \row-wise iterator to get all elements in i-th row
	 */
	MATRICE_GLOBAL_FINL _Myt_rwise_iterator rwbegin(size_t i = 0) noexcept {
		return _Myt_rwise_iterator(m_data + i * m_cols, m_rows, m_cols);
	}
	MATRICE_GLOBAL_FINL _Myt_rwise_iterator rwend() noexcept {
		return _Myt_rwise_iterator(_End(m_data, m_rows, m_cols));
	}
	MATRICE_GLOBAL_FINL const _Myt_rwise_iterator rwbegin(size_t i = 0) const noexcept {
		return _Myt_rwise_iterator(m_data + i * m_cols, m_rows, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_rwise_iterator rwend()const noexcept{
		return _Myt_rwise_iterator(_End(m_data, m_rows, m_cols));
	}
#pragma endregion

#pragma region <!-- views -->
	// \View of i-th row 
	MATRICE_GLOBAL_FINL _Myt_rview_type rview(size_t i) noexcept {
		return _Myt_rview_type(m_data + m_cols * i, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_rview_type rview(size_t i) const noexcept {
		return _Myt_rview_type(m_data + m_cols * i, m_cols);
	}
	// \View of i-th column
	MATRICE_GLOBAL_FINL _Myt_cview_type cview(size_t i) noexcept {
		return _Myt_cview_type(m_data + i, m_rows, m_cols, i);
	}
	MATRICE_GLOBAL_FINL const _Myt_cview_type cview(size_t i) const noexcept {
		return _Myt_cview_type(m_data + i, m_rows, m_cols, i);
	}

	// \View of this object.
	MATRICE_GLOBAL_INL const auto view() const noexcept {
		return _Myt_blockview_type(m_data, m_cols, m_rows);
	}
	MATRICE_GLOBAL_INL auto view() noexcept {
		return _Myt_blockview_type(m_data, m_cols, m_rows);
	}

	// \View of submatrix: x \in [x0, x1) and y \in [y0, y1)
	MATRICE_GLOBAL_INL auto block(index_t x0, index_t x1, index_t y0, index_t y1) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(x1<=m_cols, "Input var. 'x1' must not be greater than m_cols.")
		DGELOM_CHECK(y1<=m_rows, "Input var. 'y1' must not be greater than m_rows.")
#endif // _DEBUG
		return _Myt_blockview_type(m_data, m_cols, {x0, y0, x1, y1});
	}
	MATRICE_GLOBAL_INL const auto block(index_t x0, index_t x1, index_t y0, index_t y1) const {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(x1<=m_cols, "Input var. x1 must not be greater than m_cols.")
		DGELOM_CHECK(y1<=m_rows, "Input var. y1 must not be greater than m_rows.")
#endif // _DEBUG
		return _Myt_blockview_type(m_data, m_cols, { x0, y0, x1, y1 });
	}
	template<typename... _Ity, MATRICE_ENABLE_IF(sizeof...(_Ity) == 4)>
	MATRICE_GLOBAL_INL const auto block(const tuple<_Ity...>& _R)const {
		return this->block(get<0>(_R), get<1>(_R), get<2>(_R), get<3>(_R));
	}

	// \View of submatrix: x \in [x0, x1) and y \in [y0, y1)
	template<int _Extent>
	MATRICE_GLOBAL_INL auto block(index_t start, size_t num) {
		if constexpr (_Extent == ::extent_x) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(start + num - 1 <= m_rows, "Over range in the row direction.");
#endif // _DEBUG
		return this->block(0, m_cols, start, start + num);
		}
		if constexpr (_Extent == ::extent_y) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(start + num - 1 <= m_cols, "Over range in the col direction.");
#endif // _DEBUG
		return this->block(start, start + num, 0, m_rows);
		}
	}

	/** 
	 * \brief View of a square submatrix.
	 * \param [_Cx, _Cy]: central pos, _Rs: radius size 
	 */
	template<typename _Ity, MATRICE_ENABLE_IF(is_integral_v<_Ity>)>
	MATRICE_GLOBAL_INL auto block(_Ity _Cx, _Ity _Cy, size_t _Rs = 0) const {
		return this->block(_Cx - _Rs, _Cx + _Rs + 1, _Cy - _Rs, _Cy + _Rs + 1);
	}
	/**
	 * \brief View of a square submatrix.
	 * \param [_Ctr]: central pos, which type _Cty must has forward iterator begin(), _Rs: radius size
	 */
	template<typename _Cty>
	MATRICE_GLOBAL_INL auto block(const _Cty& _Ctr, int _Rs = 0)const {
		return this->block(*_Ctr.begin(), *(_Ctr.begin() + 1), _Rs);
	}
	template<typename _Ity, MATRICE_ENABLE_IF(is_integral_v<_Ity>)>
	MATRICE_HOST_INL auto block(const initlist<_Ity>& _Ctr, int _Rs)const {
		return this->block(*_Ctr.begin(), *(_Ctr.begin() + 1), _Rs);
	}

	/**
	 * \operator to get a block view of this matrix from a given range.
	 */
	template<typename _Ity, MATRICE_ENABLE_IF(is_integral_v<_Ity>)>
	MATRICE_GLOBAL_INL auto operator()(_Ity _L, _Ity _R, _Ity _U, _Ity _D) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_R<=m_cols, "Input var. _R must be no greater than m_cols.")
		DGELOM_CHECK(_D<=m_rows, "Input var. _D must be no greater than m_rows.")
#endif // _DEBUG
		return _Myt_blockview_type(m_data, m_cols, { _L, _U, _R, _D });
	}
	template<typename _Ity, MATRICE_ENABLE_IF(is_integral_v<_Ity>)>
	MATRICE_GLOBAL_INL auto operator()(_Ity _L, _Ity _R, _Ity _U, _Ity _D)const{
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_R<=m_cols, "Input var. _R must be no greater than m_cols.")
		DGELOM_CHECK(_D<=m_rows, "Input var. _D must be no greater than m_rows.")
#endif // _DEBUG
		return _Myt_blockview_type(m_data, m_cols, { _L, _U, _R, _D });
	}
	/**
	 * \operator to get a block view of this matrix from a given tupled range.
	 */
	template<typename... _Ity, MATRICE_ENABLE_IF(sizeof...(_Ity) == 4)>
	MATRICE_GLOBAL_INL auto operator()(const tuple<_Ity...>& _R)const {
		return this->operator()(get<0>(_R), get<1>(_R), get<2>(_R), get<3>(_R));
	}
#pragma endregion

	/**
	 * \size properties
	 */
	MATRICE_GLOBAL_FINL constexpr auto(rows)()const noexcept { return m_rows; }
	MATRICE_GLOBAL_FINL constexpr auto(cols)()const noexcept { return m_cols; }
	MATRICE_GLOBAL_FINL constexpr auto(size)()const noexcept { return m_rows*m_cols; }

	/**
	 * \assignment operator, fill Matrix_ from a scalar.
	 */
	MATRICE_GLOBAL_FINL _Derived& operator=(value_t _Val) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Myalloc, "This object is empty.");
#endif // MATRICE_DEBUG
		_Myalloc = (_Val);
		return (this->derived());
	}

	/**
	 * \assignment operator, fill Matrix_ from a scalar.
	 */
	MATRICE_REQUIRES(is_scalar_v<scalar_type>)
	MATRICE_GLOBAL_FINL _Derived& operator=(scalar_type _Val) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Myalloc, "This object is empty.");
#endif // MATRICE_DEBUG
		_Myalloc = value_type(_Val);
		return (this->derived());
	}

	/**
	 * \assignment operator, fill Matrix_ from initializer list
	 */
	MATRICE_GLOBAL_FINL _Derived& operator=(const_initlist _list) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Myalloc, "This object is empty.");
#endif // MATRICE_DEBUG
		_Myalloc = (const pointer)(_list.begin());
		return (this->derived());
	}
	/**
	 * \assignment operator, from nested initializer list
	 */
	template<typename _Ty>
	MATRICE_HOST_INL _Derived& operator=(nested_initlist<_Ty> _list)noexcept{
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Myalloc, "This object is empty.");
		DGELOM_CHECK(_list.size() == m_rows, "Inconsistent rows.");
		DGELOM_CHECK(_list.begin().size() == m_cols, "Inconsistent cols");
#endif // MATRICE_DEBUG
		size_t _Count = 0;
		for (const auto L : _list) {
			auto _Data = m_data + m_cols * _Count++;
			for (auto It = L.begin(); It != L.end(); ++It) {
				*(_Data++) = *It;
			}
		}
		return (this->derived());
	}
	/**
	 * \assignment operator, from row-wise iterator
	 */
	MATRICE_GLOBAL_FINL _Derived& operator=(const _Myt_rwise_iterator& _It)noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Myalloc, "This object is empty.");
#endif // MATRICE_DEBUG
		std::copy(_It.begin(), _It.end(), _Myalloc.data());
		return (this->derived());
	}
	/**
	 * \assignment operator, from column-wise iterator
	 */
	MATRICE_GLOBAL_FINL _Derived& operator=(const _Myt_cwise_iterator& _It)noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Myalloc, "This object is empty.");
#endif // MATRICE_DEBUG
		std::copy(_It.begin(), _It.end(), _Myalloc.data());
		return (this->derived());
	}
	/**
	 *\assignment operator=()
	 *\brief copy from block view with O(min(this->size(), _Bv.size())
	 */
	template<typename _Ty>
	MATRICE_GLOBAL_FINL _Derived& operator=(const _Matrix_block<_Ty>& _Bv) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Myalloc, "This object is empty.");
#endif // MATRICE_DEBUG
		_Bv.eval_to(*this);
		return (this->derived());
	}
	/**
	 * \homotype copy assignment operator
	 */
	MATRICE_GLOBAL_INL _Derived& operator=(const _Derived& _other) {
		if (this != &_other) {
			m_cols = _other.m_cols;
			m_rows = _other.m_rows;
			m_shape = _other.shape();
			if (_other._Myalloc) {
				m_data = _Myalloc = _other._Myalloc;
			}
			else {
				m_data = _other.data();
			}
#ifdef MATRICE_DEBUG
			_Mybase::_Flush_view_buf();
#endif // MATRICE_DEBUG
		}
		return (this->derived());
	}
	/**
	 * \homotype move assignment operator
	 */
	MATRICE_GLOBAL_INL _Derived& operator=(_Derived&& _other) noexcept {
		if (this != &_other) {
			m_cols = _other.m_cols;
			m_rows = _other.m_rows;
			m_shape = _other.shape();
			if (_other._Myalloc) {
				if constexpr (Size > 0) {
					//copy if they allocated on the stack
					m_data = _Myalloc = (_other._Myalloc);
				}
				else {
					m_data = _Myalloc = move(_other._Myalloc);
					_other._Myalloc.destroy();
				}
			}
			else {
				m_data = _other.data();
			}
#ifdef MATRICE_DEBUG
			_Mybase::_Flush_view_buf();
#endif // MATRICE_DEBUG
		}
		return (this->derived());
	}
	/**
	 * \try to convert managed matrix to dynamic matrix
	 */
	template<int _Rows, int _Cols,
		typename _Maty = Matrix_<value_t, _Rows, _Cols>>
	MATRICE_GLOBAL_FINL _Derived& operator=(_Maty&& _managed) {
		_Myalloc = _managed().allocator();
		m_rows = _Myalloc.rows();
		m_cols = _Myalloc.cols();
		m_data = _Myalloc.data();
#ifdef MATRICE_DEBUG
		_Mybase::_Flush_view_buf();
#endif // MATRICE_DEBUG
		return (this->derived());
	}
	/**
	 *\brief Check if this equals to _other or not
	 *\param [_other] can be any derived type of matrix/array/vector 
	 */
	MATRICE_GLOBAL_INL bool operator==(const _Myt& _other) const noexcept {
		if(size() != _other.size()) return std::false_type::value;
		if(m_data == _other.m_data) return std::true_type::value;
		for (auto _Idx = 0; _Idx < size(); ++_Idx)
			if (abs(m_data[_Idx] - _other.m_data[_Idx]) > eps)
				return std::false_type::value;
		return std::true_type::value;
	}

#pragma region <!-- Lazied Operators for Matrix Arithmetic -->
	MATRICE_MAKE_ARITHOP(+, add);
	MATRICE_MAKE_ARITHOP(-, sub);
	MATRICE_MAKE_ARITHOP(*, mul);
	MATRICE_MAKE_ARITHOP(/, div);

	template<typename _Rhs> 
	MATRICE_GLOBAL_INL auto mul(const _Rhs& _Right)const noexcept { 
		return Exp::MatBinaryExp<_Myt, _Rhs, _Xop_mat_mul>(*this, _Right);
	}
	MATRICE_GLOBAL_FINL auto sqrt()const noexcept { 
		return Exp::EwiseUnaryExp<_Myt, _Xop_ewise_sqrt>(*this); 
	}
	MATRICE_HOST_FINL auto inv()const noexcept { 
		return Exp::MatUnaryExp<_Myt, _Xop_mat_inv>(*this); 
	}
	MATRICE_HOST_FINL auto inv(const _Myt& _Right)const noexcept {
		return Exp::MatUnaryExp<_Myt, _Xop_mat_inv>(_Right, *this);
	}
	MATRICE_HOST_FINL auto transpose()const noexcept { 
		return Exp::MatUnaryExp<_Myt, _Xop_mat_trp>(*this); 
	}
	MATRICE_GLOBAL_INL auto t()const noexcept {
		return Exp::MatUnaryExp<_Myt, _Xop_mat_trp>(*this);
	}
	MATRICE_GLOBAL_FINL auto normalize(value_t _val = inf)const noexcept { 
		return ((*this)*(abs(_val) < eps ? 1 : 1 / (_val == inf ? max() : _val))); 
	}

	MATRICE_MAKE_EXP_ASSIGNOP(::Ewise);
	MATRICE_MAKE_EXP_ASSIGNOP(::Mat);
#pragma endregion

	///<brief> in-time matrix arithmetic </brief>
	MATRICE_GLOBAL_FINL auto(max)()const noexcept { 
		return (*std::max_element(begin(), end())); 
	}
	MATRICE_GLOBAL_FINL auto(min)()const noexcept {
		return (*std::min_element(begin(), end())); 
	}
	MATRICE_GLOBAL_FINL auto(sum)()const noexcept {
		return (reduce(begin(), end())); 
	}
	MATRICE_GLOBAL_FINL auto(det)() const { 
		return (det_impl(*static_cast<const _Derived*>(this))); 
	}
	MATRICE_GLOBAL_FINL auto(trace)()const noexcept {
		return (reduce(begin().stride(cols()+1), end().stride(cols()+1), 1));
	}

	/**
	 * \matrix Frobenius norm
	 */
	MATRICE_GLOBAL_FINL value_type norm_2()const noexcept {
		auto _Ans = dot(*this); 
		return (_Ans > eps ? ::sqrt(_Ans) : inf); 
	}
	/**
	 * \matrix p-norm: $[\sum_{i=1}^{m}\sum_{j=1}^{}|a_{ij}|^p]^{1/p}$
	 * Sepcial cases: $p = 0$ for $\infty$-norm, $p = 1$ for 1-norm and $p = 2$ for 2-norm
	 * \Reference: https://en.wikipedia.org/wiki/Matrix_norm
	 */
	template<size_t _P = 2> 
	MATRICE_GLOBAL_FINL value_type norm() const {
		return internal::_Matrix_norm_impl<_P>::value(*(this));
	}
	/**
	 * \dot product of this matrix with _Rhs
	 */
	template<typename _Rhs> 
	MATRICE_GLOBAL_FINL value_type dot(const _Rhs& _rhs)const noexcept {
		return this->operator*(_rhs).sum(); 
	}
	/**
	 * \contraction of two matrices (rank-2 tensors)
	 //tex: $\text{res} = \mathbf{A}:\mathbf{B}$
	 */
	template<int Rows, int Cols>
	MATRICE_GLOBAL_INL value_type contract(const Matrix_<value_type, Rows, Cols>& _Rhs)const;

	/**
	 *\brief in-place instant subtraction
	 *\param [_Right] can be scalar or any compatible types
	 */
	template<typename _Rhs>
	MATRICE_HOST_INL auto sub_inplace(const _Rhs& _Right);

	/**
	 *\brief in-place matrix-vector multiplication. Note that if the number of rows of this matrix A equals to the size of the right column vector x, this method returns (A^T)x.
	 *\param [_Right] will be unrolled to a column vector x if it is not.
	 */
	template<typename _Rhs>
	MATRICE_HOST_INL auto mv_inplace(const _Rhs& _Right) const;

	/**
	 * \brief in-place matmul with _Rhs.
	 */
	template<typename _Rhs = _Derived>
	MATRICE_HOST_INL auto mul_inplace(const _Rhs& _Right) const;

	/**
	 *\brief spread to element-wisely multiplicate with an input
	 *\param [_Right] input argument with a type of _Rhs.
	 */
	template<typename _Rhs>
	MATRICE_GLOBAL_INL _Rhs spreadmul(const _Rhs& _Right) const;

	/**
	 * \brief operate each entry via _Fn
	 */
	template<typename _Op>
	MATRICE_GLOBAL_FINL _Derived& each(_Op&& _Fn) noexcept {
		for (auto& _Val : *this) _Fn(_Val); 
		return(this->derived());
	}

	/**
	 * \brief operate each entry via _Fn in parallel.
	 */
	template<typename _Ewop>
	MATRICE_HOST_INL _Derived& parallel_each(_Ewop&& _Fn) {
#pragma omp parallel for
		for (diff_t _Idx = 0; _Idx < size(); ++_Idx) {
			m_data[_Idx] = _Fn(m_data[_Idx]);
		}
		return(this->derived());
	}

	/**
	 * \brief ref to another instance without malloc and copy.
	 */
	template<typename _Src, MATRICE_ENABLE_IF(Size==::dynamic)>
	MATRICE_GLOBAL_INL _Derived& ref(_Src& _other) noexcept {
		if (this->data() != _other.data()) {
			m_cols = _other.cols(), m_rows = _other.rows();
			m_shape = _other.shape();
			m_data = _other.data();
			_Mybase::_Flush_view_buf();
		}
		return (this->derived());
	}

	/**
	 * \brief copy from another data block
	 */
	template<typename _It>
	MATRICE_GLOBAL_FINL _Derived& from(const _It _Data) {
		using arg_type = remove_all_t<decltype(*_Data)>;
		if constexpr (std::is_convertible_v<arg_type, _Derived>&&!is_iterator_v<_It>) {
			*this = *_Data;
		}
		else {
#ifdef MATRICE_DEBUG
			DGELOM_CHECK(_Data + this->size() - 1, "Input length of _Data must be greater or equal to this->size().");
#endif
			for (diff_t _Idx = 0; _Idx < size(); ++_Idx)
				m_data[_Idx] = static_cast<value_type>(_Data[_Idx]);
		}
		return (this->derived());
	}
	/**
	 * \brief convert from another data block by function _Fn
	 */
	template<typename _It, typename _Op, 
		MATRICE_ENABLE_IF(is_iterator_v<_It >)>
	MATRICE_GLOBAL_FINL _Derived& from(const _It _Data, _Op&& _Fn) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_Data + this->size() - 1, "Input length of _Data must be greater or equal to this->size().");
#endif
#pragma omp parallel for if (size()>1000)
		for(diff_t _Idx = 0; _Idx < size(); ++_Idx)
			m_data[_Idx] = _Fn(static_cast<value_type>(_Data[_Idx]));

		return (this->derived());
	}
	/**
	 * \brief Stack from a sequence of vectors with same size, _Vecty can be any type that has members .size() and .data()
	 * \param [_L] the input _Vecty-typed vector list
	 * \param [_Dim] = 0 for row-wise stacking, = 1 for column-wise stacking 
	 */
	template<typename _Vecty>
	MATRICE_HOST_INL _Derived& stack_from(const initlist<_Vecty>& _L, size_t _Dim = 0) {
		const auto _Rows = _Dim == 0 ? _L.size() : _L.begin()->size();
		const auto _Cols = _Dim == 0 ? _L.begin()->size() : _L.size();

		if (this->empty) this->create(_Rows, _Cols);

		if (_Dim == 0) {
			for (auto _Idx = 0; _Idx < _Rows; ++_Idx)
				this->rview(_Idx) = (_L.begin() + _Idx)->data();
		}
		else if (_Dim == 1) {
			for (auto _Idx = 0; _Idx < _Cols; ++_Idx)
				this->cview(_Idx) = (_L.begin() + _Idx)->data();
		}
		else 
			DGELOM_ERROR("The _Dim value should be 0(row-wise) or 1(col-wise).");

		return (this->derived());
	}

	/**
	 *\brief convert to a vector
	 */
	MATRICE_GLOBAL_INL auto vec() const noexcept {
		Matrix_<value_type, Size, 1> _Ret(size());
		size_t _Idx = 0;
		for (auto _Col = cwbegin(); _Col != cwend(); ++_Col) {
			for (auto _It = _Col.begin(); _It != _Col.end(); ++_It) {
				_Ret(_Idx++) = *_It;
			}
		}
		return _Ret;
	}

	/**
	 *\brief Replace entries meets _Cond with _Val
	 *\param [_Cond] the condition function
	 *\param [_Val] the value for replacement
	 */
	MATRICE_GLOBAL_INL void where(std::function<bool(const value_type&)> _Cond, const value_type _Val) {
		this->each([&](auto& _My_val) {
			_My_val = _Cond(_My_val) ? _Val : _My_val; });
	}
	/**
	 *\brief Check if all elements are contained in a range
	 *\param [_Lower, _Upper] the range boundary
	 */
	MATRICE_GLOBAL_INL bool in_range(const _Myt& _Lower, const _Myt& _Upper) const {
		for (const auto idx : range(0, size())) {
			if (m_data[idx] < _Lower.m_data[idx] || m_data[idx] > _Upper.m_data[idx])
				return std::false_type::value;
		}
		return std::true_type::value;
	}

	/**
	 *\brief set to the identity matrix
	 *\param [_Size] _Size-by-_Size matrix
	 */
	MATRICE_GLOBAL_INL _Derived& identity(diff_t _Size=0) noexcept {
		if(empty) this->create(_Size, _Size, 0);
		else this->operator= ((value_type)(0));
		for (auto _Idx = 0; _Idx < rows(); ++_Idx)
			this->operator[](_Idx)[_Idx] = one<value_type>;
		return (this->derived());
	}

	/**
	 *\brief Perfom Cholesky decomposition:
	  //tex:$A = {L}{L^T}$.
	 *\return A proxy to spd instance which holds the "L" part of this matrix, so that the L part can be retrieved with "decltype (A) = A.spd()" and the inverse matrix is obtained with "auto Inv = A.spd().inv()". 
	 */
	MATRICE_GLOBAL_INL auto spd() noexcept {
		return detail::_Matrix_fact<_Derived, tag::_Linear_spd_tag>(this->derived());
	}

	/**
	 *\brief Create a zero-value filled matrix
	 *\param [_Rows, Cols] height and width of matrix, only specified for dynamic created matrix
	 */
	static MATRICE_GLOBAL_INL _Derived zeros(diff_t _Rows=0, diff_t _Cols=0) {
		_Derived _Ret(_Rows, _Cols);
		return forward<_Derived>(_Ret = zero<value_type>);
	}

	/**
	 *\brief Create a diagonal square matrix
	 *\param [_Val] diagonal element value, 
	 *\param [_Size] matrix size, only specified for dynamic case
	 */
	template<typename _Uy, MATRICE_ENABLE_IF(is_scalar_v<_Uy>)>
	static MATRICE_GLOBAL_INL _Derived diag(_Uy _Val = 1, diff_t _Size = 0) {
		_Derived _Ret(_Size, _Size);
		_Ret = value_type(0);
		const auto _Value = value_type(_Val);
		for (auto _Idx = 0; _Idx < _Ret.rows(); ++_Idx) 
			_Ret[_Idx][_Idx] = _Value;
		return forward<_Derived>(_Ret);
	}
	/**
	 *\brief creates a matrix filled by real numbers of uniform distribution.
	 *\param [_Rows, Cols] the height and width of a matrix, only needed to specify for dynamic memory alloc cases.
	 */
	static MATRICE_GLOBAL_INL _Derived rand(diff_t _Rows=0, diff_t _Cols=0) {
		_Derived _Ret(_Rows, _Cols);
		uniform_real<value_type> _Rand; mt19937 _Eng;
		return forward<_Derived>(_Ret.each([&](auto& _Val) {
			_Val = _Rand(_Eng); }));
	}

	/**
	 *\brief creates a matrix filled by random numbers of normal distribution.
	 *\param [_Pars] = {Mean, STD}
	 */
	MATRICE_REQUIRES(rows_at_compiletime>0&&cols_at_compiletime>0)
	static MATRICE_HOST_INL _Derived randn(const_initlist& _Pars = {/*mean=*/0, /*STD=*/1}) {
		_Derived _Ret;
		normal_distribution<value_type> _Rand{ *_Pars.begin(), *(_Pars.begin() + 1) };
		mt19937 _Eng;
		return forward<_Derived>(_Ret.each([&](auto& _Val) {
			_Val = _Rand(_Eng); }));
	}
	static MATRICE_HOST_INL _Derived randn(diff_t _Extent, const_initlist& _Pars = {/*mean=*/0, /*STD=*/1 }) {
		_Derived _Ret(_Extent);
		normal_distribution<value_type> _Rand{ *_Pars.begin(), *(_Pars.begin() + 1) };
		mt19937 _Eng;
		return forward<_Derived>(_Ret.each([&](auto& _Val) {
			_Val = _Rand(_Eng); }));
	}
	static MATRICE_HOST_INL _Derived randn(diff_t _Rows, diff_t _Cols, const_initlist& _Pars = {/*mean=*/0, /*STD=*/1 }) {
		_Derived _Ret(_Rows, _Cols);
		normal_distribution<value_type> _Rand{ *_Pars.begin(), *(_Pars.begin() + 1) };
		mt19937 _Eng;
		return forward<_Derived>(_Ret.each([&](auto& _Val) {
			_Val = _Rand(_Eng); }));
	}

	template<class _Mty, typename _Op, MATRICE_ENABLE_IF(is_matrix_v<_Mty>)>
	static MATRICE_GLOBAL_INL _Derived make(_Mty&& m) {
		_Derived _Ret(m.shape());
		return forward<_Derived>(_Ret.from(m.data()));
	}
	template<class _Mty, typename _Op, MATRICE_ENABLE_IF(is_matrix_v<_Mty>)>
	static MATRICE_GLOBAL_INL _Derived make(const _Mty& m, _Op&& op) {
		_Derived _Ret(m.shape());
		return forward<_Derived>(_Ret.from(m.data(), op));
	}

	MATRICE_GET(bool, empty, !bool(size()));
	MATRICE_GET(size_t, format, _Myfmt);

protected:
	using _Mybase::m_rows;
	using _Mybase::m_cols;
	using _Mybase::m_data;
	using _Mybase::m_shape;

	size_t _Myfmt = _Myalloc.fmt() | gene;
	conditional_t<(rows_at_compiletime>=::dynamic||cols_at_compiletime >=::dynamic), _Myalloc_t, typename detail::Storage_<value_type>::template Allocator<rows_at_compiletime, cols_at_compiletime>> _Myalloc;

	MATRICE_GLOBAL_INL _Myt& _Xfields(initlist<size_t> il) noexcept {
		m_shape = il;
		m_rows = m_shape.rows();
		m_cols = m_shape.cols();
		m_data = _Myalloc.data();
		_Mybase::_Flush_view_buf();
		return (*this);
	}

	MATRICE_GLOBAL_INL _Myt& _Buy(initlist<size_t> il) noexcept {
		_Xfields(il);
		m_data = _Myalloc.alloc(m_rows, m_cols).data();
		return (*this);
	}

public:
	MATRICE_HOST_INL decltype(auto) allocator() noexcept {
		return (_Myalloc);
	}
	MATRICE_HOST_INL decltype(auto) allocator()const noexcept {
		return (_Myalloc);
	}
	MATRICE_HOST_INL decltype(auto) deleter() noexcept {
		return (_Myalloc.deleter());
	}

#undef MATRICE_MAKE_EXPOP_TYPE
#undef MATRICE_LINK_PTR
#undef MATRICE_EVALEXP_TOTHIS
#undef MATRICE_MAKE_ARITHOP
#undef MATRICE_MAKE_EXP_ASSIGNOP
}; 

struct _Matrix_padding {
	template<typename _Ty, int _M, int _N>
	using _Matrix_t = detail::Matrix_<_Ty, _M, _N>;

	/**
	 *\brief Make zero padding
	 *\param <_B> templated padding size, [_In] input matrix.
	 */
	template<size_t _B, typename _Mty, typename _Ty = typename _Mty::value_t> 
	MATRICE_GLOBAL_INL static auto zero(const _Mty& _In) {
		static_assert(is_matrix_v<_Mty>, "_Mty must be a matrix type.");
		
		constexpr size_t _M = _Mty::rows_at_compiletime;
		constexpr size_t _N = _Mty::cols_at_compiletime;
		_Matrix_t<_Ty, _M + (_B << 1), _N + (_B << 1)> _Ret;
		_Ret.block(_B, size_t(_In.cols()) + _B, _B, size_t(_In.rows()) + _B) = _In;

		return forward<decltype(_Ret)>(_Ret);
	}

	template<typename _Mty, typename _Ty = typename _Mty::value_t>
	MATRICE_GLOBAL_INL static auto zero(const _Mty & _In, size_t _B) {
		static_assert(is_matrix_v<_Mty>, "_Mty must be a matrix type.");

		_Mty _Ret(_In.rows()+ (_B << 1), _In.cols()+(_B<<1), 0);
		_Ret.block(_B, _In.cols() + _B, _B, _In.rows() + _B) = _In;

		return forward<decltype(_Ret)>(_Ret);
	}

	template<typename _Mty, typename _Ty = typename _Mty::value_t>
	MATRICE_GLOBAL_INL static auto zero(const _Mty & _In, size_t _LU, size_t _RB) {
		static_assert(is_matrix_v<_Mty>, "_Mty must be a matrix type.");

		_Mty _Ret(_In.rows() + _LU+_RB, _In.cols() + _LU+_RB, 0);
		_Ret.block(_LU, _In.cols() + _LU, _LU, _In.rows() + _LU) = _In;

		return forward<decltype(_Ret)>(_Ret);
	}
};
_DETAIL_END
using matrix_padding =  detail::_Matrix_padding;

template<typename _Mty>
MATRICE_HOST_INL decltype(auto) make_matrix_deleter(const _Mty& _M) noexcept;

/**
 *\brief dgelom::Matrix_<...> factory function to create matrix object.
 *\param [params] wrap any parameters supported by dgelom::Matrix_ ctors, while except initializer list. 
 */
template<typename _Ty, int _Rows=0, int _Cols=0, typename... _Args>
MATRICE_GLOBAL_INL decltype(auto) make_matrix(_Args&&... params);

/**
 *\func dgelom::make_zero<_Ty>(_Ty&)
 *\brief Set input data to zero.
 *\param [data] can be scalar, dgelom::Matrix_, dgelom::Tensor or any container with methods data() and size().
 */
template<typename _Ty>
MATRICE_GLOBAL_INL remove_all_t<_Ty>& make_zero(_Ty& data) noexcept;

/**
 *\func dgelom::view<_Mty>(_Mty&)
 *\brief Make view of the given matrix or vector.
 *\param [_M] can be dgelom::Matrix_ or dgelom::Vec_.
 */
template<typename _Mty, MATRICE_ENABLE_IF(is_matrix_v<_Mty>||is_fxdvector_v<_Mty>)>
MATRICE_GLOBAL_FINL auto view(_Mty& _M) noexcept;

/**
 *\func dgelom::swap<_Mty>(_Mty&, _Mty&)
 *\brief Swap two matrices
 */
template<typename _Mty>
MATRICE_HOST_INL void swap(_Mty& _L, _Mty& _R) noexcept;

/**
 *\func dgelom::copy<_Mty>(_Mty&, _Mty&)
 *\brief Copy a given matrix. Always wrap a dynamic matrix with the function If a deep copy is required.
 */
template<typename _Mty, MATRICE_ENABLE_IF(is_matrix_v<_Mty>)>
MATRICE_HOST_INL _Mty copy(const _Mty& _M);

DGE_MATRICE_END

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "inl/_base.inl"