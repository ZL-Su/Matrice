/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include "_view.h"
#include "../../../addin/interface.h"
#include "../util/_type_defs.h"
#include "../util/_macro_conditions.h"
#include "../util/_exception.h"
#include "../util/_property.hpp"
#include "../core/solver.h"

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
	int nval = 1, cval = 1, hval = 1, wval = 1;
public:
	using plvt_type = tuple<int, int, std::add_pointer_t<_Ty>, bool>;

	MATRICE_GLOBAL_FINL constexpr _Basic_plane_view_base() = default;
	MATRICE_GLOBAL_FINL constexpr _Basic_plane_view_base(int _rows, int _cols, _Ty* _data = nullptr) noexcept
		: m_shape{1,1,_rows,_cols} {
		this->_Flush_view_buf();
	};
	MATRICE_GLOBAL_FINL constexpr _Basic_plane_view_base(const basic_shape_t& _Shape, _Ty* _data) noexcept
		: m_data(_data), m_shape(_Shape) {
		this->_Flush_view_buf();
	};
	MATRICE_GLOBAL_FINL constexpr _Basic_plane_view_base(const basic_shape_t& _Shape) noexcept
		: m_shape(_Shape) {
		this->_Flush_view_buf();
	};

	/**
	 * \shape of a matrix or tensor
	 * Example: 
		auto [E, D, _Rows, _Cols] = _Matrix.shape();
		auto [E, D, H, W] = _Tensor.shape();
	 */
	MATRICE_GLOBAL_FINL constexpr decltype(auto) shape() const noexcept {
		return m_shape;
	}
	template<typename _T, MATRICE_ENABLE_IF(is_scalar_v<_T>)>
	MATRICE_GLOBAL_FINL constexpr shape_t<size_t>  shape(_T _Scale) const noexcept {
		return std::make_tuple(m_rows*_Scale, m_cols*_Scale);
	}
	template<typename _T1, typename _T2 = _T1, 
		MATRICE_ENABLE_IF(is_scalar_v<_T1>&&is_scalar_v<_T2>)>
	MATRICE_GLOBAL_FINL constexpr shape_t<size_t> shape(_T1 _Rsf, _T2 _Csf) const noexcept {
		return std::make_tuple(m_rows*_Rsf, m_cols*_Csf);
	}

	/**
	 *\brief Get full dims {N,{C,{H,W}}}
	 */
	MATRICE_GLOBAL_FINL constexpr auto& dims() const noexcept {
		return (m_shape);
	}

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
#ifdef _DEBUG || MATRICE_DEBUG
		stride = m_cols * type_bytes_v<_Ty>;
		nval = m_shape.get(0);
		cval = m_shape.get(1);
		hval = m_shape.get(2);
		wval = m_shape.get(3);
#endif
	}

	basic_shape_t m_shape;
	int m_rows = m_shape.rows();
	int m_cols = m_shape.cols();
	_Ty* m_data = nullptr;
};

_TYPES_BEGIN

/*******************************************************************
	              Generic Base for Matrix Class
	    Copyright (c) : Zhilong (Dgelom) Su, since 14/Feb/2018
 ******************************************************************/
template<
	typename _Derived, 
	typename _Mytraits = matrix_traits<_Derived>, 
	typename _Valtype = typename _Mytraits::type>
class Base_ : public _Basic_plane_view_base<_Valtype>
{
#define MATRICE_LINK_PTR { \
	m_data = m_storage.data(); \
}
#define MATRICE_EVALEXP_TOTHIS { \
	m_data = m_storage.data(); \
	exp.assign(*this); \
}
#define MATRICE_EXPAND_SHAPE get<0>(_Shape), get<1>(_Shape)
#define MATRICE_MAKE_EXPOP_TYPE(DESC, NAME) typename _Exp_op::_##DESC##_##NAME<_Valtype>
#define MATRICE_MAKE_ARITHOP(OP, NAME) \
template<typename _Rhs> MATRICE_GLOBAL_INL \
auto operator##OP (const _Rhs& _Right) const { \
	return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_##NAME>(*this, _Right); \
} \
template<typename _Lhs, MATRICE_ENABLE_IF(is_scalar_v<_Lhs>)> friend \
MATRICE_GLOBAL_FINL auto operator##OP(const _Lhs& _Left, const _Derived& _Right) { \
	return Expr::EwiseBinaryExpr<_Lhs, _Derived, _Xop_ewise_##NAME>(_Left, _Right); \
}
#define MATRICE_MAKE_EXP_ASSIGNOP(NAME) \
template<typename _Lhs, typename _Rhs, typename _Op> \
MATRICE_GLOBAL_INL \
auto& operator= (const Expr##NAME##BinaryExpr<_Lhs, _Rhs, _Op>& _Ex){ \
return (*static_cast<_Derived*>(&_Ex.assign(*this))); \
} \
template<typename _Rhs, typename _Op> \
MATRICE_GLOBAL_INL \
auto& operator= (const Expr##NAME##UnaryExpr<_Rhs, _Op>& _Ex){ \
return (*static_cast<_Derived*>(&_Ex.assign(*this))); \
}

	struct _My_element {
		_My_element(_Valtype& _Val, size_t _Idx, size_t _Std) 
			: _Value(_Val), _Index(_Idx), _Stride(_Std) {}
		MATRICE_GLOBAL_INL operator _Valtype&() { return _Value; }
		MATRICE_GLOBAL_INL operator _Valtype&() const { return _Value; }
		MATRICE_GLOBAL_INL auto index() const { return _Index; }
		MATRICE_GLOBAL_INL auto x() const { return (_Index - y()*_Stride); }
		MATRICE_GLOBAL_INL auto y() const { return (_Index / _Stride); }
	private:
		_Valtype& _Value;
		size_t _Index, _Stride;
	};
	enum { _M = _Mytraits::_M, _N = _Mytraits::_N };
	using _Myalty = typename detail::Storage_<_Valtype>::template Allocator<_M, _N>;
	using _Myt = Base_;
	using _Mybase = _Basic_plane_view_base<_Valtype>;
	using _Myt_const = std::add_const_t<_Myt>;
	using _Myt_reference = std::add_lvalue_reference_t<_Myt>;
	using _Myt_const_reference = std::add_lvalue_reference_t<_Myt_const>;
	using _Myt_fwd_iterator = _Matrix_forward_iterator<_Valtype>;
	using _Myt_rwise_iterator = _Matrix_rwise_iterator<_Valtype>;
	using _Myt_cwise_iterator = _Matrix_cwise_iterator<_Valtype>;
	using _Myt_rview_type = _Matrix_rview<_Valtype>;
	using _Myt_cview_type = _Matrix_cview<_Valtype>;
	using _Myt_blockview_type = _Matrix_block<_Valtype>;
	using _Xop_ewise_add   = MATRICE_MAKE_EXPOP_TYPE(Ewise, add);
	using _Xop_ewise_sub   = MATRICE_MAKE_EXPOP_TYPE(Ewise, sub);
	using _Xop_ewise_mul   = MATRICE_MAKE_EXPOP_TYPE(Ewise, mul);
	using _Xop_ewise_div   = MATRICE_MAKE_EXPOP_TYPE(Ewise, div);
	using _Xop_ewise_sqrt  = MATRICE_MAKE_EXPOP_TYPE(Ewise, sqrt);
	using _Xop_ewise_exp   = MATRICE_MAKE_EXPOP_TYPE(Ewise, exp);
	using _Xop_ewise_log   = MATRICE_MAKE_EXPOP_TYPE(Ewise, log);
	using _Xop_ewise_log2  = MATRICE_MAKE_EXPOP_TYPE(Ewise, log2);
	using _Xop_ewise_log10 = MATRICE_MAKE_EXPOP_TYPE(Ewise, log10);
	using _Xop_mat_mul     = MATRICE_MAKE_EXPOP_TYPE(Mat, mul);
	using _Xop_mat_inv     = MATRICE_MAKE_EXPOP_TYPE(Mat, inv);
	using _Xop_mat_trp     = MATRICE_MAKE_EXPOP_TYPE(Mat, trp);
public:
	using value_t = _Valtype;
	using value_type = value_t;
	using derived_t = _Derived;
	using base_t = _Mybase;
	using pointer = std::add_pointer_t<value_t>;
	using reference = std::add_lvalue_reference_t<value_t>;
	using iterator = pointer;
	using const_iterator = std::add_const_t<iterator>;
	using const_initlist = std::add_const_t<initlist<value_t>>;
	using loctn_t = Location;
	using category = typename _Mytraits::category;
	template<typename _Xop> 
	using expr_type = Expr::Base_<_Xop>;
	using scalar_type = Scalar<value_type>;
	
	enum { options = _Myalty::location };
	enum { Size = _M*_N, CompileTimeRows = _M, CompileTimeCols = _N, };
	/**
	 *\brief for static querying memory location
	 */
	static constexpr auto location = _Myalty::location;
	/**
	 *\brief for static querying memory block rows
	 */
	static constexpr auto RowsAtCT = CompileTimeRows;
	/**
	 *\brief for static querying memory block cols
	 */
	static constexpr auto ColsAtCT = CompileTimeCols;
	/**
	 *\brief for querying infinity attribute of the value type
	 */
	static constexpr auto inf = std::numeric_limits<value_t>::infinity();
	/**
	 *\brief for querying round error attribute of the value type
	 */
	static constexpr auto eps = std::numeric_limits<value_t>::epsilon();

	MATRICE_GLOBAL_INL constexpr Base_() noexcept
		:_Mybase(_M<0?0:_M, _N<0?0:_N), m_storage() MATRICE_LINK_PTR
	MATRICE_GLOBAL_INL constexpr Base_(int _rows, int _cols) noexcept
		:_Mybase(_rows, _cols), m_storage(_rows, _cols) MATRICE_LINK_PTR
	MATRICE_GLOBAL_INL constexpr Base_(const shape_t<size_t>& _Shape) noexcept
		: Base_(MATRICE_EXPAND_SHAPE) {}
	MATRICE_GLOBAL_INL constexpr Base_(int _rows, int _cols, pointer data) noexcept
		:_Mybase(_rows, _cols, data), m_storage(_rows, _cols, data) {}
	MATRICE_GLOBAL_INL constexpr Base_(const shape_t<size_t>& _Shape, pointer _Data) noexcept
		: Base_(MATRICE_EXPAND_SHAPE, _Data) {}
	MATRICE_GLOBAL_INL constexpr Base_(int _rows, int _cols, value_t _val) noexcept
		:_Mybase(_rows, _cols), m_storage(_rows, _cols, _val) MATRICE_LINK_PTR
	MATRICE_GLOBAL_INL constexpr Base_(const shape_t<size_t>& _Shape, value_t _Val) noexcept
		:Base_(MATRICE_EXPAND_SHAPE, _Val) {}
	MATRICE_GLOBAL_INL constexpr Base_(const basic_shape_t& _Shape) noexcept
		:_Mybase(_Shape), m_storage(_Shape.rows(), _Shape.cols())
		MATRICE_LINK_PTR
	MATRICE_GLOBAL_INL constexpr Base_(const_initlist _list) noexcept
		:_Mybase((_M <= 0)?_list.size():_M, (_N<=0)?1:_N), m_storage(_list) MATRICE_LINK_PTR
	MATRICE_GLOBAL_INL constexpr Base_(const _Myt& _other) noexcept
		:_Mybase(_other.m_shape), m_storage(_other.m_storage) MATRICE_LINK_PTR
	MATRICE_GLOBAL_INL constexpr Base_(_Myt&& _other) noexcept
		:_Mybase(_other.m_shape), m_storage(move(_other.m_storage)) MATRICE_LINK_PTR
	/**
	 *\from STD vector<value_t>
	 */
	MATRICE_HOST_INL Base_(const std::vector<value_t>&_other, int _cols=1) noexcept
		:Base_(_other.size()/_cols, _cols) { from(_other.data()); }
	/**
	 *\from STD valarray<...>
	 */
	MATRICE_HOST_INL Base_(const std::valarray<value_t>& _other, int _rows = 1) noexcept 
		:_Mybase(_rows, _other.size()/_rows), m_storage(_rows, _other.size()/_rows, (pointer)std::addressof(_other[0])) MATRICE_LINK_PTR
	/**
	 *\from explicit specified matrix type
	 */
	template<int _CTR, int _CTC, MATRICE_ENABLE_IF((_CTR>0&&_CTC>0))>
	MATRICE_GLOBAL_INL constexpr Base_(const Matrix_<value_t, _CTR, _CTC>& _Oth) noexcept
		:_Mybase(_Oth.rows(), _Oth.cols()),m_storage(_Oth.allocator())
		MATRICE_LINK_PTR
	MATRICE_GLOBAL_INL constexpr Base_(const Matrix_<value_t,-1,-1>& _Oth) noexcept
		:_Mybase(_Oth.rows(), _Oth.cols()), m_storage(_Oth.allocator())
		MATRICE_LINK_PTR
	/**
	 *\from expression
	 */
	template<typename _Exp>
	MATRICE_GLOBAL_FINL Base_(const exprs::_Matrix_exp<_Exp>& exp)
		: _Mybase(m_rows = exp.rows(), m_cols = exp.cols()), 
		m_storage(m_rows, m_cols) MATRICE_EVALEXP_TOTHIS
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::EwiseBinaryExpr<_Lhs,_Rhs,_Op>& exp)
		:_Mybase(/*m_rows = exp.rows(), m_cols = exp.cols()*/exp.shape()),
		m_storage(m_rows, m_cols) MATRICE_EVALEXP_TOTHIS
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::MatBinaryExpr<_Lhs, _Rhs, _Op>& exp)
		:_Mybase(m_rows = exp.rows(), m_cols = exp.cols()), 
		m_storage(m_rows, m_cols) MATRICE_EVALEXP_TOTHIS
	template<typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::MatUnaryExpr<_Rhs, _Op>& exp)
		:_Mybase(m_rows = exp.rows(), m_cols = exp.cols()), 
		m_storage(m_rows, m_cols) MATRICE_EVALEXP_TOTHIS

	/**
	 *\from nested initializer list {{...},{...},...,{...}}
	 */
	template<typename _Ty>
	MATRICE_HOST_INL constexpr Base_(const nested_initlist<_Ty> nil)
		: Base_(nil.size(), nil.begin()->size()) {
		size_t _Count = 0;
		for (const auto It : nil) {
			auto _Data = m_data + m_cols * _Count++;
			for (auto i = 0; i < m_cols; ++i) {
				_Data[i] = It.begin()[i];
			}
		}
	}

	/**
	 *\interfaces for opencv if it is enabled
	 */
#ifdef MATRICE_USE_OPENCV
	//MATRICE_HOST_INL Base_(const ocv_view_t& mat) : _Mybase(mat.rows, mat.cols, mat.ptr<value_t>()) {}
	MATRICE_HOST_INL Base_(ocv_view_t&& mat):Base_(mat.rows, mat.cols, mat.ptr<value_t>()){
		mat.flags = 0x42FF0000; mat.dims = mat.rows = mat.cols = 0;
		mat.data = nullptr; mat.datastart = nullptr; mat.dataend = nullptr; mat.datalimit = nullptr;
		mat.allocator = nullptr;
		mat.u = nullptr;
	}
	template<typename _Fn>
	MATRICE_HOST_INL Base_(const ocv_view_t& mat, _Fn&& _Op) : Base_(mat) { each(_Op); }
	MATRICE_HOST_INL ocv_view_t cvmat() { return ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data); }
	MATRICE_HOST_INL const ocv_view_t cvmat() const { return ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data); }
#endif
public:
	/**
	 *\create a matrix with dynamic (host or device) memory allocation
	 */
	MATRICE_HOST_INL decltype(auto)create(diff_t _Rows, diff_t _Cols = (1)) {
		if constexpr (_M <= 0 && _N <= 0) 
			this->derived().__create_impl(_Rows, _Cols);
		return (*static_cast<_Derived*>(this));
	};
	template<typename _Uy, MATRICE_ENABLE_IF(is_scalar_v<_Uy>)>
	MATRICE_HOST_INL decltype(auto)create(diff_t _Rows, diff_t _Cols, _Uy _Val) {
		this->create(_Rows, _Cols);
		return (*(this) = value_type(_Val));
	};
	MATRICE_HOST_INL decltype(auto)create(const shape_t<size_t>& _Shape) {
		if constexpr (_M <= 0 && _N <= 0) 
			this->derived().__create_impl(MATRICE_EXPAND_SHAPE);
		return (*static_cast<_Derived*>(this));
	};
	template<typename _Uy, MATRICE_ENABLE_IF(is_scalar_v<_Uy>)>
	MATRICE_HOST_INL decltype(auto)create(const shape_t<size_t>& _Shape, _Uy _Val) {
		this->create(_Shape);
		return (*(this) = value_type(_Val));
	};

	/**
	 *\returns pointer to y-th row
	 *\sa ptr()
	 */
	MATRICE_GLOBAL_FINL pointer operator[](index_t y) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(y < rows(), "Matrix_ subscript out of row range.");
#endif
		return (m_data + y * m_cols); 
	}
	MATRICE_GLOBAL_FINL const pointer operator[](index_t y) const {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(y < rows(), "Matrix_ subscript out of row range.");
#endif
		return (m_data + y * m_cols); 
	}
	/**
	 *\1D index random accessor to get i-th element reference
	 */
	MATRICE_GLOBAL_FINL reference operator()(index_t i) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(i < size(), "Matrix_ subscript out of range.");
#endif
		return m_data[i]; 
	}
	MATRICE_GLOBAL_FINL const reference operator()(index_t i) const {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(i < size(), "Matrix_ subscript out of range.");
#endif
		return m_data[i]; 
	}
	/**
	 *\2D index random accessor to get element reference at r-th row and c-th col.
	 */
	MATRICE_GLOBAL_INL reference operator()(index_t r, index_t c) {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(r < rows(), "Matrix_ subscript out of row range.");
		DGELOM_CHECK(c < cols(), "Matrix_ subscript out of column range.");
#endif
		return (*this)[r][c]; 
	}
	MATRICE_GLOBAL_INL const reference operator()(index_t r, index_t c) const {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(r < rows(), "Matrix_ subscript out of row range.");
		DGELOM_CHECK(c < cols(), "Matrix_ subscript out of column range.");
#endif
		return (*this)[r][c]; 
	}

	/**
	 *\returns pointer to object memory
	 */
	MATRICE_GLOBAL_FINL constexpr pointer data() noexcept { return (m_data); }
	MATRICE_GLOBAL_FINL constexpr const pointer data() const noexcept { return (m_data); }
	/**
	 *\returns pointer to y-th row
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
	 * \returns reference to the derived object
	 */
	MATRICE_GLOBAL_FINL _Derived& eval() noexcept {
		return (this->derived());
	}
	/**
	 * \returns const reference to the derived object
	 */
	MATRICE_GLOBAL_FINL constexpr const _Derived& eval()const noexcept {
		return (this->derived());
	}

#pragma region <!-- iterators -->
	/**
	 *\returns STL-stype iterator
	 */
	MATRICE_GLOBAL_FINL iterator begin() noexcept { 
		return (m_data); 
	}
	MATRICE_GLOBAL_FINL iterator end() noexcept { 
		return (m_data + size()); 
	}
	MATRICE_GLOBAL_FINL const iterator begin() const noexcept { 
		return (m_data); 
	}
	MATRICE_GLOBAL_FINL const iterator end() const noexcept { 
		return (m_data + size()); 
	}
	/**
	 *\column iterator for accessing elements in i-th column
	 *\example:
	 *		auto _A = Matrix_<float,3,3>::rand();
	 *		auto _Fwd_col = _A.cbegin(1); //get 1-th column iterator
	 *		for(auto& _It : _Fwd_col) _It = float(0); //set this column to zero
	 */
	MATRICE_GLOBAL_FINL _Myt_fwd_iterator cbegin(size_t i) {
		return _Myt_fwd_iterator(m_data + i, m_rows, m_cols);
	}
	MATRICE_GLOBAL_FINL _Myt_fwd_iterator cend(size_t i) {
		return _Myt_fwd_iterator(_End(m_data + i, m_rows, m_cols));
	}
	MATRICE_GLOBAL_FINL const _Myt_fwd_iterator cbegin(size_t i) const {
		return _Myt_fwd_iterator(m_data + i, m_rows, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_fwd_iterator cend(size_t i) const {
		return _Myt_fwd_iterator(_End(m_data + i, m_rows, m_cols));
	}

	/**
	 *\row iterator for accessing elements in i-th row
	 */
	MATRICE_GLOBAL_FINL _Myt_fwd_iterator rbegin(size_t i) {
		return _Myt_fwd_iterator(m_data + i * m_cols, m_cols);
	}
	MATRICE_GLOBAL_FINL _Myt_fwd_iterator rend(size_t i) {
		return _Myt_fwd_iterator(_End(m_data + i * m_cols, m_cols));
	}
	MATRICE_GLOBAL_FINL const _Myt_fwd_iterator rbegin(size_t i) const {
		return _Myt_fwd_iterator(m_data + i * m_cols, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_fwd_iterator rend(size_t i) const {
		return _Myt_fwd_iterator(_End(m_data + i * m_cols, m_cols));
	}

	/**
	 *\column-wise iterator for accessing elements
	 */
	MATRICE_GLOBAL_FINL _Myt_cwise_iterator cwbegin(size_t i = 0) {
		return _Myt_cwise_iterator(m_data + i * m_rows, m_cols, m_rows);
	}
	MATRICE_GLOBAL_FINL _Myt_cwise_iterator cwend() {
		return _Myt_cwise_iterator(_End(m_data, m_cols));
	}
	MATRICE_GLOBAL_FINL const _Myt_cwise_iterator cwbegin(size_t i = 0) const {
		return _Myt_cwise_iterator(m_data + i * m_rows, m_cols, m_rows);
	}
	MATRICE_GLOBAL_FINL const _Myt_cwise_iterator cwend(size_t i = 0) const {
		return _Myt_cwise_iterator(_End(m_data, m_cols));
	}

	/**
	 * \row-wise iterator for accessing elements
	 */
	MATRICE_GLOBAL_FINL _Myt_rwise_iterator rwbegin(size_t i = 0) {
		return _Myt_rwise_iterator(m_data + i * m_cols, m_rows, m_cols);
	}
	MATRICE_GLOBAL_FINL _Myt_rwise_iterator rwend() {
		return _Myt_rwise_iterator(_End(m_data, m_rows, m_cols));
	}
	MATRICE_GLOBAL_FINL const _Myt_rwise_iterator rwbegin(size_t i = 0) const {
		return _Myt_rwise_iterator(m_data + i * m_cols, m_rows, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_rwise_iterator rwend() const {
		return _Myt_rwise_iterator(_End(m_data, m_rows, m_cols));
	}
#pragma endregion

#pragma region <!-- views -->
	// \View of i-th row 
	MATRICE_GLOBAL_FINL _Myt_rview_type rview(size_t i) {
		return _Myt_rview_type(m_data + m_cols * i, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_rview_type rview(size_t i) const {
		return _Myt_rview_type(m_data + m_cols * i, m_cols);
	}
	// \View of i-th column
	MATRICE_GLOBAL_FINL _Myt_cview_type cview(size_t i) {
		return _Myt_cview_type(m_data + i, m_rows, m_cols, i);
	}
	MATRICE_GLOBAL_FINL const _Myt_cview_type cview(size_t i) const {
		return _Myt_cview_type(m_data + i, m_rows, m_cols, i);
	}
	// \View of submatrix: x \in [x0, x1) and y \in [y0, y1)
	MATRICE_GLOBAL_INL auto block(index_t x0, index_t x1, index_t y0, index_t y1) {
#ifdef _DEBUG
		DGELOM_CHECK(x1<=m_cols, "Input var. x1 must be no greater than m_cols.")
		DGELOM_CHECK(y1<=m_rows, "Input var. y1 must be no greater than m_rows.")
#endif // _DEBUG
		return _Myt_blockview_type(m_data, m_cols, {x0, y0, x1, y1});
	}
	MATRICE_GLOBAL_INL const auto block(index_t x0, index_t x1, index_t y0, index_t y1) const {
#ifdef _DEBUG
		DGELOM_CHECK(x1<=m_cols, "Input var. x1 must be no greater than m_cols.")
		DGELOM_CHECK(y1<=m_rows, "Input var. y1 must be no greater than m_rows.")
#endif // _DEBUG
		return _Myt_blockview_type(m_data, m_cols, { x0, y0, x1, y1 });
	}
	template<typename... _Ity, MATRICE_ENABLE_IF(sizeof...(_Ity) == 4)>
	MATRICE_GLOBAL_INL const auto block(const tuple<_Ity...>& _R) const {
		return this->block(get<0>(_R), get<1>(_R), get<2>(_R), get<3>(_R));
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
#ifdef _DEBUG
		DGELOM_CHECK(_R<=m_cols, "Input var. _R must be no greater than m_cols.")
		DGELOM_CHECK(_D<=m_rows, "Input var. _D must be no greater than m_rows.")
#endif // _DEBUG
		return _Myt_blockview_type(m_data, m_cols, { _L, _U, _R, _D });
	}
	template<typename _Ity, MATRICE_ENABLE_IF(is_integral_v<_Ity>)>
	MATRICE_GLOBAL_INL auto operator()(_Ity _L, _Ity _R, _Ity _U, _Ity _D)const{
#ifdef _DEBUG
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
	MATRICE_GLOBAL_FINL _Derived& operator= (value_t _Val) noexcept {
		m_storage = _Val;
		m_data = internal::_Proxy_checked(m_storage.data());
		return (this->derived());
	}

	/**
	 * \assignment operator, fill Matrix_ from initializer list
	 */
	MATRICE_GLOBAL_FINL _Derived& operator= (const_initlist _list) noexcept {
		m_storage = _list;
		m_data = internal::_Proxy_checked(m_storage.data());
		return (this->derived());
	}
	/**
	 * \assignment operator, from nested initializer list
	 */
	template<typename _Ty>
	MATRICE_HOST_INL _Derived& operator=(nested_initlist<_Ty> _list) noexcept {
#ifdef MATRICE_DEBUG
		DGELOM_CHECK(_list.size() == m_rows, "Inconsistent rows.");
		DGELOM_CHECK(_list.begin().size() == m_cols, "Inconsistent cols");
#endif
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
	MATRICE_GLOBAL_FINL _Derived& operator= (const _Myt_rwise_iterator& _It) {
		std::copy(_It.begin(), _It.end(), m_storage.data());
		m_data = internal::_Proxy_checked(m_storage.data());
		return (this->derived());
	}
	/**
	 * \assignment operator, from column-wise iterator
	 */
	MATRICE_GLOBAL_FINL _Derived& operator=(const _Myt_cwise_iterator& _It){
		std::copy(_It.begin(), _It.end(), m_storage.data());
		m_data = internal::_Proxy_checked(m_storage.data());
		return (this->derived());
	}
	/**
	 * \homotype copy assignment operator
	 */
	MATRICE_GLOBAL_INL _Derived& operator=(const _Myt& _other) {
		m_cols = _other.m_cols, m_rows = _other.m_rows;
		m_format = _other.m_format;
		m_storage = _other.m_storage;
		m_data = internal::_Proxy_checked(m_storage.data());
		_Mybase::_Flush_view_buf();

		return (this->derived());
	}
	/**
	 * \homotype move assignment operator
	 */
	MATRICE_GLOBAL_INL _Derived& operator=(_Myt&& _other) noexcept {
		m_cols = _other.m_cols, m_rows = _other.m_rows;
		m_format = _other.m_format;
		m_storage = move(_other.m_storage);
		m_data = internal::_Proxy_checked(m_storage.data());

		std::swap(_Mybase::m_shape, _other.m_shape);
		_Mybase::_Flush_view_buf();

		_other.m_storage.free();
		_other.m_data = nullptr;
		_other.m_cols = _other.m_rows = 0;

		return (this->derived());
	}
	/**
	 * \try to convert managed matrix to dynamic matrix
	 */
	template<int _Rows, int _Cols,
		typename _Mty = Matrix_<value_t, _Rows, _Cols>,
		typename = enable_if_t<is_static_v<_Rows,_Cols>&&!is_static_v<_M,_N>>>
	MATRICE_GLOBAL_FINL _Derived& operator= (_Mty& _managed) {
		m_data = _managed.data();
		m_rows = _managed.rows(), m_cols = _managed.cols();
		m_storage.owner() = detail::Storage_<value_t>::Proxy;
		_Mybase::_Flush_view_buf();
		return (this->derived());
	}
	/**
	 *\brief Check if this equals to _other or not
	 *\param [_other] can be any derived type of matrix/array/vector 
	 */
	MATRICE_GLOBAL_INL bool operator== (const _Myt& _other) const noexcept {
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
	MATRICE_GLOBAL_INL auto mul(const _Rhs& _Right) const { 
		return Expr::MatBinaryExpr<_Myt, _Rhs, _Xop_mat_mul>(*this, _Right);
	}
	MATRICE_GLOBAL_FINL auto sqrt() const { 
		return Expr::EwiseUnaryExpr<_Myt, _Xop_ewise_sqrt>(*this); 
	}
	MATRICE_HOST_FINL auto inv() const { 
		return Expr::MatUnaryExpr<_Myt, _Xop_mat_inv>(*this); 
	}
	MATRICE_HOST_FINL auto inv(_Myt_const_reference _Right) {
		return Expr::MatUnaryExpr<_Myt, _Xop_mat_inv>(_Right, *this);
	}
	MATRICE_HOST_FINL auto transpose() const { 
		return Expr::MatUnaryExpr<_Myt, _Xop_mat_trp>(*this); 
	}
	MATRICE_GLOBAL_INL auto t() const {
		return Expr::MatUnaryExpr<_Myt, _Xop_mat_trp>(*this);
	}
	MATRICE_GLOBAL_FINL auto normalize(value_t _val = inf) const { 
		return ((*this)*(abs(_val) < eps ? 1 : 1 / (_val == inf ? max() : _val))); 
	}

	MATRICE_MAKE_EXP_ASSIGNOP(::Ewise);
	MATRICE_MAKE_EXP_ASSIGNOP(::Mat);
#pragma endregion

	///<brief> in-time matrix arithmetic </brief>
	MATRICE_GLOBAL_FINL auto(max)() const { 
		return (*std::max_element(begin(), end())); 
	}
	MATRICE_GLOBAL_FINL auto(min)() const { 
		return (*std::min_element(begin(), end())); 
	}
	MATRICE_GLOBAL_FINL auto(sum)() const { 
		return (reduce(begin(), end())); 
	}
	MATRICE_GLOBAL_FINL auto(det)() const { 
		return (det_impl(*static_cast<const _Derived*>(this))); 
	}
	MATRICE_GLOBAL_FINL auto(trace)() const { 
		return (reduce(begin(), end(), cols() + 1)); 
	}

	/**
	 * \matrix Frobenius norm
	 */
	MATRICE_GLOBAL_FINL auto norm_2()->value_type const { 
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
	MATRICE_GLOBAL_FINL value_type dot(const _Rhs& _Rhs) const {
		return this->operator*(_Rhs).sum(); 
	}

	/**
	 *\brief in-place instant subtraction
	 *\param [_Right] can be scalar or any compatible types
	 */
	template<typename _Rhs>
	MATRICE_HOST_INL decltype(auto)inplace_sub(const _Rhs& _Right);
	/**
	 *\brief instant subtraction
	 *\param [_Right] can be scalar or any compatible types
	 */
	template<typename _Rhs>
	MATRICE_HOST_INL auto sub_(const _Rhs& _Right);
	/**
	 *\brief in-place instant matrix-vector multiplication
	 *\param [_Right] can be a matrix or a vector types
	 */
	template<typename _Rhs, MATRICE_ENABLE_IF(is_fxdvector_v<_Rhs>)>
	MATRICE_HOST_INL decltype(auto)mul_(const _Rhs& _Right);
	/**
	 *\brief in-place matrix-vector multiplication. Note that if the number of rows of this matrix A equals to the size of the right column vector x, this method returns (A^T)x.
	 *\param [_Right] will be unrolled to a column vector x if it is not.
	 */
	template<typename _Rhs>
	MATRICE_HOST_INL auto inplace_mv(const _Rhs& _Right) const;

	/**
	 * \in-place matmul with _Rhs. 
	 */
	template<ttag _Ltag = ttag::N, ttag _Rtag = ttag::N, 
		typename _Rhs = _Derived, MATRICE_ENABLE_IF(is_matrix_v<_Rhs>)>
	MATRICE_HOST_FINL auto inplace_mul(const _Rhs& _Right);
	/**
	 *\brief spread to element-wisely multiplicate with an input
	 *\param [_Right] input argument with a type of _Rhs.
	 */
	template<typename _Rhs>
	MATRICE_GLOBAL_INL _Rhs spreadmul(const _Rhs& _Right) const;

	/**
	 * \operate each entry via _Fn
	 */
	template<typename _Op>
	MATRICE_GLOBAL_FINL _Derived& each(_Op&& _Fn) noexcept {
		for (auto& _Val : *this) _Fn(_Val); 
		return(this->derived());
	}

	/**
	 * \copy from another data block
	 */
	template<typename _It>
	MATRICE_GLOBAL_FINL _Derived& from(const _It _Data) {
		using arg_type = remove_all_t<decltype(*_Data)>;
		if constexpr (std::is_convertible_v<arg_type, _Derived>) {
			*this = *_Data;
		}
		else {
#ifdef _DEBUG
			DGELOM_CHECK(_Data + this->size() - 1, "Input length of _Data must be greater or equal to this->size().");
#endif // _DEBUG
			for (auto _Idx = 0; _Idx < size(); ++_Idx)
				m_data[_Idx] = static_cast<value_type>(_Data[_Idx]);
		}

		return (this->derived());
	}
	/**
	 * \convert from another data block by function _Fn
	 */
	template<typename _It, typename _Op>
	MATRICE_GLOBAL_FINL _Derived& from(const _It _Data, _Op&& _Fn) {
#ifdef _DEBUG
		DGELOM_CHECK(_Data + this->size() - 1, "Input length of _Data must be greater or equal to this->size().");
#endif // _DEBUG
		for(auto _Idx = 0; _Idx < size(); ++_Idx)
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
	 *\brief Create a zero-value filled matrix
	 *\param [_Rows, Cols] height and width of matrix, only specified for dynamic created matrix 
	 */
	static MATRICE_GLOBAL_INL auto zero(diff_t _Rows = 0, diff_t _Cols = 0) {
		_Derived _Ret;
		return forward<_Derived>(_Ret.create(_Rows, _Cols, 0));
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
	 *\brief Create a diagonal square matrix
	 *\param [_Val] diagonal element value, 
	 *\param [_Size] matrix size, only specified for dynamic case
	 */
	template<typename _Uy, MATRICE_ENABLE_IF(is_scalar_v<_Uy>)>
	static MATRICE_GLOBAL_INL auto diag(_Uy _Val = 1, diff_t _Size = 0) {
		_Derived _Ret; _Ret.create(_Size, _Size, 0);
		const auto _Value = value_type(_Val);
		for (auto _Idx = 0; _Idx < _Ret.rows(); ++_Idx) 
			_Ret[_Idx][_Idx] = _Value;
		return forward<_Derived>(_Ret);
	}
	/**
	 *\brief Create random value filled matrix
	 *\param [_Rows, Cols] height and width of matrix, only specified for dynamic case
	 */
	static MATRICE_GLOBAL_INL auto rand(diff_t _Rows=0, diff_t _Cols=0) {
		_Derived _Ret; _Ret.create(_Rows, _Cols);
		uniform_real<value_type> _Rand; mt19937 _Eng;
		return forward<_Derived>(_Ret.each([&](auto& _Val) {
			_Val = _Rand(_Eng); }));
	}
	static MATRICE_GLOBAL_INL auto randn(const_initlist& _Pars = {/*mean=*/0, /*STD=*/1}, diff_t _Rows = 0, diff_t _Cols = 0) {
		_Derived _Ret; _Ret.create(_Rows, _Cols);
		normal_distribution<value_type> _Rand{ *_Pars.begin(), *(_Pars.begin() + 1) };
		mt19937 _Eng;
		return forward<_Derived>(_Ret.each([&](auto& _Val) {
			_Val = _Rand(_Eng); }));
	}

	///<brief> properties </brief>
	__declspec(property(get=_Format_getter, put=_Format_setter))size_t format;
	MATRICE_HOST_FINL size_t _Format_getter() const { return m_format; }
	MATRICE_HOST_FINL void _Format_setter(size_t format) { m_format = rmaj|format; }

	MATRICE_GET(bool, empty, !bool(size()));

protected:
	using _Mybase::m_rows;
	using _Mybase::m_cols;
	using _Mybase::m_data;
	using _Mybase::m_shape;

	size_t m_format = rmaj|gene;

	_Myalty m_storage;

	MATRICE_GLOBAL_INL _Myt& _Reset_no_alloc(initlist<size_t> il) noexcept {
		m_shape = il;
		m_rows = m_shape.rows();
		m_cols = m_shape.cols();
		m_data = m_storage.data();
		_Mybase::_Flush_view_buf();
		return (*this);
	}

	MATRICE_GLOBAL_INL _Myt& _Reset(initlist<size_t> il) noexcept {
		_Reset_no_alloc(il);
		m_data = allocator().create(m_rows, m_cols).data();
		return (*this);
	}

public:
	MATRICE_HOST_INL _Myalty& allocator() noexcept {
		return (m_storage);
	}
	MATRICE_HOST_INL const _Myalty& allocator()const noexcept {
		return (m_storage);
	}
	MATRICE_HOST_INL decltype(auto) deleter() noexcept {
		return (m_storage.deleter());
	}

#undef MATRICE_MAKE_EXPOP_TYPE
#undef MATRICE_LINK_PTR
#undef MATRICE_EVALEXP_TOTHIS
#undef MATRICE_MAKE_ARITHOP
#undef MATRICE_MAKE_EXP_ASSIGNOP
}; _TYPES_END

_DETAIL_BEGIN
struct _Matrix_padding {
	template<typename _Ty, int _M, int _N>
	using _Matrix_t = types::Matrix_<_Ty, _M, _N>;

	/**
	 *\brief Make zero padding
	 *\param <_S> templated padding size, which should be specified for managed matrix type
	 *\param [_In] input matrix; [_Size] run-time specified padding size
	 */
	template<typename _Mty, size_t _S = 0, 
		typename _Ty = typename _Mty::value_t> 
	MATRICE_GLOBAL_INL static auto zero(const _Mty& _In, size_t _Size = _S) {
		// _Size <- max(_Size, _S)
		static_assert(is_matrix_v<_Mty>, "_Mty must be a matrix type.");
		
		constexpr auto _M = _Mty::CompileTimeRows;
		constexpr auto _N = _Mty::CompileTimeCols;
		_Matrix_t<_Ty, _M + (_S << 1), _N + (_S << 1)> _Ret;
		if constexpr (_S > 0) _Ret = { zero<_Ty> };
		else {
			_Ret.create(_In.rows()+(_Size<<1), _In.cols()+(_Size << 1), zero<_Ty>);
		}
		_Ret.block(_Size, _Size+_In.cols(), _Size, _Size+_In.rows()) = _In;

		return forward<decltype(_Ret)>(_Ret);
	}
};
_DETAIL_END
using padding =  detail::_Matrix_padding;

template<typename _Mty>
MATRICE_HOST_INL auto make_matrix_deleter(const _Mty& _M) noexcept;

/**
 *\brief dgelom::Matrix_<...> factory function to create matrix object.
 *\param [params] wrap any parameters supported by dgelom::Matrix_ ctors, while except initializer list. 
 */
template<typename _Ty, int _Rows=0, int _Cols=0, typename... _Args>
MATRICE_GLOBAL_INL types::Matrix_<_Ty, _Rows, _Cols> make_matrix(_Args&&... params);

/**
 *\brief set input data to zero
 *\param [data] can be scalar, dgelom::Matrix_, dgelom::Tensor or any container with methods data() and size().
 */
template<typename _Ty>
MATRICE_GLOBAL_INL remove_all_t<_Ty>& make_zero(_Ty& data) noexcept;

DGE_MATRICE_END

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "inl\_base.inl"