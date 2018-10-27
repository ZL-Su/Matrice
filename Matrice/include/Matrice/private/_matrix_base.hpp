/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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

#include <iterator>
#include <valarray>
#include <functional>
#include <tuple>
#include "_type_traits.h"
#include "_matrix_expr.hpp"
#include "_abstract_ops.hpp"
#include "_storage.hpp"
#include "_iterator.h"
#include "_view.h"
#include "../../../addin/interface.h"
#include "../util/_type_defs.h"
#include "../util/_macro_conditions.h"
#include "../core/solver.h"

#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)

#ifndef _HAS_CXX17
#  error Matrice library must be compiled as C++ 17 or above.
#endif

DGE_MATRICE_BEGIN
_INTERNAL_BEGIN
// \CLASS TEMPLATE to make plane view for all matrix types
template<typename _Ty, int _Type = _View_trait<_Ty>::value> 
class PlaneView_ {
	enum { MAGIC_VAL = 0x42FF0000 };
	struct Step { int buf[2] = { sizeof(_Ty), sizeof(_Ty) }; int* p = buf; };
	int type = _Type, flags = MAGIC_VAL | _Type, dims = 2; Step step;
public:
	using plvt_type = std::tuple<int, int, std::add_pointer_t<_Ty>>;
	MATRICE_GLOBAL_FINL PlaneView_() = default;
	MATRICE_GLOBAL_FINL PlaneView_(int _rows, int _cols, std::add_pointer_t<_Ty> _data = nullptr)
		: m_rows(_rows), m_cols(_cols), m_data(_data) { 
		step.buf[0] = m_cols * step.buf[1]; 
	};

	/**
	 * \shape of the matrix
	 * Example: auto [_Rows, _Cols] = _Matrix.shape();
	 */
	MATRICE_GLOBAL_FINL constexpr auto shape() const { 
		return std::tie(m_rows, m_cols); 
	}

	/**
	 * \tuple of plane view of the matrix
	 * Example: auto [_Rows, _Cols, _Data] = _Matrix.plvt();
	 */
	MATRICE_GLOBAL_FINL constexpr plvt_type plvt() {
		return std::tie(m_rows, m_cols, m_data);
	}
	MATRICE_GLOBAL_FINL constexpr plvt_type plvt() const {
		return std::tie(m_rows, m_cols, m_data);
	}

	/**
	 * \raw plane view of the matrix
	 * Example: auto& _Raw = _Matrix.raw();
	 */
	MATRICE_GLOBAL_FINL constexpr auto& raw() {
		return *this;
	}
	MATRICE_GLOBAL_FINL constexpr auto& raw() const {
		return *this;
	}

protected:
	MATRICE_GLOBAL_FINL void _Flush_view_buf() { step.buf[0] = m_cols * step.buf[1]; }

	int m_rows, m_cols;
	std::add_pointer_t<_Ty> m_data = nullptr;
};
_INTERNAL_END

_TYPES_BEGIN
template<typename _Ty> using nested_initializer_list = std::initializer_list<std::initializer_list<_Ty>>;

/*******************************************************************
	              Generic Base for Matrix Class
	    Copyright (c) : Zhilong (Dgelom) Su, since 14/Feb/2018
 ******************************************************************/
template<
	typename _Derived, 
	typename _Traits = matrix_traits<_Derived>, 
	typename _Type = typename _Traits::type>
class Base_ : public internal::PlaneView_<_Type>
{
#define _PTRLINK { m_data = m_storage.data(); }
#define _PTRLINK_EVAL { m_data = m_storage.data(); expr.assign(*this); }
#define _SHAPE std::get<0>(_Shape), std::get<1>(_Shape)
#define _EXPOP(_Desc, _Name) typename _Exp_op::_##_Desc##_##_Name<_Type>

	enum { _M = _Traits::size::rows::value, _N = _Traits::size::cols::value };
	using _Myt_storage_type = typename detail::Storage_<_Type>::template Allocator<_M, _N>;
	using _Myt = Base_;
	using _Myt_const = std::add_const_t<_Myt>;
	using _Myt_reference = std::add_lvalue_reference_t<_Myt>;
	using _Myt_const_reference = std::add_lvalue_reference_t<_Myt_const>;
	using _Myt_move_reference = std::add_rvalue_reference_t<_Myt>;
	using _Myt_fwd_iterator = _Matrix_forward_iterator<_Type>;
	using _Myt_rwise_iterator = _Matrix_rwise_iterator<_Type>;
	using _Myt_cwise_iterator = _Matrix_cwise_iterator<_Type>;
	using _Myt_rview_type = _Matrix_rview<_Type>;
	using _Myt_cview_type = _Matrix_cview<_Type>;
	using _Myt_blockview_type = _Matrix_block<_Type>;
	using _Xop_ewise_sum   = _EXPOP(Ewise, sum);
	using _Xop_ewise_min   = _EXPOP(Ewise, min);
	using _Xop_ewise_mul   = _EXPOP(Ewise, mul);
	using _Xop_ewise_div   = _EXPOP(Ewise, div);
	using _Xop_ewise_sqrt  = _EXPOP(Ewise, sqrt);
	using _Xop_ewise_exp   = _EXPOP(Ewise, exp);
	using _Xop_ewise_log   = _EXPOP(Ewise, log);
	using _Xop_ewise_log2  = _EXPOP(Ewise, log2);
	using _Xop_ewise_log10 = _EXPOP(Ewise, log10);
	using _Xop_mat_mul     = _EXPOP(Mat, mul);
	using _Xop_mat_inv     = _EXPOP(Mat, inv);
	using _Xop_mat_trp     = _EXPOP(Mat, trp);
public:
	using value_t = _Type;
	using value_type = value_t;
	using matrix_type = _Derived;
	using pointer = std::add_pointer_t<value_t>;
	using reference = std::add_lvalue_reference_t<value_t>;
	using iterator = pointer;
	using const_iterator = std::add_const_t<iterator>;
	using const_init_list = std::add_const_t<std::initializer_list<value_t>>;
	using shape_t = std::tuple<std::size_t, std::size_t>;
	using base_t = internal::PlaneView_<value_t>;
	using loctn_t = Location;

	template<typename _Xop> using expr_type = Expr::Base_<_Xop>;
	
	enum { options = _Myt_storage_type::location };
	enum { Size = _M*_N, CompileTimeRows = _M, CompileTimeCols = _N, };
	static constexpr auto inf = std::numeric_limits<value_t>::infinity();
	static constexpr auto eps = std::numeric_limits<value_t>::epsilon();

	MATRICE_GLOBAL_FINL Base_() noexcept
		:base_t(_M<0?0:_M, _N<0?0:_N), m_storage() _PTRLINK
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols) noexcept 
		:base_t(_rows, _cols), m_storage(_rows, _cols) _PTRLINK
	MATRICE_GLOBAL_FINL Base_(const shape_t& _Shape) noexcept
		:base_t(_SHAPE), m_storage(m_rows, m_cols) _PTRLINK
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols, pointer data) noexcept 
		:base_t(_rows, _cols, data), m_storage(_rows, _cols, data) {}
	MATRICE_GLOBAL_FINL Base_(const shape_t& _Shape, pointer _Data) noexcept
		:base_t(_SHAPE), m_storage(m_rows, m_cols, _Data) _PTRLINK
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols, value_t _val) noexcept 
		:base_t(_rows, _cols), m_storage(_rows, _cols, _val) _PTRLINK
	MATRICE_GLOBAL_FINL Base_(const shape_t& _Shape, value_t _Val) noexcept
		:base_t(_SHAPE), m_storage(m_rows, m_cols, _Val) _PTRLINK
	MATRICE_GLOBAL_FINL Base_(const_init_list _list) noexcept 
		:base_t(_M, _N), m_storage(_list) _PTRLINK
	MATRICE_GLOBAL_FINL Base_(_Myt_const_reference _other) noexcept 
		:base_t(_other.m_rows, _other.m_cols), m_storage(_other.m_storage) _PTRLINK
	MATRICE_GLOBAL_FINL Base_(_Myt_move_reference _other) noexcept 
		:base_t(_other.m_rows, _other.m_cols), m_storage(std::move(_other.m_storage)) _PTRLINK
	/**
	 *\from STD valarray<...>
	 */
	MATRICE_GLOBAL_FINL Base_(const std::valarray<value_t>& _other, int _rows = 1) noexcept 
		:base_t(_rows, _other.size() / _rows), m_storage(_rows, _other.size() / _rows, _other.data()) _PTRLINK
	/**
	 *\from explicit specified matrix type
	 */
	template<int _Rows, int _Cols, typename _Mty = Matrix_<value_t, _Rows, _Cols>>
	MATRICE_GLOBAL_FINL Base_(const _Mty& _other) noexcept 
		:base_t(_other.rows(), _other.cols()), m_storage(_other.m_storage) _PTRLINK
	/**
	 *\from expression
	 */
	template<typename _Exp>
	MATRICE_GLOBAL_FINL Base_(const exprs::_Matrix_exp<_Exp>& expr)
		: base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) _PTRLINK_EVAL
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Op>& expr)
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) _PTRLINK_EVAL
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::MatBinaryExpr<_Lhs, _Rhs, _Op>& expr) 
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) _PTRLINK_EVAL
	template<typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::MatUnaryExpr<_Rhs, _Op>& expr) 
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) _PTRLINK_EVAL

	///<brief> opencv interface </brief>
#ifdef __use_ocv_as_view__
	MATRICE_HOST_FINL Base_(const ocv_view_t& mat) : base_t(mat.rows, mat.cols, mat.ptr<value_t>()) {}
	MATRICE_HOST_FINL operator ocv_view_t() { return ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data); }
	MATRICE_HOST_FINL operator ocv_view_t() const { return ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data); }
#endif
public:
	/**
	 *\dynamic matrix creation
	 */
	MATRICE_HOST_ONLY void create(std::size_t rows, std::size_t cols) { 
		if constexpr (_M <= 0 && _N <= 0) static_cast<_Derived*>(this)->create(rows, cols); 
	};
	MATRICE_HOST_ONLY void create(const shape_t& _Shape) {
		if constexpr (_M <= 0 && _N <= 0) static_cast<_Derived*>(this)->create(_SHAPE);
	};

	///<brief> accessors </brief>
	MATRICE_GLOBAL_FINL pointer operator[](index_t y) { return (m_data + y * m_cols); }
	MATRICE_GLOBAL_FINL const pointer operator[](index_t y) const { return (m_data + y * m_cols); }
	MATRICE_GLOBAL_FINL reference operator()(index_t i) { return m_data[i]; }
	MATRICE_GLOBAL_FINL const reference operator()(index_t i) const { return m_data[i]; }
	template<size_t _Format = rmaj>
	MATRICE_GLOBAL_INL reference operator()(int _r, int _c) const { return (m_data + (_Format == rmaj ? m_cols : m_rows) * _r)[_c]; }

	///<brief> access methods </brief>
	MATRICE_GLOBAL_FINL pointer data() { return (m_data); }
	MATRICE_GLOBAL_FINL const pointer data() const { return (m_data); }
	MATRICE_GLOBAL_FINL pointer ptr(int y = 0) { return (m_data + m_cols * y); }
	MATRICE_GLOBAL_FINL const pointer ptr(int y = 0) const { return (m_data + m_cols * y); }

	MATRICE_GLOBAL_FINL constexpr iterator begin() { return (m_data); }
	MATRICE_GLOBAL_FINL constexpr iterator end() { return (m_data + size()); }
	MATRICE_GLOBAL_FINL constexpr const_iterator begin() const { return (m_data); }
	MATRICE_GLOBAL_FINL constexpr const_iterator end() const { return (m_data + size()); }

	/**
	 * \eval() expression, return this reference for the true matrix
	 */
	MATRICE_GLOBAL_FINL constexpr auto& eval() const { return (*static_cast<const _Derived*>(this)); }

#pragma region <!-- iterators -->
	/**
	 *\column iterator for accessing elements in i-th column
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
	MATRICE_GLOBAL_INL _Myt_blockview_type block(index_t x0, index_t x1, index_t y0, index_t y1) {
		return _Myt_blockview_type(m_data, m_cols, {x0, y0, x1, y1});
	}
	MATRICE_GLOBAL_INL const _Myt_blockview_type block(index_t x0, index_t x1, index_t y0, index_t y1) const {
		return _Myt_blockview_type(m_data, m_cols, { x0, y0, x1, y1 });
	}
	/**
	 * \View of the matrix data
	 */

#pragma endregion

#ifdef _HAS_CXX17
	/**
	 * \interface for STD valarray<...>
	 * Example: std::valarray<float> _Valarr = _M;
	 */
	MATRICE_GLOBAL_FINL operator std::valarray<value_t>() { return std::valarray<value_t>(m_data, size()); }
	MATRICE_GLOBAL_FINL operator std::valarray<value_t>() const { return std::valarray<value_t>(m_data, size()); }
#endif

	/**
	 * \size properties
	 */
	MATRICE_GLOBAL_FINL constexpr auto rows() const { return m_rows; }
	MATRICE_GLOBAL_FINL constexpr auto cols() const { return m_cols; }
	MATRICE_GLOBAL_FINL constexpr auto size() const { return m_rows*m_cols; }

	///<brief> assignment operators </brief>

	/**
	 * \assignment operator, from initializer list
	 */
	MATRICE_GLOBAL_FINL _Derived& operator= (const_init_list _list) {
#ifdef _DEBUG
		assert(size() == m_storage.size());
#endif
		m_storage = _list;
		m_data = internal::_Proxy_checked(m_storage.data());
		return (*static_cast<_Derived*>(this));
	}
	/**
	 * \assignment operator, from row-wise iterator
	 */
	MATRICE_GLOBAL_FINL _Derived& operator= (const _Myt_rwise_iterator& _It) {
		std::copy(_It.begin(), _It.end(), m_storage.data());
		m_data = internal::_Proxy_checked(m_storage.data());
		return (*static_cast<_Derived*>(this));
	}
	/**
	 * \assignment operator, from column-wise iterator
	 */
	MATRICE_GLOBAL_FINL _Derived& operator= (const _Myt_cwise_iterator& _It) {
		std::copy(_It.begin(), _It.end(), m_storage.data());
		m_data = internal::_Proxy_checked(m_storage.data());
		return (*static_cast<_Derived*>(this));
	}
	/**
	 * \homotype copy assignment operator
	 */
	MATRICE_GLOBAL_INL _Derived& operator= (_Myt_const_reference _other) {
		m_cols = _other.m_cols, m_rows = _other.m_rows;
		m_format = _other.m_format;
		m_storage = _other.m_storage;
		m_data = internal::_Proxy_checked(m_storage.data());
		base_t::_Flush_view_buf();
		return (*static_cast<_Derived*>(this));
	}
	/**
	 * \homotype move assignment operator
	 */
	MATRICE_GLOBAL_INL _Derived& operator= (_Myt_move_reference _other) {
		m_cols = _other.m_cols, m_rows = _other.m_rows;
		m_format = _other.m_format;
		m_storage = std::forward<_Myt_storage_type>(_other.m_storage);
		m_data = internal::_Proxy_checked(m_storage.data());
		base_t::_Flush_view_buf();
		_other.m_storage.free();
		_other.m_data = nullptr;
		_other.m_cols = _other.m_rows = 0;
		return (*static_cast<_Derived*>(this));
	}
	/**
	 * \try to convert managed matrix to dynamic matrix
	 */
	template<int _Rows, int _Cols,
		typename _Mty = Matrix_<value_t, _Rows, _Cols>,
		typename = typename std::enable_if_t<is_static_v<_Rows, _Cols>&&!is_static_v<_M, _N>>>
	MATRICE_GLOBAL_FINL _Derived& operator= (_Mty& _managed) {
		m_data = _managed.data();
		m_rows = _managed.rows(), m_cols = _managed.cols();
		m_storage.owner() = detail::Storage_<value_t>::Proxy;
		base_t::_Flush_view_buf();
		return (*static_cast<_Derived*>(this));
	}

#pragma region <!-- Lazied Operators for Matrix Arithmetic -->
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator+ (const _Rhs& _opd) const { return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_sum>(*this, _opd);}
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator- (const _Rhs& _opd) const { return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_min>(*this, _opd); }
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator* (const _Rhs& _opd) const { return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_mul>(*this, _opd); }
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator/ (const _Rhs& _opd) const { return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_div>(*this, _opd); }

	template<typename _Rhs> MATRICE_GLOBAL_INL auto mul(const _Rhs& rhs) const { return Expr::MatBinaryExpr<_Myt, _Rhs, _Xop_mat_mul>(*this, rhs); }
	MATRICE_GLOBAL_FINL auto sqrt() const { return Expr::EwiseUnaryExpr<_Myt, _Xop_ewise_sqrt>(*this); }
	MATRICE_HOST_FINL auto inv() const { return Expr::MatUnaryExpr<_Myt, _Xop_mat_inv>(*this); }
	MATRICE_HOST_FINL auto inv(_Myt_const_reference _rhs) { return Expr::MatUnaryExpr<_Myt, _Xop_mat_inv>(_rhs, *this); }
	MATRICE_HOST_FINL auto transpose() { return Expr::MatUnaryExpr<_Myt, _Xop_mat_trp>(*this); }
	MATRICE_GLOBAL_FINL auto normalize(value_t _val = inf) { return ((abs(_val) < eps ? value_t(1) : value_t(1) / (_val == inf ? max() : _val))*(*this)); }
	MATRICE_GLOBAL_INL Expr::MatBinaryExpr<_Myt, _Myt, _Xop_mat_mul> spread();
#pragma endregion

#pragma region <!-- Triggers for Suspended Expression -->
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const Expr::EwiseBinaryExpr<_Args...>& _Ex){ return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const Expr::EwiseUnaryExpr<_Args...>& _Ex) { return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const Expr::MatBinaryExpr<_Args...>& _Ex) { return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const Expr::MatUnaryExpr<_Args...>& _Ex) { return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
#pragma endregion

	///<brief> in-time matrix arithmetic </brief>
	MATRICE_GLOBAL_FINL auto max() const { return (*std::max_element(begin(), end())); }
	MATRICE_GLOBAL_FINL auto min() const { return (*std::min_element(begin(), end())); }
	MATRICE_GLOBAL_FINL auto sum() const { return (reduce(begin(), end())); }
	MATRICE_GLOBAL_FINL auto det() const { return (det_impl(*static_cast<const _Derived*>(this))); }
	MATRICE_GLOBAL_FINL auto trace() const { return (reduce(begin(), end(), cols() + 1)); }

	/**
	 * \matrix Frobenius norm
	 */
	MATRICE_GLOBAL_FINL auto norm_2() const { auto _Ans = dot(*this); return (_Ans > eps ? dgelom::sqrt(_Ans) : inf); }
	/**
	 * \matrix p-norm: $[\sum_{i=1}^{m}\sum_{j=1}^{}|a_{ij}|^p]^{1/p}$
	 * Sepcial cases: $p = 0$ for $\infty$-norm, $p = 1$ for 1-norm and $p = 2$ for 2-norm
	 * Reference: https://en.wikipedia.org/wiki/Matrix_norm
	 */
	template<std::size_t _P = 2> MATRICE_GLOBAL_FINL auto norm() const {
		return internal::_Matrix_norm_impl<_P>::value(*(this));
	}
	/**
	 * \dot product of this matrix with _Rhs
	 */
	template<typename _Rhs, typename = std::enable_if_t<is_matrix_v<_Rhs>>> 
	MATRICE_GLOBAL_FINL auto dot(const _Rhs& _Rhs) const { return (this->operator*(_Rhs)).sum(); }

	/**
	 * \operate each entry via _Fn
	 */
	template<typename _Fty>
	MATRICE_GLOBAL_FINL void each(_Fty _Fn) { for (auto& _Val : *this) _Fn(_Val); }
	/**
	 * \replace entries meets _Cond with _Val
	 */
	MATRICE_GLOBAL_FINL void where(std::function<bool(const value_type&)> _Cond, const value_type _Val) {
		this->each([&](auto& _My_val) {_My_val = _Cond(_My_val) ? _Val : _My_val; });
	}

	///<brief> properties </brief>
	__declspec(property(get=_Format_getter, put=_Format_setter))size_t format;
	MATRICE_HOST_FINL size_t _Format_getter() const { return m_format; }
	MATRICE_HOST_FINL void _Format_setter(size_t format) { m_format = rmaj|format; }

	__declspec(property(get = _Empty_getter)) bool empty;
	MATRICE_HOST_FINL bool _Empty_getter() const { return (size() == 0); }

protected:
	using base_t::m_rows;
	using base_t::m_cols;
	using base_t::m_data;

	size_t m_format = rmaj|gene;

public:
	_Myt_storage_type m_storage;

#undef _SHAPE
#undef _EXPOP
#undef _PTRLINK
#undef _PTRLINK_EVAL
};
_TYPES_END
DGE_MATRICE_END