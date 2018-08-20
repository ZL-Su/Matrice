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
#include "_type_traits.h"
#include "_matrix_expr.hpp"
#include "_matrix.inl.hpp"
#include "_storage.hpp"
#include "_iterator.h"
#include "_view.h"
#include "../../../addin/interface.h"
#include "../util/_type_defs.h"
#include "../core/solver.h"

#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)

#ifndef _HAS_CXX17
#  error Matrice library must be compiled as C++ 17
#endif

MATRICE_NAMESPACE_BEGIN_

#ifndef MATRICE_INT
#define MATRICE_INT                                std::ptrdiff_t
	typedef MATRICE_INT                                     int_t;
#else
	typedef int                                             int_t;
#endif

template<typename _Ty, int _Type = _View_trait<_Ty>::value> class PlaneView_
{
	enum { MAGIC_VAL = 0x42FF0000 };
	struct Step { int buf[2] = { sizeof(_Ty), sizeof(_Ty) }; int* p = buf; };
	int type = _Type, flags = MAGIC_VAL | _Type, dims = 2; Step step;
public:
	constexpr MATRICE_GLOBAL_FINL PlaneView_() = default;
	constexpr MATRICE_GLOBAL_FINL PlaneView_(int _rows, int _cols, _Ty* _data = nullptr)
		: m_rows(_rows), m_cols(_cols), m_data(_data) { step.buf[0] = m_cols * step.buf[1]; };
	constexpr MATRICE_GLOBAL_FINL void update(_Ty* data) { m_data = data; }
protected:
	constexpr MATRICE_GLOBAL_FINL void _Flush_view_buf() { step.buf[0] = m_cols * step.buf[1]; }
protected:
	int m_rows, m_cols;
	_Ty* m_data = nullptr;
};

namespace types {

template<typename _Ty> using nested_initializer_list = std::initializer_list<std::initializer_list<_Ty>>;
template<typename _Ty, int _M, int _N> class Matrix_;
template<typename _Ty, int _M, int _N> struct matrix_traits<Matrix_<_Ty, _M, _N>> {
	using type = _Ty;
	enum { M = _M, N = _N };
	struct size { struct rows { enum { value = M }; }; struct cols { enum { value = N }; }; };
};
template</*typename _Ty, int _M, int _N, */typename _Derived> class Base_;

/*******************************************************************
	              Generic Base for Matrix Class
	    Copyright (c) : Zhilong (Dgelom) Su, since 14/Feb/2018
 ******************************************************************/
//template<typename _Ty, int _M, int _N>
template<typename _Derived> 
class Base_ : public PlaneView_<typename matrix_traits<_Derived>::type>
{
	using _Myt_traits = matrix_traits<_Derived>;
	enum { _M = _Myt_traits::size::rows::value, _N = _Myt_traits::size::cols::value };
	using _Myt_storage_type = typename details::Storage_<typename _Myt_traits::type>::
			Allocator<_Myt_traits::size::rows::value, _Myt_traits::size::cols::value>;
	using _Myt = Base_;
	using _Myt_const = std::add_const_t<_Myt>;
	using _Myt_reference = std::add_lvalue_reference_t<_Myt>;
	using _Myt_const_reference = std::add_lvalue_reference_t<_Myt_const>;
	using _Myt_move_reference = std::add_rvalue_reference_t<Base_>;
public:
	using value_t = typename _Myt_traits::type;
	using value_type = value_t;
	using pointer = std::add_pointer_t<value_t>;
	using reference = std::add_lvalue_reference_t<value_t>;
	using iterator = pointer;
	using const_iterator = std::add_const_t<iterator>;
	using const_init_list = std::add_const_t<std::initializer_list<value_t>>;
	template<typename _Xop> using expr_type = Expr::Base_<_Xop>;
	using base_t = PlaneView_<value_t>;
	using loctn_t = Location;
	using _Myt_fwd_iterator = _Matrix_forward_iterator<value_t>;
	using _Myt_rwise_iterator = _Matrix_rwise_iterator<value_t>;
	using _Myt_cwise_iterator = _Matrix_cwise_iterator<value_t>;
	using _Myt_rview_type = _Matrix_rview<value_t>;
	using _Myt_cview_type = _Matrix_cview<value_t>;
	using _Myt_blockview_type = _Matrix_block<value_t>;
	using _Xop_ewise_sum = typename Expr::Op::EwiseSum<value_t>;
	using _Xop_ewise_min = typename Expr::Op::EwiseMin<value_t>;
	using _Xop_ewise_mul = typename Expr::Op::EwiseMul<value_t>;
	using _Xop_ewise_div = typename Expr::Op::EwiseDiv<value_t>;
	using _Xop_mat_mul = typename Expr::Op::MatMul<value_t>;
	using _Xop_mat_inv = typename Expr::Op::MatInv<value_t>;
	using _Xop_mat_trp = typename Expr::Op::MatTrp<value_t>;

	enum { options = _Myt_storage_type::location };
	constexpr static const value_t inf = std::numeric_limits<value_t>::infinity();
	constexpr static const value_t eps = std::numeric_limits<value_t>::epsilon();

	MATRICE_GLOBAL_FINL Base_() noexcept :base_t(_M<0?0:_M, _N<0?0:_N), m_storage() { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols) noexcept :base_t(_rows, _cols), m_storage(_rows, _cols) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols, pointer data) noexcept :base_t(_rows, _cols, data), m_storage(_rows, _cols, data) {}
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols, value_t _val) noexcept :base_t(_rows, _cols), m_storage(_rows, _cols, _val) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(const_init_list _list) noexcept :base_t(_M, _N), m_storage(_list) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(_Myt_const_reference _other) noexcept :base_t(_other.m_rows, _other.m_cols), m_storage(_other.m_storage) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(_Myt_move_reference _other) noexcept :base_t(_other.m_rows, _other.m_cols), m_storage(std::move(_other.m_storage)) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(const std::valarray<value_t>& _other, int _rows = 1) noexcept :base_t(_rows, _other.size() / _rows), m_storage(_rows, _other.size() / _rows, _other.data()) { m_data = m_storage.data(); }
	template<int _Rows, int _Cols>
	MATRICE_GLOBAL_FINL Base_(const Matrix_<value_t, _Rows, _Cols>& _other) noexcept :base_t(_other.rows(), _other.cols()), m_storage(_other.m_storage) { m_data = m_storage.data(); }
	
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Op>& expr)
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) { m_data = m_storage.data(); expr.assign(*this); }
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::MatBinaryExpr<_Lhs, _Rhs, _Op>& expr) 
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) { m_data = m_storage.data(); expr.assign(*this); }
	template<typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const Expr::MatUnaryExpr<_Rhs, _Op>& expr) 
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) { m_data = m_storage.data(); expr.assign(*this); }

	///<brief> opencv interface </brief>
#ifdef __use_ocv_as_view__
	MATRICE_HOST_FINL Base_(const ocv_view_t& mat) : base_t(mat.rows, mat.cols, mat.ptr<value_t>()) {}
	MATRICE_HOST_FINL operator ocv_view_t() { return ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data); }
	MATRICE_HOST_FINL operator ocv_view_t() const { return ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data); }
#endif
public:
	///<brief> dynamic create methods </brief>
	MATRICE_HOST_ONLY void create(size_t rows, size_t cols) { if constexpr (_M == 0) static_cast<_Derived*>(this)->create(rows, cols); };
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
	MATRICE_GLOBAL_FINL constexpr const_iterator begin() const { return (m_data); }
	MATRICE_GLOBAL_FINL constexpr const_iterator end() const { return (m_data + size()); }
	MATRICE_GLOBAL_FINL constexpr iterator begin() { return (m_data); }
	MATRICE_GLOBAL_FINL constexpr iterator end() { return (m_data + size()); }
#pragma region <!-- iterators -->
	//column iterator for accessing elements in current column
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

	//row iterator for accessing elements in current row
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

	//column-wise iterator for accessing elements
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

	//row-wise iterator for accessing elements
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
	//view of i-th row 
	MATRICE_GLOBAL_FINL _Myt_rview_type rview(size_t i) {
		return _Myt_rview_type(m_data + m_cols * i, m_cols);
	}
	MATRICE_GLOBAL_FINL const _Myt_rview_type rview(size_t i) const {
		return _Myt_rview_type(m_data + m_cols * i, m_cols);
	}
	//view of i-th column
	MATRICE_GLOBAL_FINL _Myt_cview_type cview(size_t i) {
		return _Myt_cview_type(m_data + i, m_rows, m_cols, i);
	}
	MATRICE_GLOBAL_FINL const _Myt_cview_type cview(size_t i) const {
		return _Myt_cview_type(m_data + i, m_rows, m_cols, i);
	}
	//view of submatrix [x0, x1) : [y0, y1)
	MATRICE_GLOBAL_INL _Myt_blockview_type block(index_t x0, index_t x1, index_t y0, index_t y1) {
		return _Myt_blockview_type(m_data, m_cols, {x0, y0, x1, y1});
	}
	MATRICE_GLOBAL_INL const _Myt_blockview_type block(index_t x0, index_t x1, index_t y0, index_t y1) const {
		return _Myt_blockview_type(m_data, m_cols, { x0, y0, x1, y1 });
	}
#pragma endregion
#ifdef _HAS_CXX17
	MATRICE_GLOBAL_FINL operator std::valarray<value_t>() { return std::valarray<value_t>(m_data, size()); }
	MATRICE_GLOBAL_FINL operator std::valarray<value_t>() const { return std::valarray<value_t>(m_data, size()); }
#endif
	MATRICE_GLOBAL_FINL constexpr auto rows() const { return m_rows; }
	MATRICE_GLOBAL_FINL constexpr auto cols() const { return m_cols; }
	MATRICE_GLOBAL_FINL constexpr auto size() const { return m_rows*m_cols; }


	///<brief> assignment operators </brief>
	MATRICE_GLOBAL_FINL _Derived& operator= (const_init_list _list)
	{
#ifdef _DEBUG
		assert(size() == m_storage.size());
#endif
		m_storage = _list;
		m_data = _Proxy_checked(m_storage.data());
		return (*static_cast<_Derived*>(this));
	}
	MATRICE_GLOBAL_INL _Derived& operator= (_Myt_const_reference _other) //homotype assignment operator
	{
		m_cols = _other.m_cols, m_rows = _other.m_rows;
		m_format = _other.m_format;
		m_storage = _other.m_storage;
		m_data = _Proxy_checked(m_storage.data());
		base_t::_Flush_view_buf();
		return (*static_cast<_Derived*>(this));
	}
	MATRICE_GLOBAL_INL _Derived& operator= (_Myt_move_reference _other) //homotype assignment operator
	{
		m_cols = _other.m_cols, m_rows = _other.m_rows;
		m_format = _other.m_format;
		m_storage = std::forward<_Myt_storage_type>(_other.m_storage);
		m_data = _Proxy_checked(m_storage.data());
		base_t::_Flush_view_buf();
		_other.m_storage.free();
		_other.m_data = nullptr;
		_other.m_cols = _other.m_rows = 0;
		return (*static_cast<_Derived*>(this));
	}
	template<int _Rows, int _Cols, typename = typename std::enable_if<is_static<_Rows, _Cols>::value&&!is_static<_M, _N>::value>::type> //static to dynamic
	MATRICE_GLOBAL_FINL _Derived& operator= (Matrix_<value_t, _Rows, _Cols>& _managed)
	{
		m_data = _managed.data();
		m_rows = _managed.rows(), m_cols = _managed.cols();
		m_storage.owner() = details::Storage_<value_t>::Proxy;
		base_t::_Flush_view_buf();
		return (*static_cast<_Derived*>(this));
	}

#pragma region <!-- Lazied Operators for Matrix Arithmetic -->
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator+ (const _Rhs& _opd) const 
	{ return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_sum>(*this, _opd);}
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator- (const _Rhs& _opd) const 
	{ return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_min>(*this, _opd); }
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator* (const _Rhs& _opd) const 
	{ return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_mul>(*this, _opd); }
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator/ (const _Rhs& _opd) const 
	{ return Expr::EwiseBinaryExpr<_Myt, _Rhs, _Xop_ewise_div>(*this, _opd); }
	template<typename _Rhs> MATRICE_GLOBAL_INL auto mul(const _Rhs& rhs) const 
	{ return Expr::MatBinaryExpr<_Myt, _Rhs, _Xop_mat_mul>(*this, rhs); }
	MATRICE_HOST_FINL auto inv() const 
	{ return Expr::MatUnaryExpr<_Myt, _Xop_mat_inv>(*this); }
	MATRICE_HOST_FINL auto inv(_Myt_const_reference _rhs) 
	{ return Expr::MatUnaryExpr<_Myt, _Xop_mat_inv>(_rhs, *this); }
	MATRICE_HOST_FINL auto transpose() 
	{ return Expr::MatUnaryExpr<_Myt, _Xop_mat_trp>(*this); }
	MATRICE_GLOBAL_FINL auto normalize(value_t _val = inf) 
	{ return ((abs(_val) < eps ? value_t(1) : value_t(1) / (_val == inf ? max() : _val))*(*this)); }
	MATRICE_GLOBAL_INL Expr::MatBinaryExpr<_Myt, _Myt, _Xop_mat_mul> spread();
#pragma endregion

#pragma region <!-- Triggers for Suspended Expression -->
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const Expr::EwiseBinaryExpr<_Args...>& _Ex){ return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const Expr::MatBinaryExpr<_Args...>& _Ex) { return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const Expr::MatUnaryExpr<_Args...>& _Ex) { return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
#pragma endregion

	///<brief> in-time matrix arithmetic </brief>
	MATRICE_GLOBAL_FINL value_t max() const { return (*std::max_element(m_data, m_data + size())); }
	MATRICE_GLOBAL_FINL value_t min() const { return (*std::min_element(m_data, m_data + size())); }
	MATRICE_GLOBAL_FINL value_t sum() const { return (reduce<value_t>(begin(), end())); }
	MATRICE_GLOBAL_FINL value_t det() const { return (det_impl(*static_cast<const _Derived*>(this))); }
	MATRICE_GLOBAL_FINL value_t trace() const { return (reduce(begin(), end(), cols() + 1)); }
	template<class _Rhs> MATRICE_GLOBAL_FINL value_t dot(const _Rhs& _rhs) const { return (this->operator*(_rhs)).sum(); }
	MATRICE_GLOBAL_FINL value_t norm_2() const { auto ans = dot(*this); return (ans > eps ? sqrt(ans) : inf); }

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
};
MATRICE_NAMESPACE_END_TYPES