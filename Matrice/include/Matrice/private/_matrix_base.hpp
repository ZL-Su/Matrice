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
#include "_expr_type_traits.h"
#include "_matrix_expr.hpp"
#include "_matrix.inl.hpp"
#include "_storage.hpp"
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

template<typename T> struct _View_traits { enum { value = 0x0000 }; };
template<> struct _View_traits<unsigned char> { enum { value = 0x0008 }; };
template<> struct _View_traits<int> { enum { value = 0x0016 }; };
template<> struct _View_traits<float> { enum { value = 0x0032 }; };
template<> struct _View_traits<double> { enum { value = 0x0064 }; };
template<typename _Ty, int _Type = _View_traits<_Ty>::value> class PlaneView_
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
template<typename _Ty, int _M, int _N, typename _Derived> class Base_;

/*******************************************************************
	              Generic Base for Matrix Class
	          Copyright (c) : Zhilong Su 14/Feb/2018
 ******************************************************************/
template<typename _Ty, int _M, int _N, typename _Derived = Matrix_<_Ty, _M, _N>> 
class Base_ : public PlaneView_<_Ty>
{
	typedef details::Storage_<_Ty>::Allocator<_M, _N, 
		              allocator_option<_M, _N>::value>             Storage;
	typedef Base_                                                    _Myt;
	typedef _Myt&                                                 myt_ref;
	typedef const _Myt&                                     const_myt_ref;
	typedef _Myt&&                                               myt_move;
	typedef typename exprs::Expr::Op::EwiseSum<_Ty>            EwiseSumOp;
	typedef typename exprs::Expr::Op::EwiseMin<_Ty>            EwiseMinOp;
	typedef typename exprs::Expr::Op::EwiseMul<_Ty>            EwiseMulOp;
	typedef typename exprs::Expr::Op::EwiseDiv<_Ty>            EwiseDivOp;
	typedef typename exprs::Expr::Op::MatMul<_Ty>                MatMulOp;
	typedef typename exprs::Expr::Op::MatInv<_Ty>                MatInvOp;
	typedef typename exprs::Expr::Op::MatTrp<_Ty>                MatTrpOp;
public:
	typedef PlaneView_<_Ty>                                        base_t;
	typedef typename details::Storage_<_Ty>::value_t              value_t;
	typedef typename details::Storage_<_Ty>::pointer              pointer;
	typedef typename details::Storage_<_Ty>::reference          reference;
	typedef typename details::Storage_<_Ty>::idx_t                  idx_t;
	typedef typename Location                                     loctn_t;
	typedef typename pointer                                     iterator;
	typedef const iterator                                 const_iterator;
	typedef const std::initializer_list<value_t>          const_init_list;
	typedef enum StorageFormat { RowMajor = 101, ColMajor = 102 }format_t;
	typedef enum StorageType { GDs, GTd, GBd, GSp, DSm, BSm, SSm } type_t;
	enum { options = Expr::OpFlag::undef | Storage::location };
	typedef struct Column_View
	{
		Column_View(const _Derived& ref, int_t _c) : myref(ref), c(_c) {};
		int_t c = 0;
		const _Derived& myref;
		reference operator() (int_t r) { return myref(r, c); }
		reference operator() (int_t r) const { return myref(r, c); }
	} cview_t;
	constexpr static const _Ty inf = std::numeric_limits<_Ty>::infinity();
	constexpr static const _Ty eps = std::numeric_limits<_Ty>:: epsilon();

	MATRICE_GLOBAL_FINL Base_() noexcept :base_t(_M<0?0:_M, _N<0?0:_N), m_storage() { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols) noexcept :base_t(_rows, _cols), m_storage(_rows, _cols) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols, pointer data) noexcept :base_t(_rows, _cols, data), m_storage(_rows, _cols, data) {}
	MATRICE_GLOBAL_FINL Base_(int _rows, int _cols, value_t _val) noexcept :base_t(_rows, _cols), m_storage(_rows, _cols, _val) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(const std::initializer_list<value_t> _list) noexcept :base_t(_M, _N), m_storage(_list) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(const_myt_ref _other) noexcept :base_t(_other.m_rows, _other.m_cols), m_storage(_other.m_storage) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(myt_move _other) noexcept :base_t(_other.m_rows, _other.m_cols), m_storage(std::move(_other.m_storage)) { m_data = m_storage.data(); }
	MATRICE_GLOBAL_FINL Base_(const std::valarray<value_t>& _other, int _rows = 1) noexcept :base_t(_rows, _other.size() / _rows), m_storage(_rows, _other.size() / _rows, _other.data()) { m_data = m_storage.data(); }
	template<int _Rows, int _Cols>
	MATRICE_GLOBAL_FINL Base_(const Matrix_<value_t, _Rows, _Cols>& _other) noexcept :base_t(_other.rows(), _other.cols()), m_storage(_other.m_storage) { m_data = m_storage.data(); }
	
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const exprs::Expr::EwiseBinaryExpr<_Lhs, _Rhs, _Op>& expr)
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) { m_data = m_storage.data(); expr.assign(*this); }
	template<typename _Lhs, typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const exprs::Expr::MatBinaryExpr<_Lhs, _Rhs, _Op>& expr) 
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) { m_data = m_storage.data(); expr.assign(*this); }
	template<typename _Rhs, typename _Op>
	MATRICE_GLOBAL_FINL Base_(const exprs::Expr::MatUnaryExpr<_Rhs, _Op>& expr) 
		:base_t(m_rows = expr.rows(), m_cols = expr.cols()), m_storage(m_rows, m_cols) { m_data = m_storage.data(); expr.assign(*this); }

	///<brief> opencv interface </brief>
#ifdef __use_ocv_as_view__
	MATRICE_HOST_FINL Base_(const ocv_view_t& mat) : base_t(mat.rows, mat.cols, mat.ptr<value_t>()) {}
	MATRICE_HOST_FINL operator ocv_view_t() { return ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data); }
	MATRICE_HOST_FINL operator ocv_view_t() const { return ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data); }
#endif
public:
	///<brief> dynamic create methods </brief>
	MATRICE_HOST_ONLY void create(int_t rows, int_t cols) { if constexpr (_M == 0) static_cast<_Derived*>(this)->create(rows, cols); };
	///<brief> accessors </brief>
	MATRICE_GLOBAL_FINL pointer   operator[](int_t y) { return (m_data + y * m_cols); }
	MATRICE_GLOBAL_FINL const pointer operator[](int_t y) const { return (m_data + y * m_cols); }
	MATRICE_GLOBAL_FINL reference operator()(int_t i) { return m_data[i]; }
	MATRICE_GLOBAL_FINL const reference operator()(int_t i) const { return m_data[i]; }
	template<format_t _Fmt = RowMajor>
	MATRICE_GLOBAL_INL reference operator()(int _r, int _c) const { return (m_data + (_Fmt == RowMajor ? m_cols : m_rows) * _r)[_c]; }
	///<brief> access methods </brief>
	MATRICE_GLOBAL_FINL pointer data() { return (m_data); }
	MATRICE_GLOBAL_FINL const pointer data() const { return (m_data); }
	MATRICE_GLOBAL_FINL pointer ptr(int y = 0) { return (m_data + m_cols * y); }
	MATRICE_GLOBAL_FINL const pointer ptr(int y = 0) const { return (m_data + m_cols * y); }
	MATRICE_GLOBAL_FINL constexpr const_iterator begin() const { return (m_data); }
	MATRICE_GLOBAL_FINL constexpr const_iterator end() const { return (m_data + size()); }
	MATRICE_GLOBAL_FINL constexpr iterator begin() { return (m_data); }
	MATRICE_GLOBAL_FINL constexpr iterator end() { return (m_data + size()); }
	MATRICE_GLOBAL_FINL constexpr int_t rows() const { return m_rows; }
	MATRICE_GLOBAL_FINL constexpr int_t cols() const { return m_cols; }
	MATRICE_GLOBAL_FINL constexpr std::size_t size() const { return m_rows*m_cols; }
	MATRICE_GLOBAL_FINL constexpr cview_t col(int_t c) const { return cview_t(*this, c); }
	MATRICE_GLOBAL_FINL operator std::valarray<value_t>() { return std::valarray<value_t>(m_data, size()); }
	MATRICE_GLOBAL_FINL operator std::valarray<value_t>() const { return std::valarray<value_t>(m_data, size()); }
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
	MATRICE_GLOBAL_FINL _Derived& operator= (const_myt_ref _other) //homotype assignment operator
	{
		m_cols = _other.m_cols, m_rows = _other.m_rows;
		m_format = _other.m_format;
		m_storage = _other.m_storage;
		m_data = _Proxy_checked(m_storage.data());
		base_t::_Flush_view_buf();
		return (*static_cast<_Derived*>(this));
	}
	MATRICE_GLOBAL_FINL _Derived& operator= (myt_move _other) //homotype assignment operator
	{
		m_cols = _other.m_cols, m_rows = _other.m_rows;
		m_format = _other.m_format;
		m_storage = std::forward<Storage>(_other.m_storage);
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
	{ return exprs::Expr::EwiseBinaryExpr<_Myt, _Rhs, EwiseSumOp>(*this, _opd);}
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator- (const _Rhs& _opd) const 
	{ return exprs::Expr::EwiseBinaryExpr<_Myt, _Rhs, EwiseMinOp>(*this, _opd); }
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator* (const _Rhs& _opd) const 
	{ return exprs::Expr::EwiseBinaryExpr<_Myt, _Rhs, EwiseMulOp>(*this, _opd); }
	template<typename _Rhs> MATRICE_GLOBAL_INL auto operator/ (const _Rhs& _opd) const 
	{ return exprs::Expr::EwiseBinaryExpr<_Myt, _Rhs, EwiseDivOp>(*this, _opd); }
	template<typename _Rhs> MATRICE_GLOBAL_INL auto mul(const _Rhs& rhs) const 
	{ return exprs::Expr::MatBinaryExpr<_Myt, _Rhs, MatMulOp>(*this, rhs); }
	MATRICE_HOST_FINL auto inv() const 
	{ return exprs::Expr::MatUnaryExpr<_Myt, MatInvOp>(*this); }
	MATRICE_HOST_FINL auto inv(const_myt_ref _rhs) 
	{ return exprs::Expr::MatUnaryExpr<_Myt, MatInvOp>(_rhs, *this); }
	MATRICE_HOST_FINL auto transpose() 
	{ return exprs::Expr::MatUnaryExpr<_Myt, MatTrpOp>(*this); }
	MATRICE_GLOBAL_FINL auto normalize(value_t _val = inf) 
	{ return ((std::abs(_val) < eps ? value_t(1) : value_t(1) / (_val == inf ? max() : _val))*(*this)); }
	MATRICE_GLOBAL_INL exprs::Expr::MatBinaryExpr<_Myt, _Myt, MatMulOp> spread();
#pragma endregion

#pragma region <!-- Triggers for Suspended Expression -->
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const exprs::Expr::EwiseBinaryExpr<_Args...>& _Ex){ return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const exprs::Expr::MatBinaryExpr<_Args...>& _Ex) { return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
	template<typename... _Args> MATRICE_GLOBAL_INL auto& operator= (const exprs::Expr::MatUnaryExpr<_Args...>& _Ex) { return (*static_cast<_Derived*>(&_Ex.assign(*this))); }
#pragma endregion

	///<brief> in-time matrix arithmetic </brief>
	MATRICE_GLOBAL_FINL value_t max() const { return (*std::max_element(m_data, m_data + size())); }
	MATRICE_GLOBAL_FINL value_t min() const { return (*std::min_element(m_data, m_data + size())); }
	MATRICE_GLOBAL_FINL value_t sum() const { return (reduce<value_t>(begin(), end())); }
	MATRICE_GLOBAL_FINL value_t det() const { return (det_impl(*static_cast<const _Derived*>(this))); }
	MATRICE_GLOBAL_FINL value_t trace() const { return (reduce(begin(), end(), cols() + 1)); }
	template<class _Rhs> MATRICE_GLOBAL_FINL value_t dot(const _Rhs& _rhs) const { return (this->operator*(_rhs)).sum(); }
	MATRICE_GLOBAL_FINL value_t norm_2() const { auto ans = dot(*this); return (ans > eps ? std::sqrt(ans) : inf); }
	template<typename _Rhs, typename _Ret = _Rhs> MATRICE_HOST_INL _Ret& solve(_Rhs& b) { typename Solver_<value_t>::Linear<_M, _N, AUTO> solver(*static_cast<_Derived*>(this)); return solver.solve(b); }
	MATRICE_HOST_INL auto solve() { typename Solver_<value_t>::Linear<_M, _N, SVD> solver; return solver(*static_cast<_Derived*>(this)); }

	///<brief> properties </brief>
	__declspec(property(get = _Prop_format_getter)) format_t format;
	constexpr MATRICE_HOST_FINL format_t _Prop_format_getter() const { return m_format; }
	__declspec(property(get = _Prop_empty_getter)) bool empty;
	constexpr MATRICE_HOST_FINL bool _Prop_empty_getter() const { return (size() == 0); }
	constexpr MATRICE_GLOBAL_FINL type_t& storage_type() { return m_type; }
	constexpr MATRICE_GLOBAL_FINL type_t storage_type() const { return m_type; }

protected:
	using base_t::m_rows;
	using base_t::m_cols;
	using base_t::m_data;
	format_t m_format = RowMajor;
	type_t m_type = type_t::GDs;
public:
	Storage m_storage;
};
MATRICE_NAMESPACE_END_TYPES