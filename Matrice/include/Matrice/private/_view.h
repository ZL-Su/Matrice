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
#include <valarray>
#include "../util/utils.h"

DGE_MATRICE_BEGIN
_TYPES_BEGIN
#define _VIEW_EWISE_COPY_N(_LEFT, _N)\
DGELOM_CHECK(_N <= size(), "matrix size over range of current view."); \
size_t _Size = (_N==0) ? size() : _N; \
for(size_t _Idx = 0; _Idx < _Size; ++_Idx) { \
	_LEFT(_Idx) = this->operator()(_Idx); \
}
template<typename _Ty, int _Rows, int _cols> class Matrix_;
/**********************************************************************
						       Matrix view base class 
	    Copyright (c) : Zhilong (Dgelom) Su, since 25/Jul/2018
 **********************************************************************/
template<typename _Ty, typename _Derived> class _View_base 
{
#define _VIEW_EWISE_OP(_OPERATION) \
for (difference_type i = 0; i < size(); ++i) {\
	static_cast<_Derived*>(this)->operator()(i) = _OPERATION;\
} return (*this)
#define _MATRICE_DEFVIEW_ARITHOP(OP, NAME) \
template<typename _Rhs> MATRICE_GLOBAL_FINL auto operator##OP(const _Rhs& _Right) { \
	return Expr::EwiseBinaryExpr<_Derived, _Rhs, _Exp_op::_Ewise_##NAME<value_t>>(*static_cast<_Derived*>(this), _Right); \
} \
template<typename _Lhs> friend MATRICE_GLOBAL_FINL auto operator##OP(const _Lhs& _Left, const _Derived& _Right) { \
	return Expr::EwiseBinaryExpr<_Lhs, _Derived, _Exp_op::_Ewise_##NAME<value_t>>(_Left, _Right); \
}
	using _Myt = _View_base;
public:
	using value_t = _Ty;
	using value_type = value_t;
	using difference_type = std::ptrdiff_t;
	using pointer = std::add_pointer_t<value_type>;
	using reference = std::add_lvalue_reference_t<value_type>;
	enum { CompileTimeRows = 0, CompileTimeCols = 0 };
	struct range_type
	{
		// _Rang = {from_x, from_y, end_x, end_y} : [from_x, end_x), [from_y, end_y)
		template<typename _Idx, MATRICE_ENABLE_IF(is_integral_v<_Idx>)>
		MATRICE_GLOBAL_FINL range_type(const initlist<_Idx> _Rang)
			:_My_from_x(*_Rang.begin()), _My_from_y(*(_Rang.begin()+1)),
			 _My_end_x(*(_Rang.begin()+2)), _My_end_y(*(_Rang.begin()+3)) {}

		template<typename _Idx, MATRICE_ENABLE_IF(is_integral_v<_Idx>)>
		MATRICE_GLOBAL_FINL range_type& operator= (const initlist<_Idx> _Rang) {
			_My_from_x = *_Rang.begin(), _My_from_y = *(_Rang.begin() + 1);
			_My_end_x = *(_Rang.begin() + 2), _My_end_y = *(_Rang.begin() + 3);
		}

		MATRICE_GLOBAL_FINL auto& begin_x() { return _My_from_x; }
		MATRICE_GLOBAL_FINL auto& begin_y() { return _My_from_y; }
		MATRICE_GLOBAL_FINL auto& end_x() { return _My_end_x; }
		MATRICE_GLOBAL_FINL auto& end_y() { return _My_end_y; }
		MATRICE_GLOBAL_FINL const auto& begin_x() const { return _My_from_x; }
		MATRICE_GLOBAL_FINL const auto& begin_y() const { return _My_from_y; }
		MATRICE_GLOBAL_FINL const auto& end_x() const { return _My_end_x; }
		MATRICE_GLOBAL_FINL const auto& end_y() const { return _My_end_y; }
		MATRICE_GLOBAL_FINL size_t size() const { 
			return(_My_end_x - _My_from_x)*(_My_end_y - _My_from_y); 
		}
		
		difference_type _My_from_x, _My_from_y;
		difference_type _My_end_x, _My_end_y;
	};

	MATRICE_GLOBAL_INL _View_base(pointer _Ptr, size_t _Size, size_t _Stride, size_t _Offset)
		:_My_data(_Ptr), _My_size(_Size), _My_stride(_Stride), _My_offset(_Offset) {}

	MATRICE_GLOBAL_FINL reference operator[] (size_t i) {
		return _My_data[i*_My_stride];
	}
	MATRICE_GLOBAL_FINL const reference operator[] (size_t i) const {
		return _My_data[i*_My_stride];
	}
	MATRICE_GLOBAL_FINL reference operator() (size_t i) {
		return _My_data[i*_My_stride];
	}
	MATRICE_GLOBAL_FINL const reference operator() (size_t i) const {
		return _My_data[i*_My_stride];
	}

	MATRICE_GLOBAL_FINL size_t size() const { return (static_cast<const _Derived*>(this)->size()); }
	MATRICE_GLOBAL_FINL size_t rows() const { return (static_cast<const _Derived*>(this)->rows()); }
	MATRICE_GLOBAL_FINL size_t cols() const { return (static_cast<const _Derived*>(this)->cols()); }
	MATRICE_GLOBAL_FINL constexpr auto shape() const { 
		return basic_shape_t({ 1, 1, rows(), cols() });
	}
	MATRICE_GLOBAL_FINL constexpr auto dims() const { return basic_shape_t({ 1,1,rows(), cols() }); }
	MATRICE_GLOBAL_FINL void create(size_t, size_t) {}
	MATRICE_GLOBAL_FINL value_t sum() const {
		value_t _Ret = 0;
#pragma omp parallel if (size() > 100)
		{
#pragma omp for reduction (+ : _Ret)
			for (ptrdiff_t _Idx = 0; _Idx < size(); ++_Idx) {
				_Ret += static_cast<const _Derived*>(this)->operator()(_Idx);
			}
		}
		return (_Ret);
	}

	/**
	 *\brief Copy a scalar value to memory that the view maps to
	 *\param [_Val] an input scalar
	 */
	template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)>
	MATRICE_GLOBAL_INL auto& operator= (const _Ty _Val) {
		_VIEW_EWISE_OP(static_cast<value_type>(_Val));
	}
	/**
	 *\brief Copy from another view
	 *\param [_Oth] another view input
	 */
	MATRICE_GLOBAL_FINL auto& operator= (const _Myt& _Oth) {
		_VIEW_EWISE_OP(_Oth(i));
	}
	/**
	 *\brief Fill view memory from a given pointer
	 *\param [_Data] an input pointer, the size of the pointer pointed memory should not less than the size of the view 
	 */
	MATRICE_GLOBAL_FINL auto& operator= (const pointer _Data) {
		_VIEW_EWISE_OP(_Data[i]);
	}
	/**
	 *\brief Fill view memory from a initializer_list
	 *\param [_L] the size of the list should not less than the size of the view
	 */
	MATRICE_GLOBAL_FINL auto& operator= (const initlist<value_type> _L) {
		_VIEW_EWISE_OP(*(_L.begin() + i));
	}
	/**
	 *\brief Fill view memory from a customer class type
	 *\param [_M] _Mty should have element accessor ::operator(i)
	 */
	template<typename _Mty, MATRICE_ENABLE_IF(is_class_v<_Mty>)>
	MATRICE_GLOBAL_INL auto& operator= (const _Mty& _M) {
		_VIEW_EWISE_OP(_M(i));
	}
	/**
	 *\brief Evaluation from an expression
	 *\param [_Ex] input expression
	 */
	template<typename _Arg> 
	MATRICE_GLOBAL_INL auto& operator= (const Expr::Base_<_Arg>& _Ex) { 
		return (_Ex.assign(*static_cast<_Derived*>(this))); 
	}
	/**
	 *\brief Evaluation this view to an input container.
	 *\param [_Rt] given templated type with a linear accessor: operator(i).
	 *\note The first min(_Rt.size(), size()) elements will be evaluated if these two sizes are not identity.
	 */
	template<typename _Arg>
	MATRICE_GLOBAL_INL void eval_to(_Arg& _Rt) const {
#if _DEBUG
		DGELOM_CHECK(_Rt.size() != 0, "Input _Rt is empty.");
#endif
		const auto _Size = min(_Rt.size(), size());
		const auto _Deri = static_cast<const _Derived*>(this);
		for (difference_type _Idx = 0; _Idx < _Size; ++_Idx) {
			_Rt(_Idx) = _Deri->operator()(_Idx);
		}
	}

	_MATRICE_DEFVIEW_ARITHOP(+, add)
	_MATRICE_DEFVIEW_ARITHOP(-, sub)
	_MATRICE_DEFVIEW_ARITHOP(*, mul)
	_MATRICE_DEFVIEW_ARITHOP(/, div)

protected:
	pointer _My_data;
	size_t  _My_size;
	size_t  _My_stride;
	size_t  _My_offset;

#undef _VIEW_EWISE_OP
};

/**********************************************************************
						      Row view for Matrix 
	    Copyright (c) : Zhilong (Dgelom) Su, since 25/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_arithmetic_v<_Ty>)>
class _Matrix_rview MATRICE_NONHERITABLE : public _View_base<_Ty, _Matrix_rview<_Ty>>
{
	using _Base = _View_base<_Ty, _Matrix_rview<_Ty>>;
	using _Base::_My_data;
	using _Base::_My_size;
	using _Base::_My_stride;
	using _Base::_My_offset;
public:
	using typename _Base::pointer;
	using typename _Base::reference;
	using typename _Base::value_t;
	using _Base::operator+;
	using _Base::operator-;
	using _Base::operator*;
	using _Base::operator/; 
	using _Base::operator=;

	MATRICE_GLOBAL_FINL _Matrix_rview(pointer _Ptr, size_t _Size, size_t _Stride = 1, size_t _Offset = 1)
		:_Base(_Ptr, _Size, _Stride, _Offset) {}

	MATRICE_GLOBAL_FINL constexpr auto rows() const { return 1; }
	MATRICE_GLOBAL_FINL constexpr auto cols() const { return _My_size; }
	MATRICE_GLOBAL_FINL constexpr auto size() const { return _My_size; }

	/**
	 *\Retrieve the first _N element into a true static row-matrix.
	 */
	template<size_t _N=0> 
	MATRICE_GLOBAL_FINL auto eval() const {
		Matrix_<value_t, min(_N, 1), _N> _Ret(1, _My_size);
		_VIEW_EWISE_COPY_N(_Ret, _N)
		return forward<decltype(_Ret)>(_Ret);
	}
};

/**********************************************************************
						     Column view for Matrix 
	    Copyright (c) : Zhilong (Dgelom) Su, since 25/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_arithmetic_v<_Ty>)>
class _Matrix_cview MATRICE_NONHERITABLE : public _View_base<_Ty, _Matrix_cview<_Ty>>
{
	using _Base = _View_base<_Ty, _Matrix_cview<_Ty>>;
	using _Base::_My_data;
	using _Base::_My_size;
	using _Base::_My_stride;
	using _Base::_My_offset;
public:
	using typename _Base::pointer;
	using typename _Base::reference;
	using typename _Base::value_t;
	using _Base::operator+;
	using _Base::operator-;
	using _Base::operator*;
	using _Base::operator/;
	using _Base::operator=;

	MATRICE_GLOBAL_FINL _Matrix_cview(pointer _Ptr, size_t _Size, size_t _Stride = 1, size_t _Offset = 1)
		:_Base(_Ptr, _Size, _Stride, _Offset) {}

	MATRICE_GLOBAL_FINL constexpr auto rows() const { return _My_size; }
	MATRICE_GLOBAL_FINL constexpr auto cols() const { return 1; }
	MATRICE_GLOBAL_FINL constexpr auto size() const { return _My_size; }

	/**
	 *\Retrieve the first _N element into a true static column-matrix.
	 */
	template<size_t _N=0> 
	MATRICE_GLOBAL_FINL auto eval() const {
		Matrix_<value_t, _N, min(_N, 1)> _Ret(_My_size, 1);
		_VIEW_EWISE_COPY_N(_Ret, _N)
		return std::forward<decltype(_Ret)>(_Ret);
	}
};

/**********************************************************************
						      Block view for Matrix 
	    Copyright (c) : Zhilong (Dgelom) Su, since 25/Jul/2018
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_arithmetic_v<_Ty>)>
class _Matrix_block MATRICE_NONHERITABLE : public _View_base<_Ty, _Matrix_block<_Ty>>
{
	using _Base = _View_base<_Ty, _Matrix_block<_Ty>>;
	using typename _Base::range_type;
	using typename _Base::difference_type;
	using _Base::_My_data;   //begin of this block data
	using _Base::_My_size;   //cols of this block
	using _Base::_My_stride; //cols of source matrix
	using _Base::_My_offset; //offset relative to original matrix data
public:
	using typename _Base::pointer;
	using typename _Base::reference;
	using typename _Base::value_t;
	using _Base::operator+;
	using _Base::operator-;
	using _Base::operator*;
	using _Base::operator/;
	using _Base::operator=;

	MATRICE_GLOBAL_FINL _Matrix_block(pointer _Ptr, size_t _Cols, const range_type _Range)
		: _Base(_Ptr + _Range.begin_x() + _Range.begin_y()*_Cols, 
			_Range.end_x() - _Range.begin_x(), _Cols, 
			_Range.begin_x() + _Range.begin_y()*_Cols),
			_My_origin(_Ptr), _My_range(_Range) {}

	//i zero-based local row index
	MATRICE_GLOBAL_FINL pointer operator[] (difference_type i) {
		return (_My_data + i * _My_stride);
	}
	MATRICE_GLOBAL_FINL const pointer operator[] (difference_type i) const {
		return (_My_data + i * _My_stride);
	}
	//i zero-based local linear index
	MATRICE_GLOBAL_INL reference operator() (difference_type i) {
		auto _Row = i / _My_size;
		return this->operator[](_Row)[i - _My_size * _Row];
	}
	MATRICE_GLOBAL_INL const reference operator() (difference_type i) const {
		auto _Row = i / _My_size;
		return this->operator[](_Row)[i - _My_size * _Row];
	}

	MATRICE_GLOBAL_FINL auto rows() const { return _My_range.end_y() - _My_range.begin_y(); }
	MATRICE_GLOBAL_FINL auto cols() const { return _My_size; }
	MATRICE_GLOBAL_FINL auto size() const { return rows()*cols(); }

	/**
	 * \Evaluate a view to a true matrix.
	 */
	template<int _M = 0, int _N = _M>
	MATRICE_GLOBAL_FINL auto eval() const {
		Matrix_<value_t, _M, _N> _Ret(rows(), cols());
		_VIEW_EWISE_COPY_N(_Ret, _Ret.size());
		return forward<decltype(_Ret)>(_Ret);
	}
private:
	pointer _My_origin;
	range_type _My_range;
};

/**********************************************************************
								    Tensor CHW view
		 Copyright (c) : Zhilong (Dgelom) Su, since 24/Jan/2019
 **********************************************************************/
template<typename _Ty, MATRICE_ENABLE_IF(is_arithmetic_v<_Ty>)>
class _Chw_view MATRICE_NONHERITABLE : public _View_base<_Ty, _Chw_view<_Ty>>
{
	using _Myt = _Chw_view;
	using _Mybase = _View_base<_Ty, _Chw_view>;
	using _Mybase::_My_data;   //begin of this block data
	using _Mybase::_My_size;   //cols of this block
	using _Mybase::_My_stride; //cols of source matrix
	using _Mybase::_My_offset; //offset relative to original matrix data
public:
	using typename _Mybase::range_type;
	using typename _Mybase::difference_type;
	using typename _Mybase::reference;
	using typename _Mybase::pointer;
	using typename _Mybase::value_t;

	MATRICE_GLOBAL_FINL _Chw_view() noexcept {}

	/**
	 *\brief Get the view for matrix at (_Myidx, _C)
	 *\param [_C] channel index
	 */
	MATRICE_GLOBAL_FINL auto operator()(size_t _C) const {
		range_type _Range(_C*_Mywsize, 0, (_C +1)*_Mywsize, _Myhsize);
		return _Matrix_block<value_t>(_My_data, _My_stride, _Range);
	}

	/**
	 *\brief Evaluate to a tensor
	 */
	MATRICE_GLOBAL_FINL auto eval() const {
		detail::_Tensor<value_t> _Ret{ 1, _Mycsize, _Myhsize, _Mywsize };
		return std::move(_Ret);
	}
private:
	size_t _Myidx, _Mycsize;
	size_t _Myhsize, _Mywsize;
};

#undef _VIEW_EWISE_COPY_N
_TYPES_END
DGE_MATRICE_END
