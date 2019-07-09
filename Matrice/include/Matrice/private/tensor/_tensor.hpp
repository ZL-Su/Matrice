/*********************************************************************
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
***********************************************************************/
#pragma once
#include "../_matrix_base.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Ty>
struct index<_Ty, tag::_Tensor_tag> {
	static_assert(is_integral_v<_Ty>, "_Ty must be an integer type");

	MATRICE_GLOBAL_INL index(_Ty d, _Ty h, _Ty w) noexcept {
		data[0] = d, data[1] = h, data[2] = w;
	}
	MATRICE_GLOBAL_INL _Ty value() noexcept {
		return 0;
	}

	std::array<_Ty, 3> data;
};
/**
 *\brief CLASS TEMPLATE dgelom::tensor prototype
 *\param <_Ty> data type
 *\param <_Depth> tensor depth which grows vertically
 */
template<typename _Ty, size_t _Depth>
class _Tensor 
	: public types::Base_<_Tensor<_Ty, _Depth>, tensor_traits<_Tensor<_Ty, _Depth>>>
{
	using _Myt = _Tensor;
	using _Mytraits = tensor_traits<_Myt>;
	using _Mybase = types::Base_<_Myt, _Mytraits>;
public:
	enum { Size = 0, CompileTimeRows = 0, CompileTimeCols = 0 };
	using typename _Mybase::value_type;
	//using _Mybase::Base_;
	using _Mybase::operator=;
	using _Mybase::operator();

	/**
	 *\brief default constructor
	 */
	_Tensor() noexcept 
		: _Mybase() {
	}
	/**
	 *\brief constructor
	 *\param [h, w] rows and cols of each tensor cell
	 */
	_Tensor(size_t h, size_t w) noexcept
		:_Mybase(_Depth*h, w) {
		_Mybase::_Myshape = { 1, _Depth, h, w };
		_Mybase::_Flush_view_buf();
	}
	/**
	 *\brief constructor with a value initialization
	 *\param [h, w] rows and cols of each tensor cell
	 *\param [_Val] initial value 
	 */
	_Tensor(size_t h, size_t w, value_type _Val) noexcept
		:_Tensor(h, w) {
		_Mybase::operator=(_Val);
	}
	/**
	 *\brief copy constructor
	 *\param [oth] an other tensor
	 */
	_Tensor(const _Myt& oth) noexcept
		:_Mybase(oth),
		m_height(oth.m_height),
		m_width(oth.m_width) {
	}
	/**
	 *\brief move constructor
	 *\param [oth] an other tensor
	 */
	_Tensor(_Myt&& oth) noexcept
		:_Mybase(move(oth)),
		m_height(oth.m_height),
		m_width(oth.m_width) {
	}
	/**
	 *\brief template constructor
	 *\param [args...] argument(s) with any supported type(s) 
	 */
	template<typename... _Args>
	_Tensor(_Args&&... args) noexcept
		:_Mybase(forward<_Args>(args)...) {
		m_width = this->cols();
		m_height = this->rows() / m_depth;
	}

	MATRICE_HOST_INL const value_type& operator()(size_t d, size_t r, size_t c) const noexcept {
#if defined _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < _Depth, "depth index over range.");
#endif
		const auto inner_size = m_width * m_height;
		return (_Mybase::operator[](d*inner_size+r*m_width)[c]);
	}

	MATRICE_HOST_INL value_type& operator()(size_t d, size_t r, size_t c) noexcept {
#if defined _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < _Depth, "depth index over range.");
#endif
		const auto inner_size = m_width * m_height;
		return (_Mybase::operator[](d*inner_size + r * m_width)[c]);
	}

	MATRICE_HOST_INL decltype(auto) slice(size_t d) const noexcept {
#if defined _DEBUG || MATRICE_DEBUG
		DGELOM_CHECK(d < _Depth, "depth index over range.");
#endif
		const auto w = _Mybase::_Myshape.get(3);
		const auto h = _Mybase::_Myshape.get(2);
		const auto r0 = d*h, r1 = (d+1) * h;
		const auto c0 = w, c1 = w;
		return _Mybase::block(c0, c1, r0, r1);
	}

	MATRICE_HOST_INL void __create_impl(size_t h, size_t w) {
		m_width = w, m_height = h;
		_Mybase::_Myshape = { _Depth, 1, m_height, m_width };
		_Mybase::m_cols = m_width;
		_Mybase::m_rows = _Depth * m_height;
		_Mybase::allocator().create(_Mybase::m_rows, _Mybase::m_cols);
		_Mybase::m_data = _Mybase::allocator().data();
		_Mybase::_Flush_view_buf();
	}

private:
	size_t m_width, m_height;
	size_t m_depth = _Mytraits::depth;
};

template<typename _Ty, size_t _Depth>
struct tensor_traits<_Tensor<_Ty, _Depth>> {
	using type = _Ty;
	using category = tag::_Tensor_tag;
	static constexpr auto depth = _Depth;
	static constexpr auto _M = 0, _N = 0;
	static constexpr bool Is_base = std::false_type::value;
};

_DETAIL_END
DGE_MATRICE_END