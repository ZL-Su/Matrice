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

#include "../../../core"

DGE_MATRICE_BEGIN

/**
 *\brief description
 *\param <_N> dimensions of multi_array object.
 */
template<typename _Ty, size_t _N> class multi_array { 
	static_assert(is_scalar_v<_Ty>, "_Ty must be a scalar type.");

	using _Myt = multi_array;
	using _Mybuf = Matrix_<_Ty, 0, 0>;
public:
	static constexpr auto dimension = _N;
	using pointer = typename _Mybuf::pointer;
	using value_type = typename _Mybuf::value_type;
	using index_type = array_n<size_t, dimension>;
	using reference = std::add_lvalue_reference_t<value_type>;

	multi_array() noexcept {
	}
	multi_array(const index_type& _Shape) noexcept {
		_Alloc(_Shape);
	}
	multi_array(const std::array<size_t, dimension>& _Shape) noexcept {
		_Alloc((typename index_type::pointer)_Shape.data());
	}

	MATRICE_HOST_INL _Myt& resize(const index_type& _Shape) {
		_Alloc(_Shape);
		return (*this);
	}
	MATRICE_HOST_INL _Myt& resize(const std::array<size_t, dimension>& _Shape) {
		_Alloc((typename index_type::pointer)_Shape.data());
		return (*this);
	}
	MATRICE_HOST_INL size_t size() const noexcept {
		return _Data.size();
	}
	MATRICE_HOST_INL pointer data() noexcept { return _Data.data(); }
	MATRICE_HOST_INL const pointer data() const noexcept { return _Data.data(); }
	MATRICE_HOST_INL reference operator()(const index_type& _idx) {
		return _Data.data()[_Index(_idx)];
	}
	MATRICE_HOST_INL const reference operator()(const index_type& _idx) const {
		return _Data.data()[_Index(_idx)];
	}
	MATRICE_HOST_INL reference operator[](size_t _idx) noexcept {
		return _Data.data()[_idx];
	}
	MATRICE_HOST_INL const reference operator()(size_t _idx) const noexcept {
		return _Data.data()[_idx];
	}

	MATRICE_HOST_INL typename _Mybuf::iterator begin() noexcept {
		return _Data.begin(); 
	}
	MATRICE_HOST_INL typename _Mybuf::const_iterator begin() const noexcept {
		return _Data.begin(); 
	}
	MATRICE_HOST_INL typename _Mybuf::iterator end() noexcept {
		return _Data.end(); 
	}
	MATRICE_HOST_INL typename _Mybuf::const_iterator end() const noexcept { 
		return _Data.end(); 
	}

	/**
	 *\brief retrieve data buffer by lvalue reference.
	 */
	MATRICE_HOST_INL _Mybuf& all() & {
		return _Data; 
	}
	/**
	 *\brief retrieve data buffer by rvalue reference.
	 */
	MATRICE_HOST_INL _Mybuf all() && { 
		return move(_Data); 
	}

private:
	MATRICE_HOST_INL void _Alloc(const index_type& _shape) {
		size_t _Size = 1;
		for (int d = dimension - 1; d >= 0; --d) {
			_Stride(d) = _Size;
			_Size *= _shape(d);
		}
		if constexpr (dimension == 2)
			_Data.create(_shape(0), _shape(1), zero<value_type>);
		else
			_Data.create(_Size, zero<value_type>);
	}
	MATRICE_HOST_INL size_t _Index(const index_type& _index)const noexcept {
		return ((_Stride*_index).sum());
	}

private:
	array_n<int, dimension> _Stride;
	_Mybuf _Data;
};

DGE_MATRICE_END