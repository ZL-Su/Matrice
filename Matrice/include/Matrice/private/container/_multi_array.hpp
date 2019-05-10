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

template<typename _Ty, size_t _N>
class multi_array {
	using _Myt = multi_array;
	using _Myidx = Matrix_<size_t, _N, 1>;
	using _Mybuf = Matrix_<_Ty, 0, 0>;
public:
	static constexpr auto dimension = _N;
	using value_type = typename _Mybuf::value_type;
	using pointer = typename _Mybuf::pointer;
	using index_type = _Myidx;
	multi_array() {
	}
	multi_array(const _Myidx& _Shape) {
		_Alloc(_Shape);
	}

	MATRICE_HOST_INL _Myt& resize(const _Myidx& _Shape) {
		_Alloc(_Shape);
	}
	MATRICE_HOST_INL size_t size() const {
		return _Data.size();
	}
	pointer data() { return _Data.data(); }
	const pointer data() const { return _Data.data(); }
	value_type& operator()(const _Myidx& _idx) { 
		return _Data(_Index(_idx)); 
	}
	const value_type& operator()(const _Myidx& _idx) const {
		return _Data(_Index(_idx));
	}
	value_type& operator[](size_t _idx) {
		return _Data(_idx);
	}
	const value_type& operator()(size_t _idx) const {
		return _Data(_idx);
	}

private:
	MATRICE_HOST_INL void _Alloc(const _Myidx& _shape) {
		size_t _Size = 1;
		for (int d = dimension - 1; d >= 0; --d) {
			_Stride(d) = _Size;
			_Size *= _shape(d);
		}
		if constexpr (dimension == 2)
			_Data.create(_shape(0), _shape(1), zero<value_type>);
		else {
			_Data.create(_Size, 1, zero<value_type>);
		}
	}
	MATRICE_HOST_INL size_t _Index(const _Myidx& _index) const {
		return ((_Stride*_index).sum());
	}

	Matrix_<int, dimension, 1> _Stride;
	_Mybuf _Data;
};

DGE_MATRICE_END