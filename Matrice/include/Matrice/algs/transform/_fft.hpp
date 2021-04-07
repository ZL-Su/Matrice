/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
3D Vision and Photo-Mechanics.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#pragma once

#include "core.hpp"

MATRICE_ALG_BEGIN()
_DETAIL_BEGIN
template<typename _Ty>
class _Fft_descriptor {
	using _Myt = _Fft_descriptor;
public:
	using value_type = _Ty;
	using pointer = add_pointer_t<value_type>;

	struct status_type
	{
		long value;
	};

	struct options_type
	{
		bool inplace = true;
		bool padding = false; // false for no zero padding

		value_type fwdscale = 1;
		value_type bwdscale = 1;

		size_t distance[2]{ 0, 0 };

		size_t threads = 1; // up to 27
	};

	_Fft_descriptor(const pointer _Src, size_t _Size) noexcept 
		:_Mysrc(_Src), _Mysize(_Size) {

	}

	template<typename _Cont>
	_Fft_descriptor(const _Cont& _Src) noexcept
		:_Mysrc(_Src.data()), _Mysize(_Src.size()) {

	}

	MATRICE_HOST_INL decltype(auto) options() const noexcept {
		return (_Myoptions);
	}
	MATRICE_HOST_INL decltype(auto) options() noexcept {
		return (_Myoptions);
	}

private:
	status_type _Forward() noexcept;
	status_type _Backward() noexcept;

protected:
	options_type _Myoptions;

	size_t  _Mysize;
	pointer _Mysrc = nullptr;

	pointer _Myreal = nullptr;
	pointer _Myimag = nullptr;

};
_DETAIL_END

template<typename _Ty>
using fft_t = detail::_Fft_descriptor<_Ty>;

MATRICE_ALG_END()