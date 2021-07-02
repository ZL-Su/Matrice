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
#include "util/_macros.h"
#include "_type_traits.h"
#include "_shape.hpp"

DGE_MATRICE_BEGIN
template<typename _Mty> class Ref {
public:
	using element_type = _Mty;
	Ref() noexcept = delete;
	MATRICE_GLOBAL_FINL Ref(element_type& inst, bool t=false)noexcept
		:_Myinst(inst), _Myreqt(t) {
	}

	MATRICE_GLOBAL_FINL operator element_type&() noexcept {
		return _Myinst;
	}
	MATRICE_GLOBAL_FINL operator bool() noexcept {
		return _Myreqt;
	}
	MATRICE_GLOBAL_FINL decltype(auto) data() noexcept {
		return _Myinst.data();
	}
	MATRICE_GLOBAL_FINL element_type& get() const noexcept {
		return _Myinst;
	}
	MATRICE_GLOBAL_FINL element_type eval() noexcept {
		if (_Myreqt) {
			return _Myinst.t().eval();
		}
		else {
			return _Myinst;
		}
	}
private:
	element_type& _Myinst;
	bool _Myreqt;
};
template<class... T>
struct is_ref<Ref<T...>> : std::true_type {};

/// <summary>
/// \brief HELPER CLASS, wrap data pointer and shape of a matrix or tensor.
/// </summary>
/// <typeparam name="_Ty">Any scalar type</typeparam>
template<typename _Ty, MATRICE_ENABLE_IF(is_scalar_v<_Ty>)> 
class Ref_ {
public:
	using value_type = _Ty;
	using pointer = value_type*;
	using reference = value_type&;

	template<class _Mty, 
		MATRICE_ENABLE_IF(is_matrix_v<remove_all_t<_Mty>>)>
	explicit Ref_(_Mty& _M) noexcept
		: _Mydata(_M.data()), _Myshape(_M.shape()) {
	}

	MATRICE_GLOBAL_INL decltype(auto)data() const noexcept {
		return _Mydata;
	}
	MATRICE_GLOBAL_INL decltype(auto)data() noexcept {
		return _Mydata;
	}
	MATRICE_GLOBAL_INL decltype(auto)shape() const noexcept {
		return (_Myshape);
	}

private:
	pointer _Mydata;
	const shape_t<3>& _Myshape;
};

template<typename _Mty>
MATRICE_GLOBAL_FINL decltype(auto) ref(_Mty& cont, bool t=false)noexcept {
	return Ref<_Mty>(cont, t);
}
template<typename _Mty>
MATRICE_GLOBAL_FINL decltype(auto) ref_t(_Mty& cont) noexcept {
	return ref(cont, true);
}
DGE_MATRICE_END