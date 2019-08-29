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

template<typename _Mty>
MATRICE_GLOBAL_FINL decltype(auto) ref(_Mty& cont, bool t=false)noexcept {
	return Ref<_Mty>(cont, t);
}
template<typename _Mty>
MATRICE_GLOBAL_FINL decltype(auto) ref_t(_Mty& cont) noexcept {
	return ref(cont, true);
}
DGE_MATRICE_END