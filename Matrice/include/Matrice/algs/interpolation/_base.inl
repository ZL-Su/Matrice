/*********************************************************************
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
***********************************************************************/
#pragma once

MATRICE_ALGS_BEGIN
template<typename _Ty, typename _Derived> template<typename... _Args>
MATRICE_GLOBAL_FINL InterpBase_<_Ty, _Derived>::InterpBase_(const _Args&... args) {
	this->_Bspline_coeff(args...);
}

template<typename _Ty, typename _Derived> template<typename... _Args>
MATRICE_GLOBAL_FINL void InterpBase_<_Ty, _Derived>::_Bspline_coeff(const _Args&... args) {
	static_cast<derived_t*>(this)->_Bspline_coeff(args...);
}
MATRICE_ALGS_END