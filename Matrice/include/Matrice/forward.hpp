/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for
more detail.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

namespace dgelom {
namespace types {
	template<typename _Ty, int _M, int _N> class Matrix_;
	template<typename _Ty, int _Dim> class Vec_;
	template<typename _Ty> class Vec2_;
}

template<typename _Ty, int _RowsAtCompileTime, int _ColsAtCompileTime=_RowsAtCompileTime> 
using Matrix_ = types::Matrix_<_Ty, _RowsAtCompileTime, _ColsAtCompileTime>;

template<typename _Ty, int _Ndims>
using Vec_ = types::Vec_<_Ty, _Ndims>;

template<typename _Ty>
using Vec2_ = types::Vec2_<_Ty>;

template<typename _Ty, int _Size>
using Array_ = types::Matrix_<_Ty, _Size, (_Size > 0 ? 1 : _Size)> ;
}