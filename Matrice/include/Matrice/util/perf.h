/**********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2022, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#pragma once

namespace dgelom {
namespace perf {

template<typename _Ty, int _M, int _N, 
	template<typename _Ty, int _M, int _N> class Matrix_>
auto gflops(const Matrix_<_Ty, _M, _N>& x) noexcept {
	return 2 * x.size() * 1.E-9 * (sizeof(_Ty) / 4.f);
}

} // namespace perf
} // namespace dgelom
