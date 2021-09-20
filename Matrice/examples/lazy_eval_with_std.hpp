/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
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
**********************************************************************/
#pragma once
#include <ranges>
#include <core/matrix.h>

DGE_MATRICE_BEGIN
namespace example {
void _Lazy_eval_with_std_ranges() {
	dgelom::Matrix<int> _Vec = { 1, 2, 3, 4, 5, 6, 7, 8 };
	auto&& _Exp = _Vec
		| std::ranges::views::filter([](const auto& x) {return x % 2 == 0; }) // 2, 4, 6, 8
		| std::ranges::views::transform([](auto& x) {return x += 2; }) // 4, 6, 8, 10
		| std::ranges::views::take(3); // 4, 6, 8
	decltype(_Vec) _Res(3, 1);
	for (auto _It = _Exp.begin(); _It != _Exp.end(); ++_It) {
		const auto _Idx = _It - _Exp.begin();
		_Res(_Idx) = *_It;
	}
}
}
DGE_MATRICE_END