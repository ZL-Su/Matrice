/***************************************************************************
This file is part of Matrice, an effcient and elegant C++ library for SC.
      Copyright(C) 2018, Zhilong (Dgelom) Su (su-zl@seu.edu.cn), 
		                   all rights reserved.

This program is free software : you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/
#pragma once
#include <type_traits>
namespace dgelom {

template<typename _Pixty, typename = typename std::enable_if<std::is_arithmetic<_Pixty>::value>::type>
void _Grayscale_stretch(_Pixty* Img, std::size_t rows, std::size_t cols)
{
	typedef _Pixty                         pixel_t;
	typedef pixel_t*                      iterator;

	const std::size_t N = rows * cols;
	iterator _Begin = Img, _End = Img + N;

	pixel_t _Max = pixel_t(0), _Min = pixel_t(255);
	for (; _Begin != _End; ++_Begin)
	{
		_Max = _Max > *_Begin ? _Max : *_Begin;
		_Min = _Min < *_Begin ? _Min : *_Begin;
	}

	float_t _Scal = 255.f / float_t(_Max - _Min);
	for (_Begin -= N; _Begin != _End; ++_Begin)
		*_Begin = pixel_t(_Scal * (*_Begin - _Min));
}

}
