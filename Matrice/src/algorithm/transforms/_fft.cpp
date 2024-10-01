/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library for
Geometric Optical Sensing and Visual Intelligence.
Copyright(C) 2018-2024, Zhilong(Dgelom) Su, all rights reserved.

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
#include "algs/transform/_fft.hpp"
#if MATRICE_MATH_KERNEL == MATRICE_USE_MKL
#include <fftw3_mkl.h>
#endif

MATRICE_ALG_BEGIN()
_DETAIL_BEGIN

template<typename _Ty>
_Fft_descriptor<_Ty>::status_type 
_Fft_descriptor<_Ty>::_Forward() noexcept
{
	return status_type();
}

template<typename _Ty>
_Fft_descriptor<_Ty>::status_type
_Fft_descriptor<_Ty>::_Backward() noexcept
{
	return status_type();
}

template<> _Fft_descriptor<float>::status_type _Fft_descriptor<float>::_Forward() noexcept;
template<> _Fft_descriptor<float>::status_type _Fft_descriptor<float>::_Backward() noexcept;
template<> _Fft_descriptor<double>::status_type _Fft_descriptor<double>::_Forward() noexcept;
template<> _Fft_descriptor<double>::status_type _Fft_descriptor<double>::_Backward() noexcept;

_DETAIL_END
MATRICE_ALG_END()