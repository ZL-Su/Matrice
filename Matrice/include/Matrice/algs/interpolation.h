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

#include "interpolation\_interpolation.h"

MATRICE_NAMESPACE_BEGIN_
enum {
	bilinear = algs::INTERP | algs::BILINEAR,
	bcspline = algs::INTERP | algs::BSPLINE | algs::BICUBIC,
	bqspline = algs::INTERP | algs::BSPLINE | algs::BIQUINTIC,
	bsspline = algs::INTERP | algs::BSPLINE | algs::BISEPTIC,
}; //interpolation algorithms

/*******************************************************************
	              Unified Interpolation Interface
	    Copyright (c) : Zhilong (Dgelom) Su, since 31/Jul/2018
 ******************************************************************/
template<typename _Ty, size_t _Options = bcspline, typename = std::enable_if_t<std::is_scalar_v<_Ty>>>
using interpolation = algs::Interpolation<_Ty, _Options>;
_MATRICE_NAMESPACE_END
