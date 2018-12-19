/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#include <stdexcept>
#include "../../include/Matrice/private/_decl_dev_funcs.h"

#if (defined __enable_cuda__ && !defined __disable_cuda__)
#include <cuda_runtime.h>

#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)

MATRICE_PRIVATE_BEGIN

template<int _Opt> void _Device_sync()
{
	cudaError_t sts;
	switch (_Opt){
	case 0: sts = cudaDeviceSynchronize(); break;
	default: break;
	}

	if(sts != cudaError_t::success) 
		throw std::runtime_error("Fail to device thread synchronization.");
}

template void _Device_sync<0>();

MATRICE_PRIVATE_END
#endif