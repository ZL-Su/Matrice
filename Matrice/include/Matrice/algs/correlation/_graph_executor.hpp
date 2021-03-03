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
***********************************************************************/
#pragma once

#include <vector>
#include "_correlation_traits.h"

MATRICE_ALG_BEGIN(corr)
_DETAIL_BEGIN

/// <summary>
/// \brief CLASS, graph based correlation pipeline. 
/// </summary>
class _Graph_executor {
	using _Myt = _Graph_executor;
public:
	enum class computing_mode
	{
		sequence = 0,
		parallel = 1,
		instant  = 3,
		host     = 10,
		device   = 11
	};

		//

		// Specify computing mode.
		computing_mode exec_mode = computing_mode::sequence;
	};

	struct node_type
	{
		// data

		// operator
	};

private:
	std::vector<node_type> _Mynodes;
	
};

_DETAIL_END
MATRICE_ALG_END(corr)