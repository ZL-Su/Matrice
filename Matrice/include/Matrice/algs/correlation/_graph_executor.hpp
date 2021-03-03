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
	enum computing_type
	{
		sequence = 0,
		parallel = 1,
		instant  = 3,
	};

	enum device_type
	{
		cpu = 0,
		gpu = 1
	};

	struct options_type
	{
		// Specify capacity to prefetch node resources. 
		size_t preload_capacity = 1;

		// Specify computing paradigm.
		size_t exec_mode = computing_type::sequence | device_type::cpu;

		// Specify to update reference layer.
		bool adjacent = false;
	};

	class layer_type
	{
		struct status
		{
			bool processed = false;
			bool active = true;
		};

		// FIELD, data

		// FIELD, operator

		// FIELD, status
		status _Mystatus;

		// FIELD, parallel computing index
		size_t _Myidx = 0;
	public:
		/**
		 * \brief METHOD, detach this layer from the computational graph.
		 */
		void detach() noexcept {
			_Mystatus.active = false;
		}

		/**
		 * \brief METHOD, join this layer into the computational graph.
		 */
		void join() noexcept {
			_Mystatus.active = true;
		}

		/**
		 * \brief METHOD, check if this layer is detached or not.
		 */
		bool is_active() const noexcept {
			return _Mystatus.active;
		}

		/**
		 * \brief METHOD, check if this layer is processed or not.
		 */
		bool is_processed() const noexcept {
			return _Mystatus.processed;
		}
		
	};

	/**
	 * \brief METHOD, add a computing node.
	 */
	_Myt& add_layer(const layer_type& node) {
		_Mylayers.push_back(node);
		return (*this);
	}

private:
	std::vector<layer_type> _Mylayers;
	
};

_DETAIL_END
MATRICE_ALG_END(corr)