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

// Forward declarations
template<typename _Ty, class _Ip> class _Image;

/// <summary>
/// \brief CLASS, graph based correlation pipeline. 
/// </summary>
template<typename _Ty>
class _Graph_executor {
	using _Myt = _Graph_executor;
public:
	using value_type = _Ty;

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

	/// <summary>
	/// \brief Set the options to control the graph's computation. 
	/// </summary>
	struct options_type
	{
		// Specify window size of every node can view.
		size_t node_field = 15;

		// Specify capacity to prefetch node resources. 
		size_t preload_capacity = 1;

		// Specify computing paradigm.
		size_t exec_mode = computing_type::sequence | device_type::cpu;

		// Specify to update reference layer.
		bool adjacent = false;
	};

	/// <summary>
	/// \brief CLASS layer_type wraps an image to be processed.
	/// Each layer contains a set of nodes, each of them 
	/// </summary>
	class layer_type
	{
		struct status
		{
			bool processed = false;
			bool active = true;
		};

		// FIELD, smooth image data
		_Image<value_type, void> _Mydata;

		// FIELD, operator to 

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
	 * \brief METHOD, add a computing layer.
	 */
	_Myt& add_layer(const layer_type& node) {
		_Mylayers.push_back(node);
		return (*this);
	}

	/**
	 * \brief METHOD, insert a computing layer.
	 */
	_Myt& insert_layer(size_t pos, layer_type&& node) {
		decltype(auto) _Where = _Mylayers.cbegin() + pos;
		_Mylayers.insert(_Where, node);
		return (*this);
	}

private:
	std::vector<layer_type> _Mylayers;
	
};

_DETAIL_END
MATRICE_ALG_END(corr)