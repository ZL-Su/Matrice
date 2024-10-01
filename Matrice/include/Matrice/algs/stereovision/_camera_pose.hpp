/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2022, Zhilong(Dgelom) Su, all rights reserved.

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
#pragma once

#include "core/vector.h"
#include "algs/geometry.hpp"

MATRICE_ALG_BEGIN(vision)

/// <summary>
/// \brief CLASS TEMPLATE camera_pose
/// encapsulates the camera rotation and translation relative to a
/// reference coordinate frame.
/// </summary>
/// <typeparam name="_Ty">double or float type</typeparam>
template<typename _Ty = double_t> 
class camera_pose {
public:
	using value_type = _Ty;
	template<uint8_t _Dim>
	using vector = auto_vector_t<value_type, _Dim>;

	explicit camera_pose() noexcept = default;
	camera_pose(initlist<value_type> rt) noexcept {
		_Myr.from(rt.begin());
		if (rt.size() == 6) {
			_Myt.from(rt.begin() + 3);
		}
	}

	/// <summary>
	/// \brief Move the camera to a new position '_Myt + dt'.
	/// </summary>
	/// <param name="dt"> Incremental translation vector to be moved. </param>
	MATRICE_GLOBAL_INL decltype(auto) move(const vector<3>& dt) noexcept {
		_Myt = _Myt + dt;
		return (*this);
	}

	/// <summary>
	/// \brief Get 3x3 rotation matrix.
	/// </summary>
	MATRICE_GLOBAL_INL auto R() const noexcept {
		return rodrigues(_Myr);
	}

	/// <summary>
	/// \brief Get 3x1 translation vector.
	/// </summary>
	MATRICE_GLOBAL_INL decltype(auto)t()const noexcept {
		return (_Myt);
	}

	/// <summary>
	/// \brief Get 3x4 transformation matrix
	/// </summary>
	/// <typeparam name="_Ty"></typeparam>
	MATRICE_GLOBAL_INL auto T() const noexcept {
		return concat(R(), _Myt);
	}

private:
	vector<3> _Myr; //relative rotation vector
	vector<3> _Myt; //relative translation vector
};

MATRICE_ALG_END(vision)