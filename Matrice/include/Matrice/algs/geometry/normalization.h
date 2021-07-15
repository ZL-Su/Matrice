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

#include "core/matrix.h"
#ifdef MATRICE_SIMD_ARCH
#include "arch/simd.h"
#endif // MATRICE_SIMD_ARCH

DGE_MATRICE_BEGIN

template<typename _Data_type, 
	typename value_t = typename _Data_type::value_t>
class normalization MATRICE_NONHERITABLE
{
	using data_t = _Data_type;
#ifdef MATRICE_SIMD_ARCH
	using packet_t = simd::Packet_<value_t, 4>;
#endif
public:
	MATRICE_HOST_INL normalization(const data_t& _Data) : m_data(_Data) {}

	MATRICE_HOST_INL auto operator() () const {
		if (_My_future.valid()) _My_future.get();

		auto _Sx = sqrt(2/m_params(1)), _Sy = sqrt(2/m_params(3));
		return detail::Matrix_<value_t, 3, 3>{_Sx, 0, -_Sx * m_params(0),
			0., _Sy, -_Sy * m_params(2), 0., 0., 1.};
	}
	MATRICE_HOST_INL auto operator() (data_t& _Data) const {
		if (_My_future.valid()) _My_future.get();

		auto _Mx = m_params(0), _My = m_params(2);
		auto _Sx = sqrt(2/m_params(1)), _Sy = sqrt(2/m_params(3));

		const auto _Size = m_data.cols();
		std::transform(m_data[0], m_data[0] + _Size, _Data[0], [&](auto _X)->value_t { return (_Sx*(_X - _Mx)); });
		std::transform(m_data[1], m_data[1] + _Size, _Data[1], [&](auto _Y)->value_t { return (_Sy*(_Y - _My)); });
	}
private:
	const data_t& m_data;
	mutable detail::Matrix_<value_t, 2, 2> m_params;
	mutable std::future<void> _My_future = std::async(std::launch::async, [&] {
		const auto _N = m_data.cols();
		m_params[0][0] = reduce<value_t>(m_data[0], m_data[0] + _N) / _N;
		m_params[0][1] = reduce<std::minus>(m_data[0], m_data[0] + _N, m_params(0),
			[&](auto _Val)->value_t {return (_Val*_Val); }) / _N;
	
		m_params[1][0] = reduce<value_t>(m_data[1], m_data[1] + _N) / _N;
		m_params[1][1] = reduce<std::minus>(m_data[1], m_data[1] + _N, m_params(2),
			[&](auto _Val)->value_t {return (_Val*_Val); }) / _N;
	});
};

DGE_MATRICE_END
