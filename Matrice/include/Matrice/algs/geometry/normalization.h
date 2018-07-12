#pragma once

#include "../../core/matrix.h"
#include "../../arch/ixpacket.h"

MATRICE_NAMESPACE_BEGIN_

template<typename _Data_type, typename value_t = typename _Data_type::value_t> 
class normalization MATRICE_NONHERITABLE
{
	using data_t = _Data_type;
	using packet_t = simd::Packet_<value_t, 4>;
public:
	MATRICE_HOST_INL normalization(const data_t& _Data) : m_data(_Data) {}

	MATRICE_HOST_INL auto operator() () const {
		if (_My_future.valid()) _My_future.get();

		auto _Sx = std::sqrt(2/m_params(1)), _Sy = std::sqrt(2/m_params(3));
		return types::Matrix_<value_t, 3, 3>{_Sx, 0, -_Sx * m_params(0),
			0., _Sy, -_Sy * m_params(2), 0., 0., 1.};
	}
	MATRICE_HOST_INL auto operator() (data_t& _Data) const {
		if (_My_future.valid()) _My_future.get();

		auto _Mx = m_params(0), _My = m_params(2);
		auto _Sx = std::sqrt(2/m_params(1)), _Sy = std::sqrt(2/m_params(3));

		const auto _Size = m_data.cols();
		std::transform(m_data[0], m_data[0] + _Size, _Data[0], [&](auto _X)->value_t { return (_Sx*(_X - _Mx)); });
		std::transform(m_data[1], m_data[1] + _Size, _Data[1], [&](auto _Y)->value_t { return (_Sy*(_Y - _My)); });
	}
private:
	const data_t& m_data;
	mutable types::Matrix_<value_t, 2, 2> m_params;
	mutable std::future<void> _My_future = std::async(std::launch::async, [&] {
		const auto _N = m_data.cols();
		m_params[0][0] = reduce<value_t>(m_data[0], m_data[0] + _N) / _N;
		m_params[0][1] = reduce<std::minus>(m_data[0], m_data[0] + _N, m_params(0), [&](auto _Val)->value_t {return (_Val*_Val); }) / _N;
	
		m_params[1][0] = reduce<value_t>(m_data[1], m_data[1] + _N) / _N;
		m_params[1][1] = reduce<std::minus>(m_data[1], m_data[1] + _N, m_params(2), [&](auto _Val)->value_t {return (_Val*_Val); }) / _N;
	});
};

_MATRICE_NAMESPACE_END
