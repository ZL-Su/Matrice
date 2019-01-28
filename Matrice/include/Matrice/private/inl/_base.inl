/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../_matrix_base.hpp"
#include "../_range.h"
#include "../../thread/_thread.h"

DGE_MATRICE_BEGIN
_TYPES_BEGIN
template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs>
auto& Base_<_Derived, _Traits, _Type>::inplace_sub(const _Rhs& _Right) {
	using packet_t = simd::Packet_<value_type>;
	constexpr auto step = packet_t::size;
	const auto plen = simd::vsize<step>(this->size());

	if constexpr (is_scalar_v<_Rhs>) {
		packet_t pb(_Right);
		for (const auto i : range(0, plen, step)) {
			(packet_t(m_data + i) - pb).unpack(m_data + i);
		}
		for (auto i = plen; i < this->size(); ++i) {
			m_data[i] -= _Right;
		}
	}
	else if constexpr (is_matrix_v<_Rhs>) {
		DGELOM_CHECK(this->dims() != _Right.dims(), "Inconsistent shapes.");
		const auto rdata = _Right.data();
		for (const auto i : range(0, plen, step)) {
			(packet_t(m_data+i)-packet_t(rdata+i)).unpack(m_data+i);
		}
		for (auto i = plen; i < this->size(); ++i) {
			m_data[i] -= rdata[i];
		}
	}
	else {
		parallel_nd([&](const auto& i) { m_data[i] -= _Right(i); });
	}

	return (*this);
}
_TYPES_END
DGE_MATRICE_END