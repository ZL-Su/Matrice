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
template<typename _Rhs> inline
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
		DGELOM_CHECK(this->dims()==_Right.dims(), "Inconsistent shapes.");
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

template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs, typename> inline
auto& Base_<_Derived, _Traits, _Type>::mul_(const _Rhs& _Right) {
	using packet_t = simd::Packet_<value_type>;
	constexpr auto step = packet_t::size;
	if constexpr (_Rhs::_ctrs_ < step) {
		packet_t pb(_Right.data());
		for (auto r = 0; r < this->rows(); ++r) {
			_Right.data()[r] = simd::reduce(packet_t((*this)[r])*pb);
		}
	}
	else {
		DGELOM_ERROR("Undefined operation.");
	}

	return (_Right);
}

template<typename _Derived, typename _Traits, typename _Type>
template<ttag _Ltag, ttag _Rtag, typename _Rhs, typename> inline
auto Base_<_Derived, _Traits, _Type>::inplace_mul(const _Rhs& _Right) {
	Matrix_<value_type, _Myt::_ctrs, _Rhs::_ctcs> _Ret(rows(), _Right.cols());
	detail::_Blas_kernel_impl<value_type>::mul<_Ltag, _Rtag>(this->plvt(), _Right.plvt(), _Ret.plvt());
	return std::move(_Ret);
}
_TYPES_END
DGE_MATRICE_END