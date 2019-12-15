/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../_plain_base.hpp"
#include "../_range.h"
#include "../../thread/_thread.h"
#include "private/nonfree/_lnalge.h"
#include "private/math/kernel_wrapper.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs> MATRICE_HOST_INL
decltype(auto) Base_<_Derived, _Traits, _Type>::sub_inplace(const _Rhs& _Right) {
#ifdef MATRICE_SIMD_ARCH
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
#else
	for (auto _Idx = 0; _Idx < this->size(); ++_Idx) m_data[_Idx] -= _Right(_Idx);
#endif
	return forward<_Myt>(*this);
}

template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs, typename> inline
decltype(auto) Base_<_Derived, _Traits, _Type>::mul_(const _Rhs& _Right) {
#ifdef MATRICE_SIMD_ARCH
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
#else

#endif
	return forward<_Rhs>(_Right);
}

template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs> MATRICE_HOST_INL
auto Base_<_Derived, _Traits, _Type>::mul_inplace(const _Rhs& _Right) const{
	Matrix_<value_type, _Myt::rows_at_compiletime, _Rhs::cols_at_compiletime>_Ret(m_rows, _Right.cols());
	return forward<decltype(_Ret)>(blas_kernel_t::gemm(*this, _Right, _Ret));
}

template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs> MATRICE_GLOBAL_INL
_Rhs Base_<_Derived, _Traits, _Type>::spreadmul(const _Rhs& _Right)const {
	_Rhs _Ret(_Right.shape());

	// spread each entry along row and element-wisely mul. with _Right
	if (m_cols == 1 || size() == _Right.rows()) {
		for (const auto _r : range(0, _Ret.rows())) {
			_Ret.rview(_r) = m_data[_r] * _Right.rview(_r);
		}
	}
	// spread each entry along column and element-wisely mul. with _Right
	else if (m_rows == 1 || size() == _Right.cols()) {
		for (const auto _c : range(0, _Ret.cols())) {
			_Ret.cview(_c) = m_data[_c] * _Right.cview(_c);
		}
	}
	else {
		DGELOM_ERROR("Only one-dimension array spread is supported.");
	}
	return forward<_Rhs>(_Ret);
}
_DETAIL_END
template<typename _Mty>
MATRICE_HOST_INL decltype(auto) make_matrix_deleter(const _Mty& _M) noexcept {
	return _M.deleter();
};

template<typename _Ty, int _Rows, int _Cols, typename... _Args>
MATRICE_GLOBAL_INL decltype(auto) make_matrix(_Args&&... params) {
	return detail::Matrix_<_Ty, _Rows, _Cols>(forward<_Args>(params)...);
};

template<typename _Ty>
MATRICE_GLOBAL_INL remove_all_t<_Ty>& make_zero(_Ty& data) noexcept {
	using trivial_type = remove_all_t<_Ty>;

	if constexpr (is_scalar_v<trivial_type>) {
		data = zero<trivial_type>;
	}
	else if constexpr (is_matrix_v<trivial_type>) {
		data = zero<typename trivial_type::value_type>;
	}
	else {
		auto dp = data.data();
		for (const auto& idx : range(0, data.size())) {
			dp[idx] = zero<typename trivial_type::value_type>;
		}
	}

	return (data);
}
DGE_MATRICE_END