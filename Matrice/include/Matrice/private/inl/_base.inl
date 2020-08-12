/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include "../_plain_base.hpp"
#include "../_range.h"
#include "../../thread/_thread.h"
#include "private/nonfree/blas_lapack_kernel.h"
#include "private/math/kernel_wrapper.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs> MATRICE_HOST_INL
auto Base_<_Derived, _Traits, _Type>::sub_inplace(const _Rhs& _Right) {
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
	for (auto _Idx = 0; _Idx < this->size(); ++_Idx) 
		m_data[_Idx] -= _Right(_Idx);
#endif
	return forward<_Myt>(*this);
}

template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs> MATRICE_GLOBAL_INL
auto Base_<_Derived, _Traits, _Type>::mv_inplace(const _Rhs& _Right)const {
	Matrix_<value_type, _Myt::rows_at_compiletime, 1> _Ret(m_rows);
	return forward<decltype(_Ret)>(blas_kernel_t::gemv(*this, _Right, _Ret));
}

template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs> MATRICE_HOST_INL
auto Base_<_Derived, _Traits, _Type>::mul_inplace(const _Rhs& _Right) const{
	Matrix_<value_type, _Myt::rows_at_compiletime, _Rhs::cols_at_compiletime>_Ret(m_rows, _Right.cols());
	return forward<decltype(_Ret)>(blas_kernel_t::gemm(*this, _Right, _Ret));
}

template<typename _Derived, typename _Traits, typename _Type>
template<int Rows, int Cols> MATRICE_GLOBAL_INL 
_Type Base_<_Derived, _Traits, _Type>::contract(const Matrix_<_Type, Rows, Cols>& other) const
{
#ifdef MATRICE_DEBUG
	DGELOM_CHECK(rows() == other.rows() && cols() == other.cols(),
		"Rows and columns must be matched in Base_<_Derived, _Traits, _Type>::contract(_Rhs).");
#endif
	return (this->const_derived()* other).sum();
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

/**
 *\brief Factory function returns the deleter of a given dgelom::Matrix_ instance.
 *\param [_M] instance of dgelom::Matrix_;
 */
template<typename _Mty>
MATRICE_HOST_INL decltype(auto) make_matrix_deleter(const _Mty& _M)noexcept {
	return _M.deleter();
};

/**
 *\brief Factory function to make a dgelom::Matrix_ with given parameters.
 *\param [params...] variadic argument supports any inputs supported by the ctors of dgelom::Matrix_;
 */
template<typename _Ty, int _Rows, int _Cols, typename... _Args>
MATRICE_GLOBAL_INL decltype(auto) make_matrix(_Args&&... params) {
	return detail::Matrix_<_Ty, _Rows, _Cols>(forward<_Args>(params)...);
};

/**
 *\brief Make the given argument data := zero.
 *\param [data] can be a scalar, matrix or any types with the member data() and operator[].
 */
template<typename _Ty>
MATRICE_GLOBAL_INL remove_all_t<_Ty>& make_zero(_Ty& data) noexcept {
	using trivial_type = remove_all_t<_Ty>;

	if constexpr (is_scalar_v<trivial_type>) {
		data = zero<trivial_type>;
	}
	else if constexpr (is_matrix_v<trivial_type>) {
		data = zero<typename trivial_type::value_type>;
	}
	else if constexpr (has_data_v<_Ty>) {
		auto dp = data.data();
		for (const auto& idx : range(0, data.size())) {
			dp[idx] = zero<typename trivial_type::value_type>;
		}
	}
	else {
		// do nothing
	}
	return (data);
}

/**
 * \brief Make a copy from the given matrix _M
 */
template<typename _Mty, MATRICE_ENABLE_IF(is_matrix_v<_Mty>)>
MATRICE_HOST_INL _Mty copy(const _Mty& _M) {
	return { _M };
}

/**
 * \brief Make a full view to the given matrix or vector _M.
 */
template<typename _Mty, 
	MATRICE_ENABLE_IF(is_matrix_v<_Mty> || is_fxdvector_v<_Mty>)>
MATRICE_GLOBAL_FINL auto view(_Mty& _M) noexcept {
	return detail::_Matrix_block<typename _Mty::value_type>(_M.data(), _M.cols(), _M.rows());
}

/**
 * \brief Swap matrice _L and _R
 */
template<typename _Mty>
MATRICE_HOST_INL void swap(_Mty& _L, _Mty& _R) noexcept
{
	auto _Tmp = move(_L);
	_L = move(_R);
	_R = move(_Tmp);
}
DGE_MATRICE_END