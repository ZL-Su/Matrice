/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once

#include "private/_plain_base.hpp"
#include "private/_range.h"
#include "private/nonfree/blas_lapack_kernel.h"
#include "private/math/kernel_wrapper.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
template<typename _Derived, typename _Traits, typename _Type>
template<typename _Rhs> MATRICE_HOST_INL
decltype(auto) Base_<_Derived, _Traits, _Type>::sub_inplace(const _Rhs& _Right) {
#ifdef MATRICE_SIMD_ARCH_
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
		DGELOM_CHECK(this->shape()==_Right.shape(), "Inconsistent shapes.");
		const auto rdata = _Right.data();
		for (const auto i : range(0, plen, size_t(step))) {
			(packet_t(m_data+i)-packet_t(rdata+i)).unpack(m_data+i);
		}
		for (auto i = plen; i < this->size(); ++i) {
			m_data[i] -= rdata[i];
		}
	}
	else {
		throw("No implemented error");
	}
#else
	for (auto _Idx = 0; _Idx < this->size(); ++_Idx)
		m_data[_Idx] -= _Right(_Idx);
#endif
	return this->derived();
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
template<index_t Rows, index_t Cols> MATRICE_GLOBAL_INL
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
 * \brief Make a copy from the given matrix _M.
 */
template<typename _Mty, typename>
MATRICE_HOST_INL _Mty copy(const _Mty& _M) {
	return { _M };
}

/**
 * \brief Swap matrice _L and _R.
 */
template<typename _Mty>
MATRICE_HOST_INL void swap(_Mty& _L, _Mty& _R) noexcept
{
	auto _Tmp = move(_L);
	_L = move(_R);
	_R = move(_Tmp);
}

/**
 *\func dgelom::stack<_Axis>(initlist<const _Mty&>)
 *\brief Stack a matrix or vector list to make a new matrix.
 *\param list Argument holds the given matrix or vector list with std::initializer_list
 */
template<axis _Axis, typename _Mty, typename>
MATRICE_HOST_FINL auto stack(const initlist<_Mty>& list)
{
	// Stack along the row dimension
	if constexpr(_Axis == axis::y) {
		static constexpr auto ctcols = _Mty::cols_at_compiletime;
		const auto rows = list.size() * list.begin()->rows();
		const auto cols = list.begin()->cols();
		auto _Ret = detail::Matrix_<_Mty::value_type, ::dynamic, ctcols>(rows, cols);
		auto off = size_t(0);
		for (auto _It = list.begin(); _It != list.end(); ++_It) {
			_Ret.block<0>(off, off+_It->rows()) = *_It;
			off += _It->rows();
		}
		return _Ret;
	}
	// Stack along the column dimension
	if constexpr(_Axis == axis::x) {
		static constexpr auto ctrows = _Mty::rows_at_compiletime;
		const auto cols = list.size() * list.begin()->cols();
		const auto rows = list.begin()->rows();
		auto _Ret = detail::Matrix_<_Mty::value_type, ctrows, ::dynamic>(rows, cols);
		auto off = size_t(0);
		for (auto _It = list.begin(); _It != list.end(); ++_It) {
			_Ret.block<1>(off, off+_It->cols()) = *_It;
			off += _It->cols();
		}
		return _Ret;
	}
}
DGE_MATRICE_END