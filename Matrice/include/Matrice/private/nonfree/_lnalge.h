#pragma once
#include <tuple>
#include "../_type_traits.h"

DGE_MATRICE_BEGIN _DETAIL_BEGIN

template<typename _Ty, 
	typename = std::enable_if_t<is_floating_point_v<_Ty>>> 
struct _Blas_kernel_impl {
	static_assert("Oops, unsupported data type _Ty in _Blas_kernel_impl<_Ty, void>.");

	template<typename... _Args> static constexpr auto dot(const _Args&...) {}
	template<typename... _Args> static constexpr auto mul(const _Args&...) {}
};

template<typename _Ty, 
	typename = std::enable_if_t<is_floating_point_v<_Ty>>> 
struct _Lapack_kernel_impl {
	static_assert("Oops, unsupported data type _Ty in _Lapack_kernel_impl<_Ty, void>.");

	template<typename... _Args> static constexpr auto svd(const _Args&...) {}
	template<typename... _Args> static constexpr auto spd(const _Args&...) {}
	template<typename... _Args> static constexpr auto lud(const _Args&...) {}
	template<typename... _Args> static constexpr auto slv(const _Args&...) {}
};

_DETAIL_END DGE_MATRICE_END

#include "inl\_lnalge.inl"