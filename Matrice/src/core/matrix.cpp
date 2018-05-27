#include "..\..\include\core\matrix.h"

namespace Dgelo {
namespace types {

template<typename _Ty, class _Storage>
MatrixBase_<_Ty, _Storage>::MatrixBase_(const int_t _m, const int_t _n) noexcept
	: m_storage(_m, _n)
{
	m_ptr = m_storage.ptr();
#ifdef __use_container_view__
	*static_cast<View*>(this) = View(_m, _n, viewty_cast<value_t>::type, m_ptr);
#endif // __use_container_view__
};
template<typename _Ty, class _Storage>
MatrixBase_<_Ty, _Storage>::MatrixBase_(const std::size_t _m, const std::size_t _n) noexcept
	: m_storage(_m, _n)
{
	m_ptr = m_storage.ptr();
#ifdef __use_container_view__
	*static_cast<View*>(this) = View(_m, _n, viewty_cast<value_t>::type, m_ptr);
#endif // __use_container_view__
};


template MatrixBase_<float, details::Storage_<float>::SharedAllocator<>>::MatrixBase_(const int_t, const int_t) noexcept;
template MatrixBase_<float, details::Storage_<float>::DynamicAllocator<>>::MatrixBase_(const int_t, const int_t) noexcept;
}}