#include "../../include/Matrice/core/matrix.h"

namespace dgelom {
namespace types {

#pragma region Special Allocator
template<typename _Ty, class _Storage>
MatrixBase_<_Ty, _Storage>::MatrixBase_(const int_t _m, const int_t _n) noexcept
	: m_storage(_m, _n), m_rows(m_storage.rows()), m_cols(m_storage.cols())
{
	m_data = m_storage.data();
#ifdef __use_ocv_as_view__
	*static_cast<ocv_view_t*>(this) = ocv_view_t(_m, _n, ocv_view_t_cast<value_t>::type, m_data);
#endif // __use_ocv_as_view__
};

template class MatrixBase_<float, details::Storage_<float>::SharedAllocator<>>;
template class MatrixBase_<double, details::Storage_<double>::SharedAllocator<>>;
template class MatrixBase_<float, details::Storage_<float>::DynamicAllocator<>>;
template class MatrixBase_<double, details::Storage_<double>::DynamicAllocator<>>;

template class MatrixBase_<float, details::Storage_<float>::ManagedAllocator<2, 2>>;
template class MatrixBase_<double, details::Storage_<double>::ManagedAllocator<2, 2>>;
template class MatrixBase_<float, details::Storage_<float>::ManagedAllocator<2, 1>>;
template class MatrixBase_<double, details::Storage_<double>::ManagedAllocator<2, 1>>;
template class MatrixBase_<float, details::Storage_<float>::ManagedAllocator<1, 2>>;
template class MatrixBase_<double, details::Storage_<double>::ManagedAllocator<1, 2>>;
template class MatrixBase_<float, details::Storage_<float>::ManagedAllocator<3, 3>>;
template class MatrixBase_<double, details::Storage_<double>::ManagedAllocator<3, 3>>;
template class MatrixBase_<float, details::Storage_<float>::ManagedAllocator<3, 1>>;
template class MatrixBase_<double, details::Storage_<double>::ManagedAllocator<3, 1>>;
template class MatrixBase_<float, details::Storage_<float>::ManagedAllocator<1, 3>>;
template class MatrixBase_<double, details::Storage_<double>::ManagedAllocator<1, 3>>;
template class MatrixBase_<float, details::Storage_<float>::ManagedAllocator<3, 4>>;
template class MatrixBase_<double, details::Storage_<double>::ManagedAllocator<3, 4>>;
template class MatrixBase_<float, details::Storage_<float>::ManagedAllocator<4, 4>>;
template class MatrixBase_<double, details::Storage_<double>::ManagedAllocator<4, 4>>;
#pragma endregion
}}