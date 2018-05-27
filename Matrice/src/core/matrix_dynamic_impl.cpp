#include "../../include/Matrice/core/matrix.h"

MATRICE_NAMESPACE_BEGIN_TYPES
template<typename _Ty>
void Matrix_<_Ty, 0, 0>::create(int_t rows, int_t cols)
{
	base_t::m_storage = std::move(details::Storage_<value_t>::Allocator<0, 0, -1>(rows, cols));
	m_data = base_t::m_storage.data(), m_rows = rows, m_cols = cols;
#ifdef __use_ocv_as_view__
	//*static_cast<ocv_view_t*>(this) = ocv_view_t(m_rows, m_cols, ocv_view_t_cast<value_t>::type, m_data);
#endif
}
template void Matrix_<int, 0, 0>::create(int_t, int_t);
template void Matrix_<char, 0, 0>::create(int_t, int_t);
template void Matrix_<bool, 0, 0>::create(int_t, int_t);
template void Matrix_<float, 0, 0>::create(int_t, int_t);
template void Matrix_<double, 0, 0>::create(int_t, int_t);
template void Matrix_<unsigned char, 0, 0>::create(int_t, int_t);

MATRICE_NAMESPACE_END_TYPES