#include "../../include/Matrice/core/matrix.h"

MATRICE_NAMESPACE_BEGIN_TYPES
template<typename _Ty>
void Matrix_<_Ty, 0, 0>::create(int_t rows, int_t cols)
{
	base_t::m_storage = std::move(details::Storage_<value_t>::Allocator<0, 0, allocator_trait<0, 0>::value>(rows, cols));
	m_rows = rows, m_cols = cols;
	m_data = _Proxy_checked(base_t::m_storage.data());
	base_t::base_t::_Flush_view_buf();
}
template void Matrix_<int, 0, 0>::create(int_t, int_t);
template void Matrix_<char, 0, 0>::create(int_t, int_t);
template void Matrix_<bool, 0, 0>::create(int_t, int_t);
template void Matrix_<float, 0, 0>::create(int_t, int_t);
template void Matrix_<double, 0, 0>::create(int_t, int_t);
template void Matrix_<unsigned char, 0, 0>::create(int_t, int_t);
template<typename _Ty>
void Matrix_<_Ty, 0, 0>::create(int_t rows, int_t cols, value_t _val)
{
	base_t::m_storage = std::move(details::Storage_<value_t>::Allocator<0, 0, allocator_trait<0,0>::value>(rows, cols));
	base_t::m_storage = { _val };
	m_rows = rows, m_cols = cols;
	m_data = _Proxy_checked(base_t::m_storage.data());
	base_t::base_t::_Flush_view_buf();
}
template void Matrix_<int, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<char, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<bool, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<float, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<double, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<unsigned char, 0, 0>::create(int_t, int_t, value_t);
MATRICE_NAMESPACE_END_TYPES