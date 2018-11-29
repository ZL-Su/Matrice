#include "../../include/Matrice/core/matrix.h"

MATRICE_NAMESPACE_BEGIN_TYPES
template<typename _Ty>
void Matrix_<_Ty, 0, 0>::create(int_t rows, int_t cols)
{
	_Mybase::m_storage = std::move(detail::Storage_<value_t>::Allocator<0, 0>(rows, cols));
	_Mybase::m_rows = rows, _Mybase::m_cols = cols;
	_Mybase::m_data = internal::_Proxy_checked(_Mybase::m_storage.data());
	_Mybase::base_t::_Flush_view_buf();
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
	_Mybase::m_storage = std::move(detail::Storage_<value_t>::Allocator<0, 0>(rows, cols));
	_Mybase::m_storage = { _val };
	_Mybase::m_rows = rows, _Mybase::m_cols = cols;
	_Mybase::m_data = internal::_Proxy_checked(_Mybase::m_storage.data());
	_Mybase::base_t::_Flush_view_buf();
}
template void Matrix_<int, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<char, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<bool, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<float, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<double, 0, 0>::create(int_t, int_t, value_t);
template void Matrix_<unsigned char, 0, 0>::create(int_t, int_t, value_t);
MATRICE_NAMESPACE_END_TYPES