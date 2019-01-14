#include "../../include/Matrice/core/matrix.h"

MATRICE_NAMESPACE_BEGIN_TYPES

#define MATRICE_INSTANTIATE_METHOD_CREATE(TYPE) \
template void Matrix_<TYPE, 0,0>::__create_impl(diff_t, diff_t);

template<typename _Ty>
void Matrix_<_Ty, 0, 0>::__create_impl(diff_t rows, diff_t cols) {
	_Mybase::m_storage = std::move(detail::Storage_<value_t>::Allocator<0, 0>(rows, cols));
	_Mybase::m_rows = rows, _Mybase::m_cols = cols;
	_Mybase::m_data = internal::_Proxy_checked(_Mybase::m_storage.data());
	_Mybase::base_t::_Flush_view_buf();
}
MATRICE_INSTANTIATE_METHOD_CREATE(int)
MATRICE_INSTANTIATE_METHOD_CREATE(bool)
MATRICE_INSTANTIATE_METHOD_CREATE(char)
MATRICE_INSTANTIATE_METHOD_CREATE(float)
MATRICE_INSTANTIATE_METHOD_CREATE(double)
MATRICE_INSTANTIATE_METHOD_CREATE(uint8_t)

MATRICE_NAMESPACE_END_TYPES