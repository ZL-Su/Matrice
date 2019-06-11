#include "../../include/Matrice/core/matrix.h"

MATRICE_NAMESPACE_BEGIN_TYPES

#define MATRICE_INSTANTIATE_METHOD_CREATE(TYPE) \
template void Matrix_<TYPE, 0, 0>::__create_impl(size_t, size_t);\
template void Matrix_<TYPE,-1, 0>::__create_impl(size_t, size_t);\
template void Matrix_<TYPE,-1,-1>::__create_impl(size_t, size_t);

template<typename _Ty>
void Matrix_<_Ty, 0, 0>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	_Mybase::m_rows = rows, _Mybase::m_cols = cols;
	_Mybase::_Myshape = { 1,1,rows,cols };
	_Mybase::m_data = internal::_Proxy_checked(_Mybase::m_storage.data());
	_Mybase::base_t::_Flush_view_buf();
}
template<typename _Ty>
void Matrix_<_Ty, -1, 0>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	_Mybase::m_rows = rows, _Mybase::m_cols = cols;
	_Mybase::_Myshape = { 1,1,rows,cols };
	_Mybase::m_data = internal::_Proxy_checked(_Mybase::m_storage.data());
	_Mybase::base_t::_Flush_view_buf();
}
template<typename _Ty>
void Matrix_<_Ty, -1, -1>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	_Mybase::m_rows = rows, _Mybase::m_cols = cols;
	_Mybase::_Myshape = { 1,1,rows,cols };
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