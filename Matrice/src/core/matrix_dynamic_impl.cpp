#include "../../include/Matrice/core/matrix.h"

MATRICE_NAMESPACE_BEGIN_TYPES

#define MATRICE_INSTANTIATE_METHOD_CREATE(TYPE) \
template void Matrix_<TYPE, 0, 0>::__create_impl(size_t, size_t);

template<typename _Ty>
void Matrix_<_Ty, 0, 0>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	this->_Reset_no_alloc({ rows, cols });
}

MATRICE_INSTANTIATE_METHOD_CREATE(int)
MATRICE_INSTANTIATE_METHOD_CREATE(bool)
MATRICE_INSTANTIATE_METHOD_CREATE(char)
MATRICE_INSTANTIATE_METHOD_CREATE(float)
MATRICE_INSTANTIATE_METHOD_CREATE(double)
MATRICE_INSTANTIATE_METHOD_CREATE(uint8_t)

#ifdef MATRICE_ENABLE_CUDA
template<typename _Ty>
void Matrix_<_Ty, -1, 0>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	this->_Reset_no_alloc({ rows, cols });
}
template<typename _Ty>
void Matrix_<_Ty, -1, -1>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	this->_Reset_no_alloc({ rows, cols });
}

template void Matrix_<int, -1, 0>::__create_impl(size_t, size_t);
template void Matrix_<int, -1, -1>::__create_impl(size_t, size_t);
template void Matrix_<bool, -1, 0>::__create_impl(size_t, size_t);
template void Matrix_<bool, -1, -1>::__create_impl(size_t, size_t);
template void Matrix_<char, -1, 0>::__create_impl(size_t, size_t);
template void Matrix_<char, -1, -1>::__create_impl(size_t, size_t);
template void Matrix_<float, -1, 0>::__create_impl(size_t, size_t);
template void Matrix_<float, -1, -1>::__create_impl(size_t, size_t);
template void Matrix_<double, -1, 0>::__create_impl(size_t, size_t);
template void Matrix_<double, -1, -1>::__create_impl(size_t, size_t);
template void Matrix_<uint8_t, -1, 0>::__create_impl(size_t, size_t);
template void Matrix_<uint8_t, -1, -1>::__create_impl(size_t, size_t);
#endif

MATRICE_NAMESPACE_END_TYPES