#include "../../include/Matrice/core/matrix.h"

MATRICE_NAMESPACE_BEGIN_TYPES

#define MATRICE_INSTANTIATE_METHOD_CREATE(TYPE, ALTYPE) \
template void Matrix_<TYPE, ALTYPE>::__create_impl(size_t, size_t);

template<typename _Ty>
void Matrix_<_Ty, 0, 0>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	this->_Reset_no_alloc({ rows, cols });
}

MATRICE_INSTANTIATE_METHOD_CREATE(int, ::dynamic)
MATRICE_INSTANTIATE_METHOD_CREATE(bool, ::dynamic)
MATRICE_INSTANTIATE_METHOD_CREATE(char, ::dynamic)
MATRICE_INSTANTIATE_METHOD_CREATE(float, ::dynamic)
MATRICE_INSTANTIATE_METHOD_CREATE(double, ::dynamic)
MATRICE_INSTANTIATE_METHOD_CREATE(uint8_t, ::dynamic)

#ifdef MATRICE_ENABLE_CUDA
template<typename _Ty>
void Matrix_<_Ty, ::global>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	this->_Reset_no_alloc({ rows, cols });
}
template<typename _Ty>
void Matrix_<_Ty, ::device>::__create_impl(size_t rows, size_t cols) {
	_Mybase::m_storage.create(rows, cols);
	this->_Reset_no_alloc({ rows, cols });
}

MATRICE_INSTANTIATE_METHOD_CREATE(int, ::global)
MATRICE_INSTANTIATE_METHOD_CREATE(bool, ::global)
MATRICE_INSTANTIATE_METHOD_CREATE(char, ::global)
MATRICE_INSTANTIATE_METHOD_CREATE(float, ::global)
MATRICE_INSTANTIATE_METHOD_CREATE(double, ::global)
MATRICE_INSTANTIATE_METHOD_CREATE(uint8_t, ::global)
MATRICE_INSTANTIATE_METHOD_CREATE(int, ::device)
MATRICE_INSTANTIATE_METHOD_CREATE(bool, ::device)
MATRICE_INSTANTIATE_METHOD_CREATE(char, ::device)
MATRICE_INSTANTIATE_METHOD_CREATE(float, ::device)
MATRICE_INSTANTIATE_METHOD_CREATE(double, ::device)
MATRICE_INSTANTIATE_METHOD_CREATE(uint8_t, ::device)
#endif

MATRICE_NAMESPACE_END_TYPES