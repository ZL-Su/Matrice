#include <memory>
#include <complex>
#include <stdexcept>
#include "../../include/Matrice/private/_memory.h"

namespace dgelom { namespace privt {

template<typename ValueType, typename IntegerType>
ValueType * aligned_malloc(IntegerType size)
{
	try 
	{
		void* raw_ptr = std::malloc(size * sizeof(ValueType) + MATRICE_ALIGN_BYTES);
		std::size_t space = reinterpret_cast<std::size_t>(raw_ptr);
		space = space &~(std::size_t(MATRICE_ALIGN_BYTES - 1));
		void* aligned_ptr = reinterpret_cast<void*>(space + MATRICE_ALIGN_BYTES);
		*(reinterpret_cast<void**>(aligned_ptr) - 1) = raw_ptr;

		return (reinterpret_cast<ValueType*>(aligned_ptr));
	}
	catch (std::bad_alloc)
	{
		std::exception("Bad memory allocation.");
	};
}
template<typename ValueType>
void aligned_free(ValueType * aligned_ptr) noexcept
{
	if (aligned_ptr)
		std::free(*(reinterpret_cast<void**>(reinterpret_cast<void*>(aligned_ptr)) - 1));
}
template<typename ValueType> 
bool is_aligned(ValueType * aligned_ptr) noexcept
{
	return !(reinterpret_cast<std::size_t>(reinterpret_cast<void*>(aligned_ptr)) % MATRICE_ALIGN_BYTES);
}

template int* aligned_malloc<int>(int);
template int* aligned_malloc<int>(std::size_t);
template int* aligned_malloc<int>(std::ptrdiff_t);
template char* aligned_malloc<char>(int);
template char* aligned_malloc<char>(std::size_t);
template char* aligned_malloc<char>(std::ptrdiff_t);
template bool* aligned_malloc<bool>(int);
template bool* aligned_malloc<bool>(std::size_t);
template bool* aligned_malloc<bool>(std::ptrdiff_t);
template float* aligned_malloc<float>(int);
template float* aligned_malloc<float>(std::size_t);
template float* aligned_malloc<float>(std::ptrdiff_t);
template double* aligned_malloc<double>(int);
template double* aligned_malloc<double>(std::size_t);
template double* aligned_malloc<double>(std::ptrdiff_t);
template unsigned char* aligned_malloc<unsigned char>(int);
template unsigned char* aligned_malloc<unsigned char>(std::size_t);
template unsigned char* aligned_malloc<unsigned char>(std::ptrdiff_t);
//template std::complex<float>* aligned_malloc<std::complex<float>>(int);
//template std::complex<float>* aligned_malloc<std::complex<float>>(std::size_t);
//template std::complex<float>* aligned_malloc<std::complex<float>>(std::ptrdiff_t);
//template std::complex<double>* aligned_malloc<std::complex<double>>(int);
//template std::complex<double>* aligned_malloc<std::complex<double>>(std::size_t);
//template std::complex<double>* aligned_malloc<std::complex<double>>(std::ptrdiff_t);

template void aligned_free(int*) noexcept;
template void aligned_free(char*) noexcept;
template void aligned_free(bool*) noexcept;
template void aligned_free(float*) noexcept;
template void aligned_free(double*) noexcept;
template void aligned_free(unsigned char*) noexcept;
template bool is_aligned(int*) noexcept;
template bool is_aligned(float*) noexcept;
template bool is_aligned(double*) noexcept;
template bool is_aligned(unsigned char*) noexcept;
}}