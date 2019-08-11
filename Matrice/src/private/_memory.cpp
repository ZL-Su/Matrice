#include <memory>
#include <complex>
#include <util/_exception.h>
#include <private/_memory.h>
#ifdef MATRICE_ENABLE_CUDA
#include <cuda_runtime.h>
#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)
#endif

DGE_MATRICE_BEGIN
namespace internal { namespace impl {
#ifdef MATRICE_ENABLE_CUDA
template<typename _Ty>
MATRICE_GLOBAL _Ty* _Malloc(size_t rows, size_t& cols, device_alloc_tag) {
	_Ty* Data = nullptr;
	::cudaError_t Sts;
	if (cols == 1) {
		Sts = ::cudaMalloc(&Data, rows*cols * sizeof(_Ty));
	}
	else {
		size_t pitch = 0;
		Sts = cudaMallocPitch(&Data, &pitch, cols * sizeof(_Ty), rows);
		cols = pitch;
	}
#if (defined _DEBUG || defined MATRICE_DEBUG)
	DGELOM_CHECK(Sts == ::cudaSuccess, ::cudaGetErrorString(Sts));
#endif
	return Data;
}

template<typename _Ty>
MATRICE_GLOBAL decltype(auto) _Malloc(size_t rows, size_t cols, global_alloc_tag) {
	_Ty* Data = nullptr;
	const auto Sts = ::cudaMallocManaged(&Data, rows*cols*sizeof(_Ty));
#if (defined _DEBUG || defined MATRICE_DEBUG)
	DGELOM_CHECK(Sts == ::cudaSuccess, ::cudaGetErrorString(Sts));
#endif
	return Data;
}
#define MATRICE_CUDA_MEMFUNC_INST(TYPE) \
template TYPE* _Malloc(size_t, size_t&, device_alloc_tag); \
template TYPE* _Malloc(size_t, size_t, global_alloc_tag);

MATRICE_CUDA_MEMFUNC_INST(int)
MATRICE_CUDA_MEMFUNC_INST(char)
MATRICE_CUDA_MEMFUNC_INST(bool)
MATRICE_CUDA_MEMFUNC_INST(float)
MATRICE_CUDA_MEMFUNC_INST(double)
MATRICE_CUDA_MEMFUNC_INST(uint8_t)
MATRICE_CUDA_MEMFUNC_INST(size_t)

#undef MATRICE_CUDA_MEMFUNC_INST

#endif //!MATRICE_ENABLE_CUDA

}}
DGE_MATRICE_END

namespace dgelom { namespace privt {

template<typename _Ty, typename _Ity>
_Ty* aligned_malloc(_Ity size) {
	try {
		auto raw_ptr = std::malloc(size*sizeof(_Ty)+MATRICE_ALIGN_BYTES);
		auto space = reinterpret_cast<size_t>(raw_ptr);
		space = space &~(size_t(MATRICE_ALIGN_BYTES - 1));
		auto aligned_ptr = reinterpret_cast<void*>(space+MATRICE_ALIGN_BYTES);
		*(reinterpret_cast<void**>(aligned_ptr) - 1) = raw_ptr;

		return (reinterpret_cast<_Ty*>(aligned_ptr));
	}
	catch (std::bad_alloc) {
		std::exception("Bad memory allocation.");
	};
}
template<typename _Ty>
void aligned_free(_Ty* aligned_ptr) noexcept {
	if (aligned_ptr) {
		std::free(*(reinterpret_cast<void**>(reinterpret_cast<void*>(aligned_ptr)) - 1));
		aligned_ptr = nullptr;
	}
}
template<typename _Ty>
bool is_aligned(_Ty* aligned_ptr) noexcept {
	return !(reinterpret_cast<size_t>(reinterpret_cast<void*>(aligned_ptr)) % MATRICE_ALIGN_BYTES);
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
template size_t* aligned_malloc<size_t>(int);
template size_t* aligned_malloc<size_t>(std::size_t);
template size_t* aligned_malloc<size_t>(std::ptrdiff_t);

template void aligned_free(int*) noexcept;
template void aligned_free(char*) noexcept;
template void aligned_free(bool*) noexcept;
template void aligned_free(float*) noexcept;
template void aligned_free(double*) noexcept;
template void aligned_free(unsigned char*) noexcept;
template void aligned_free(std::size_t*) noexcept;
template void aligned_free(std::ptrdiff_t*) noexcept;
template bool is_aligned(int*) noexcept;
template bool is_aligned(float*) noexcept;
template bool is_aligned(double*) noexcept;
template bool is_aligned(unsigned char*) noexcept;
template bool is_aligned(std::size_t*) noexcept;
template bool is_aligned(std::ptrdiff_t*) noexcept;
}}