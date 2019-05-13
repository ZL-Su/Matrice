
#include "../../include/Matrice/util/_exception.h"
#include "../../include/Matrice/private/_unified_memory.h"

#if (defined MATRICE_ENABLE_CUDA && !defined __disable_cuda__)
#include <cuda_runtime.h>
#pragma warning(disable: 4715 4661 4224 4267 4244 4819 4199)

using std::size_t;
using uchar = unsigned char;

#pragma region <!-- unified_sync class implementation -->
MATRICE_PRIVATE_BEGIN
template<typename _Scalar, Location _Host>
_Scalar* unified_sync<_Scalar, _Host, OnDevice, LINEAR>::op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _1)
{
	const auto _Stat = cudaMemcpy(_Dst, _Src, _Rows*_Cols * sizeof(_Scalar), ::cudaMemcpyHostToDevice);
#ifdef _DEBUG
	DGELOM_CHECK(_Stat == cudaSuccess, cudaGetErrorString(_Stat));
#endif
	return (_Dst);
}
template int* unified_sync<int, OnStack, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnStack, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnStack, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnStack, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template uchar* unified_sync<uchar, OnStack, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnStack, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnStack, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnHeap, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnHeap, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnHeap, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnHeap, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template uchar* unified_sync<uchar, OnHeap, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnHeap, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnHeap, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);

template<typename _Scalar, Location _Host>
_Scalar* unified_sync<_Scalar, _Host, OnDevice, PITCHED>::op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _Pytes)
{
	const auto _Stat = cudaMemcpy2D(_Dst, _Pytes, _Src, _Cols * sizeof(_Scalar), _Cols * sizeof(_Scalar), _Rows, ::cudaMemcpyHostToDevice);
#ifdef _DEBUG
	DGELOM_CHECK(_Stat == cudaSuccess, cudaGetErrorString(_Stat));
#endif
	return (_Dst);
}
template uchar* unified_sync<uchar, OnStack, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnStack, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnStack, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnStack, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnStack, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnStack, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnStack, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template uchar* unified_sync<uchar, OnHeap, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnHeap, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnHeap, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnHeap, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnHeap, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnHeap, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnHeap, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);

template<typename _Scalar, Location _Host>
_Scalar* unified_sync<_Scalar, OnDevice, _Host, LINEAR>::op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _1)
{
	auto _Stat = cudaMemcpy(_Dst, _Src, _Rows *_Cols * sizeof(_Scalar), ::cudaMemcpyDeviceToHost);
#ifdef _DEBUG
	DGELOM_CHECK(_Stat == cudaSuccess, cudaGetErrorString(_Stat));
#endif
	return (_Dst);
}
template uchar* unified_sync<uchar, OnDevice, OnStack, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnDevice, OnStack, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnDevice, OnStack, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template uchar* unified_sync<uchar, OnDevice, OnHeap, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnDevice, OnHeap, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnDevice, OnHeap, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnDevice, OnStack, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnDevice, OnStack, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnDevice, OnStack, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnDevice, OnStack, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnDevice, OnHeap, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnDevice, OnHeap, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnDevice, OnHeap, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnDevice, OnHeap, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);

template<typename _Scalar, Location _Host>
_Scalar* unified_sync<_Scalar, OnDevice, _Host, PITCHED>::op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _Pytes)
{
	auto _Stat = cudaMemcpy2D(_Dst, _Cols * sizeof(_Scalar), _Src, _Pytes, _Cols * sizeof(_Scalar), _Rows, ::cudaMemcpyDeviceToHost);
#ifdef _DEBUG
	DGELOM_CHECK(_Stat == cudaSuccess, cudaGetErrorString(_Stat));
#endif
	return (_Dst);
}

template uchar* unified_sync<uchar, OnDevice, OnStack, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnDevice, OnStack, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnDevice, OnStack, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template uchar* unified_sync<uchar, OnDevice, OnHeap, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnDevice, OnHeap, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnDevice, OnHeap, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnDevice, OnStack, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnDevice, OnStack, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnDevice, OnStack, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnDevice, OnStack, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnDevice, OnHeap, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnDevice, OnHeap, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnDevice, OnHeap, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnDevice, OnHeap, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);

template<typename _Scalar>
_Scalar* unified_sync<_Scalar, OnDevice, OnDevice, LINEAR>::op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _1)
{
	auto _Stat = cudaMemcpy(_Dst, _Src, _Rows *_Cols * sizeof(_Scalar), ::cudaMemcpyDeviceToDevice);
#ifdef _DEBUG
	DGELOM_CHECK(_Stat == cudaSuccess, cudaGetErrorString(_Stat));
#endif
	return (_Dst);
}
template uchar* unified_sync<uchar, OnDevice, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnDevice, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnDevice, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnDevice, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnDevice, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnDevice, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnDevice, OnDevice, LINEAR>::op(pointer, const_pointer, size_t, size_t, size_t);

template<typename _Scalar>
_Scalar* unified_sync<_Scalar, OnDevice, OnDevice, PITCHED>::op(pointer _Dst, const_pointer _Src, size_t _Rows, size_t _Cols, size_t _Pytes)
{
	auto _Stat = cudaMemcpy2D(_Dst, _Pytes, _Src, _Pytes, _Cols * sizeof(_Scalar), _Rows, ::cudaMemcpyDeviceToDevice);
#ifdef _DEBUG
	DGELOM_CHECK(_Stat == cudaSuccess, cudaGetErrorString(_Stat));
#endif
	return (_Dst);
}
template uchar* unified_sync<uchar, OnDevice, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template float* unified_sync<float, OnDevice, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template double* unified_sync<double, OnDevice, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template int* unified_sync<int, OnDevice, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template char* unified_sync<char, OnDevice, OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template bool* unified_sync<bool, OnDevice, Location::OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
template size_t* unified_sync<size_t, OnDevice, Location::OnDevice, PITCHED>::op(pointer, const_pointer, size_t, size_t, size_t);
MATRICE_PRIVATE_END
#pragma endregion

#endif