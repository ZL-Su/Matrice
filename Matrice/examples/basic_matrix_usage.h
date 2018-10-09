/*************************************************************
	This is an example illustrating the use of matrix type in
	the Matrice C++ library.
 *************************************************************/
#pragma once

#include "../include/Matrice/core/matrix.h"

using default_type = dgelom::default_type;

int main() try
{
#pragma region \how to create a matrix
	/* 
	The matrix type has a uniform interface like:
		template<typename T, int _Rows, int _Cols> class Matrix_<T, _Rows, _Cols>.
	Template parameter "T" must be a scalar type, which gives the data type of 
		the matrix; "_Rows" and "_Cols" give the compile time size of the matrix.
		The matrix is created on different memories according to the values of 
		"_Rows" and "_Cols".
	*/

	// \create a 3x3 matrix on the stack without giving initial values
	dgelom::Matrix_<default_type, 3, 3> _Mat33;

	// \create a 3x3 zero-valued matrix on the stack
	dgelom::Matrix_<default_type, 3, 3> _Zero33{ dgelom::zero<default_type>::value };

	// \create a 3x3 identity matrix on the stack
	dgelom::Matrix_<default_type, 3, 3> _Iden33{ 1., 0., 0., 0., 1., 0., 0., 0., 1. };

	// \create a 3x3 matrix from a pointer
	default_type _Ptr[9] = {/*....*/ };
	dgelom::Matrix_<default_type, 3, 3> _Mat33_from_ptr(_Ptr);

	// \create a stack matrix from an existing matrix
	dgelom::Matrix_<default_type, 3, 3> _Mat33_copy_1(_Mat33); //use copy constructuor
	dgelom::Matrix_<default_type, 3, 3> _Mat33_copy_2 = _Mat33;            //use assignment operator

	// \create a dynamic matrix on host
	dgelom::Matrix_<default_type, 0, 0> _Mat_1; //empty matrix
	dgelom::Matrix_<default_type, 0, 0> _Mat_2(3, 3); // 3x3 matrix
	dgelom::Matrix_<default_type, 0, 0> _Mat_2(3, 3, dgelom::zero<default_type>::value); // 3x3 matrix, initialized by zeros

#pragma endregion
}
catch (std::exception& _E) { throw _E; }