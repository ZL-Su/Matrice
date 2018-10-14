/*************************************************************
	This is an example illustrating the use of matrix type in
	the Matrice C++ library.
 *************************************************************/
#pragma once

#include "../include/Matrice/core/matrix.h"

using default_type = dgelom::default_type;

int main() try
{
#pragma region <!-- How to create a matrix -->
	{
	/* 
	The matrix type has a uniform interface as follows:
		template<typename T, int _Rows, int _Cols> class Matrix_<T, _Rows, _Cols>.
	Template parameter "T" must be a scalar type, which gives the data type of 
		the matrix; "_Rows" and "_Cols" give the compile time size of the matrix.
		The matrix is created on different memories according to the values of 
		"_Rows" and "_Cols".
	*/

	// \create a 3x3 matrix on the stack by giving _Rows > 0 and _Cols > 0
	dgelom::Matrix_<default_type, 3, 3> _Mat33;

	// \create a 3x3 zero-valued matrix on the stack
	dgelom::Matrix_<default_type, 3, 3> _Zero33{ dgelom::zero<default_type>::value };

	// \create a 3x3 identity matrix on the stack
	dgelom::Matrix_<default_type, 3, 3> _Iden33{ 1., 0., 0., 0., 1., 0., 0., 0., 1. };

	// \create a 3x3 matrix from a pointer, one example:
	default_type _Ptr[9] = {/*....*/ };
	dgelom::Matrix_<default_type, 3, 3> _Mat33_from_ptr(_Ptr);

	// \create a stack matrix from an existing matrix
	dgelom::Matrix_<default_type, 3, 3> _Mat33_copy_1(_Mat33); //use copy constructuor
	dgelom::Matrix_<default_type, 3, 3> _Mat33_copy_2 = _Mat33;            //use assignment operator

	// \create a dynamic matrix on host by specifying _Rows = 0 and _Cols = 0
	dgelom::Matrix_<default_type, 0, 0> _Mat_1; //empty matrix
	dgelom::Matrix_<default_type, 0, 0> _Mat_2(3, 3); // 3x3 matrix
	dgelom::Matrix_<default_type, 0, 0> _Mat_3(3, 3, dgelom::zero<default_type>::value); // 3x3 matrix, initialized by zeros
	dgelom::Matrix_<default_type, 0, 0> _Mat_4(3, 3, _Ptr); // from a pointer
	// ...or, use the more compact alias type
	dgelom::Matrix<default_type> _Mat_5/*(...)*/;

	// \even, we can create matrix like this
	dgelom::Matrix_<default_type, 3, 3> _Mat33_from_dyn(_Mat_3);
	dgelom::Matrix_<default_type, 0, 0> _Mat_from_sta(_Zero33);
	// ...or, the matrix can be created like this
	dgelom::Matrix_<default_type, 3, 3> _Mat33_from_dyn_1 = _Mat_3;
	dgelom::Matrix_<default_type, 0, 0> _Mat_from_sta_1 = _Zero33;
	}
#pragma endregion

#pragma region <!-- Basic matrix operations -->
	/*
	What you should know before using matrix operations:
		(1) All operations are almost applied to all matrix types,
		(2) Almost all matrix operations are implemented with the lazy 
			 evaluation technology, aiming to avoid the unnecessary cell
			 computation and data copy,
		(3) Be careful to use the 'auto' keyword for matrix operations 
			 except that you absolutely know what you are doing.
	*/
	{
		// Given two matrices A and B:
		using matrix_t = dgelom::Matrix_<default_type, 3, 3>;
		matrix_t A{ 1., 0., 0., 0., 1., 0., 0., 0., 1. }, B{ 0.25 };

		// \Element-wise operations: +, -, *, /
		matrix_t A_plus_B = A + B;
		matrix_t A_minus_B = A - B;
		matrix_t A_dotm_B = A * B;
		matrix_t A_over_B = A / B;

		// ...or, we just forming expressions but not executing the operations:
		auto A_plus_B_exp = A + B;
		auto A_minus_B_exp = A - B;
		auto A_dotm_B_exp = A * B;
		auto A_over_B_exp = A / B;

		// Of cource, Matrice lib has supported several ways to evaluate an expression,
		// ...invoke evaluation by passing a expression to a matrix:
		matrix_t Ewise_sum_1(A_plus_B_exp); 
		// or...
		matrix_t Ewise_sum_2 = A_plus_B_exp; 

		// or use the 'eval()' method of the expression...
		auto Ewise_sum_3 = A_plus_B_exp.eval(); 
		// in the above way, the returned matrix type is inferred automatically, which generally is the optimal matrix type. If you wanna use the return type that you specified, the templated version of 'eval()' method should be called as follows
		auto Ewise_sum_4 = A_plus_B_exp.eval<dgelom::Matrix<default_type>>();

		// \Metrix-level operations: matmul, transpose, inverse...
		// ...expresstion
		auto A_mul_B_exp = A.mul(B);
		// ...direct evaluation
		matrix_t A_mul_B = A.mul(B);
		// ...or matrix_t A_mul_B = A_mul_B_exp;
		// ...or auto A_mul_B = A_mul_B_exp.eval();
	}
#pragma endregion
}
catch (std::exception& _E) { throw _E; }