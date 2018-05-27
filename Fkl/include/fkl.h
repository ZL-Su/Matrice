/**************************************************************************
This file is C/C++ interface of FKL, an effcient Modern Fortran library.
Copyright(C) 2015-2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

#ifndef __FKL__
#define __FKL__
#endif

/*********************************************************
 *DESCRIP: Fortran kernel library for numeric computation
 *INCLUDE: 
		1. bla -- basic linear algebra operation
		2. img -- digital image processing
		3. spe -- special functions
		4. dms -- data modeling subroutines
		5. cgs -- computation geometry subroutines
 *Copyright(C) 2015-2018, Zhilong(Dgelom) Su
 *********************************************************
 *********************************************************
 * *Variable Name Convention:
		@1. pointer variable with a underscore prefix "_" is input-only, and modifiable variable has no such prefix.
		@2. 
 * *Function Name Convention:
		function is named in format: _<Type><name>[s|v|m].
		@<Type> is the data type of varaibles such as "s" and "d" for single and double precision floating points respectively. 
		@<Name> is functio name, which represents full or abbreviated name of algorithm or expression that the function to do. For each abbreviated name, its function comment gives detail description.
		@[s|v|m] is optional suffix, to indicate which type is processed, such as "s" for scalar, "v" and "m" for rank-one and rank-two arrays respectively. If no such suffix means the function is suitable for both rank-1 and rank-2 arrays.
 * *Example: somethings like this
		_ddotv is a function which takes two double typed rank-1 arrays as input to calculate their dot product 
***********************************************************/
namespace fkl
{
#pragma region DEFINED TYPES of FKL LIBRARY
	typedef double  Dfp, dp_t;       //Double f.point
	typedef float   Sfp, sp_t;       //Single f.point
	typedef double* dptr;            //Double f.point pointer
	typedef float*  sptr;            //Single f.point pointer
	typedef unsigned char* ucptr;
	typedef int* iptr;
#ifndef DFP
#define DFP Dfp
#endif
#ifndef SFP
#define SFP Sfp
#endif
#pragma endregion

#ifdef __cplusplus
extern "C" {
#endif

namespace bla{
	/*
	 @Purp: determinent of square matrix a [det(a)]
	 @Type: single precision floating point[DFP]
	 @Para:
	   a -- read-only pointer to array
	   n -- dimension of a(n, n)
	*/
	sp_t _sdetm(const sptr _a, int n);
	/*
	 @Purp: determinent of square matrix a [det(a)]
	 @Type: double precision floating point[DFP]
	 @Para:
	   a -- read-only pointer to array
	   n -- dimension of a(n, n)
	*/
	dp_t _ddetm(const dptr _a, int n);

	/*
	 @Purp: dot(inner) product of two vectors
	 @Type: single precision floating point[SFP]
	 @Para:
		_x,_y -- read-only pointer to input vectors
		n -- dimension of vectors _x and _y
	*/
	sp_t _sdotv(const sptr _x, const sptr _y, int n);
	/*
	 @Purp: dot(inner) product of two vectors
	 @Type: double precision floating point[DFP]
	 @Para:
		_x,_y -- read-only pointer to input vectors
		n -- dimension of vectors _x and _y
	*/
	dp_t _ddotv(const dptr _x, const dptr _y, int n);

	/*
	 @Purp: assign a scalar to array [a(m,n) = s]
	 @Type: single precision floating point[SFP]
	 @Para:
		a -- read-write pointer to array
		s -- scalar will be assigned to a(m,n)
		m,n -- dimension of a(m, n)
	 */
	void _sufa(sptr* a, sp_t s, int m, int n);
	/*
	 @Purp: assign a scalar to array [a(m,n) = s]
	 @Type: double precision floating point[DFP]
	 @Para:
		a -- read-write pointer to array
		s -- scalar will be assigned to a(m,n)
		m,n -- dimension of array a(m, n)
	 */
	void _dufa(dptr a, dp_t s, int m, int n);

	/*
	 @Purp: scalar-array multiplication[a:= s¡Áa(m, n)]
	 @Type: single precision floating point[32f]
	 @Para:
		s -- the scalar factor
		a -- pointer to array, overwritten by the result
		m,n -- dimension of array a(m, n)
	*/
	void _ssxa(Sfp s, Sfp* a, int m, int n);
	/*
	 @Purp: scalar-array multiplication[a:= s¡Áa(m, n)]
	 @Type: double precision floating point[DFP]
	 @Para:
	 	s -- the scalar factor
	 	a -- pointer to array, overwritten by the result
		m,n -- dimension of array a(m, n)
	*/
	void _dsxa(Dfp s, Dfp* a, int m, int n);

	/*
	 @Purp: small matrix-vector multiplication: x := A(m,n)x(n)
	 @Type: single precision floating point[SFP]
	 @Para:
	     A -- pointer to the matrix with row-major format
	     x -- pointer to  the vector, overwritten by the result
	     m,n -- dimension of matrix a(m, n) and vector x(n)
	*/
	void _smv(const Sfp* a, Sfp* x, int m, int n);
	/*
	 @Purp: small matrix-vector multiplication: x := A(m,n)x(n)
	 @Type: double precision floating point[DFP]
	 @Para:
	     A -- pointer to the matrix with row-major format
	     x -- pointer to  the vector, overwritten by the result
	     m,n -- dimension of matrix a(m, n) and vector x(n)
	*/
	void _dmv(const Dfp* a, Dfp* x, int m, int n);

	/*
	 @Purp: calculate y(m)A(m,n)x(n) [Hasn't been implemented]
	 @Type: double precision floating point[DFP]
	 @Para:
	     y -- pointer to the left-handed vector
	     A -- pointer to the matrix with row-major format
	     x -- pointer to the right-handed vector
	     m,n -- dimension of matrix a(m, n), vector y(m), x(n)
	*/
	sp_t _svmv(const sptr y, const sptr A, const sptr x, int m, int n);
	/*
	 @Purp: calculate y(m)A(m,n)x(n)
	 @Type: double precision floating point[DFP]
	 @Para:
	     y -- pointer to the left-handed vector
	     A -- pointer to the matrix with row-major format
	     x -- pointer to the right-handed vector
	     m,n -- dimension of matrix a(m, n), vector y(m), x(n)
	*/
	dp_t _dvmv(const dptr y, const dptr A, const dptr x, int m, int n);

	/*
	 @Purp: small matrix multiplication: C(m,n) := A(m,k)B(k,n)
	 @Type: single precision floating point[SP]
	 @Para:
		c, a, b -- pointer to matrices C, A and B, row-major format
		m, k, n -- dimension of matrices C, A and B
	*/
	void _smm(sptr c, const sptr a, const sptr b, int m, int k, int n);
	/*
	 @Purp: small matrix multiplication: C(m,n) := A(m,k)B(k,n)
	 @Type: double precision floating point[DP]
	 @Para:
	    c, a, b -- pointer to matrices C, A and B, row-major format
	    m, k, n -- dimension of matrices C, A and B
	*/
	void _dmm(dptr c, const dptr a, const dptr b, int m, int k, int n);

	/*
	 @Purp: rectangle matrix multiplication: S(n,n) := (A^T)A
	 @Type: double precision floating point[DP]
	 @Para:
		a, s -- pointer to matrices A and S, row-major format
		m, n -- dimension of matrix A
	*/
	void _dmtm(dptr s, const dptr a, int m, int n);
	/*
	 @Purp: rectangle matrix multiplication: S(m,m) := AA^T
	 @Type: double precision floating point[DP]
	 @Para:
		a, s -- pointer to matrices A and S, row-major format
		m, n -- dimension of matrix A
	*/
	void _dmmt(dptr s, const dptr a, int m, int n);
	/*
	 @Purp: rectangle matrix multiplication: S(n,n) := (A^T)A
	 @Type: single precision floating point[SP]
	 @Para:
		a, s -- pointer to matrices A and S, row-major format
		m, n -- dimension of matrix A
	*/
	void _smtm(sptr s, const sptr a, int m, int n);
	/*
	 @Purp: rectangle matrix multiplication: S(m,m) := AA^T
	 @Type: single precision floating point[SP]
	 @Para:
		a, s -- pointer to matrices A and S, row-major format
		m, n -- dimension of matrix A
	*/
	void _smmt(sptr s, const sptr a, int m, int n);

	/*
	 @Purp: triple square matrix multiplication: R(n,n) := A(n,n)B(n,n)C(n,n)
	 @Type: double precision floating point[DP]
	 @Para:
	    r, c, a, b -- pointer to matrices C, A and B, row-major format
	    n -- dimension of matrices
	*/
	void _dtplsm(dptr r, const dptr a, const dptr b, const dptr c, int n);

	/*
	 @Purp: sum of two scalar-array multiplication [x := p*x(m, n) + q*y(m, n)]
	 @Type: double precision floating point[DFP]
	 @Para:
		p, q -- the two scalar factors
		x    -- pointer to first array, which will be overwritten by the results
		_y   -- read-only pointer to the sencond array
		m,n  -- dimension of two arrays x and y
	*/
	void _dpxqy(Dfp p, Dfp* x, Dfp q, Dfp* _y, int m, int n);

	/*
	 @Purp: transpose a square matrix A [A := trans(A)] in place
	 @Type: double precision floating point[DFP]
	 @Para:
		a -- pointer to matrix, will be overwritten by it's transposition
		n -- dimension of matrix a(n,n)
	 */
	void _dtrpsm(Dfp* a, int n);
	void _strpsm(Sfp* a, int n);

	/*
	 @Purp: transpose a general matrix a [b = trans(a)]
	 @Type: double precision floating point[DFP]
	 @Para:
	    _a -- pointer to the source matrix
	     b -- pointer to the transpose
	    m,n -- dimension of matrix a(m,n), b(n,m)
	*/
	void _dtrpm(const Dfp* _a, Dfp* b, int m, int n);
	void _strpm(const Sfp* _a, Sfp* b, int m, int n);

	/*
	@Purp: array contraction or Frobenius inner product. Do entrywise product after A minus m1 and B minus m2, return summation of all entries of the resulted array: [res = SUM{(A-m1)*(B-m2)}]
	@Type: double precision floating point[DFP]
	@Para:
	_a -- pointer to the source matrix A
	_b -- pointer to the source matrix B
	m1,m2 -- two scalar will be subtracted from A and B
	m,n -- dimension of matrix A(m,n), B(m,n)
	 k -- flag with default value 2, return only SUM(A-m1) if k is assigned to 1. 
	*/
	Dfp _dctr(Dfp* _a, DFP m1, Dfp* _b, DFP m2,
		int m, int n, int k = 2);
	/*
	@Purp: array contraction or Frobenius inner product.
	Do entrywise product after A minus m1 and B minus m2, 
	then return summation of all entries of the resulted array: [res = SUM{(A-m1)*(B-m2)}]
	@Type: single precision floating point[SFP]
	@Para:
	_a -- pointer to the source matrix A
	_b -- pointer to the source matrix B
	m1,m2 -- two scalar will be subtracted from A and B
	m,n -- dimension of matrix A(m,n), B(m,n)
	k -- flag with default value 2, return only SUM(A-m1) if k is assigned to 1.
	*/
	Dfp _sctr(Sfp* _a, DFP m1, Sfp* _b, DFP m2,
		int m, int n, int k = 2);
	/*
	@Purp: array contraction or Frobenius inner product.
	entrywise product after A minus v.
	return summation of all entries of the resulted array: [res = SUM{A[*A]}]
	@Type: vector with DFP
	@Para:
	_A -- pointer to the source matrix A = {v1, v2,..., vn}
	v -- pointer to vector will be subtracted form A
	m,n -- m is the dim of vector, n is the number of vectors
	k -- flag with default value 2, return only SUM(A-v) if k is assigned to 1.
	*/
	Dfp _dctrv(Dfp* _A, DFP* v, int m, int n, int k = 2);

	/*
	@Purp: 2-norm(Frobenius) of a array: [res = SQRT(SUM(A*A))]
	@Type: single precision floating point[SFP]
	@Para:
	_a -- pointer to the source matrix A
	m,n -- dimension of matrix A(m,n)
	*/
	Sfp _snorm_2(Sfp* _a, int m, int n);
	/*
	@Purp: 2-norm(Frobenius) of a array: [res = SQRT(SUM(A*A))]
	@Type: double precision floating point[DFP]
	@Para:
	_a -- pointer to the source matrix A
	m,n -- dimension of matrix A(m,n)
	*/
	Dfp _dnorm_2(Dfp* _a, int m, int n);

}
/* ******
 * Data Solve linear equation system
 * ******/
namespace les{
	/* @Purp: inverse of A(n,n)
	 * @Type: double precision floating point[DFP]
	 * @Param:
		1. pai -- pointer to matrix A(n,n), which will be overwritten by the transpose of its inverse A^{-1}.
	 * @Return: 1 for success, 0 for singular A(n,n), -i(i = 1, 2, ..., n) for the ith diagnal entry equals to 0.
	 */
	int _dginv(Dfp* pai, int n);
	/* @Purp: inverse of A(n,n)
	 * @Type: single precision floating point[SFP]
	 * @Param:
		1. pai -- pointer to matrix A(n,n), which will be overwritten by the transpose of its inverse A^{-1}.
	 * @Return: 1 for success, 0 for singular A(n,n), -i(i = 1, 2, ..., n) for the ith diagnal entry equals to 0.
	 */
	int _sginv(Sfp* pai, int n);

	/* @Purp: Gauss method to solve A(n,n)x = b with pivoting
	 * @Type: double precision floating point[DFP]
	 * @Param:
		1. pAi -- pointer to coeff. matrix A(n,n), which will be overwritten by the transpose of its inverse A^{-T}.
		2. pbx -- pointer to right hand column matrix b, which  will  be overwritten by the solution.
	 * @Return: 1 for success, 0 for singular A(n,n), -i(i = 1, 2, ..., n) for the ith diagnal entry equals to 0.  
	 */
	int _dgesv(dptr pAi, dptr pbx, int n);
	/* @Purp: Gauss method to solve A(n,n)x = b with pivoting
	 * @Type: single precision floating point[SFP]
	 * @Param:
		1. pAi -- pointer to coeff. matrix A(n,n), which will be overwritten by the transpose of its inverse A^{-T}.
		2. pbx -- pointer to right hand column matrix b, which  will  be overwritten by the solution.
	 * @Return: 1 for success, 0 for singular A(n,n), -i(i = 1, 2, ..., n) for the ith diagnal entry equals to 0.  
	 */
	int _sgesv(sptr pAi, sptr pbx, int n);

	/* @Purp: Gauss method to solve overdetermined (m >= n) A(m,n)x(n) = b(m) with pivoting
	 * @Type: double precision floating point[DFP]
	 * @Param:
		1. pAi -- pointer to coeff. matrix A(m,n).
		2. pbx -- pointer to right hand column matrix b, which first n element will  be overwritten by the solution.
	 * @Return: 1 for success, 0 for singular A(m,n), -i(i = 1, 2, ..., n) for the ith diagnal entry equals to 0.  
	 */
	int _dgeosv(dptr pAi, dptr pbx, int m, int n);
	/* @Purp: Gauss method to solve overdetermined (m >= n) A(m,n)x = b with pivoting
	 * @Type: single precision floating point[SFP]
	 * @Param:
		1. pAi -- pointer to coeff. matrix A(m,n).
		2. pbx -- pointer to right hand column matrix b, which first n elements will  be overwritten by the solution.
	 * @Return: 1 for success, 0 for singular A(m,n), -i(i = 1, 2, ..., n) for the ith diagnal entry equals to 0.  
	 */
	int _sgeosv(sptr pAi, sptr pbx, int m, int n);

	/* @Purp: General SVD for computing U, ¦² and Vt of input matrix A (row-major).
	   @Type: double/single precision floating point[DFP]
	   @Param:
		1. pA -- pointer to coeff. matrix A(m,n), which will be overwritten by U(m,n).
		2. pS -- pointer to singular value array S(n) with descending order.
		3. pVt -- pointer to Vt(n,n)
		4. m,n -- dimension of A.
	   @Return: 0 for success, -2 for misconvergence, -4 for dimension inconformity w.r.t. n */
	int _dgesvd(dptr pA, dptr p¦², dptr pVt, int m, int n); int _sgesvd(sptr pA, sptr p¦², sptr pVt, int m, int n);

	/* @Purp: General SVD for computing U, ¦² and Vt of input matrix A. A is stored in row-major if it's square matrix, else A is passed in column-major with preserving initial shape[will be decrepe
	 * @Type: double precision floating point[DFP]
	 * @Param:
		1. pA -- pointer to coeff. matrix A(m,n), which will be overwritten by U(m,n).
		2. pS -- pointer to singular value array S(n) with descending order.
		3. pVt -- pointer to Vt(n,n)
		4. m,n -- dimension of A. The pass order of m and n should be exchangded if A is not a square matrix.
	 * @Return: 1 for success, 0 for misconvergence, -4 for dimension inconformity w.r.t. n */
	int _dgsvd(Dfp* pA, Dfp* pS, Dfp* pVt, int m, int n);

	/* @Purp: Fast SVD A = U¦²Vt without computing U matrix, where A is stored in row-major.
	   @Type: double(_d)/single(_s) precision floating point[DFP]
	   @Pars:
		1. pA -- pointer to matrix A(m,n), which will be overwritten by a bidiagnal matrix;
		2. p¦² -- pointer to singular value array ¦²(n) with descending order;
		3. pVt -- pointer to Vt(n,n);
		4. m,n -- dimension of A.
	   @Return: 0 for success, -2 for misconvergence, -4 for dimension inconformity w.r.t. n*/
	int _dfgsvd(dptr pA, dptr p¦², dptr pVt, int m, int n); int _sfgsvd(sptr pA, sptr p¦², sptr pVt, int m, int n);

	/* @Purp: Fast SVD A = U¦²Vt without computing U matrix. A is stored in row-major if it's square matrix, else A is passed in column-major with preserving initial shape
	 * @Type: double precision floating point[DFP]
	 * @Param:
		1. pA -- pointer to coeff. matrix A(m,n), which will be overwritten by a bidiagnal matrix.
		2. pS -- pointer to singular value array S(n) with descending order.
		3. pVt -- pointer to Vt(n,n)
		4. m,n -- dimension of A. The pass order of m and n should be exchangded if A is not a square matrix.
	 * @Return: 1 for success, 0 for misconvergence, -4 for dimension inconformity w.r.t. n*/
	int _dfsvd(dptr pA, dptr pS, dptr pVt, int m, int n);
	/* @Purp: Fast SVD A = U¦²Vt. A is stored in row-major if it's square matrix, else A is passed in column-major with preserving initial shape
	 * @Type: single precision floating point[SFP]
	 * @Param:
		1. pA -- pointer to coeff. matrix A(m,n), which will be overwritten by a bidiagnal matrix.
		2. pS -- pointer to singular value array S(n) with descending order.
		3. pVt -- pointer to Vt(n,n)
		4. m,n -- dimension of A. The pass order of m and n should be exchangded if A is not a square matrix.
	 * @Return: 1 for success, 0 for misconvergence, -4 for dimension inconformity w.r.t. n*/
	int _sfsvd(sptr pA, sptr pS, sptr pVt, int m, int n);
	/* @Purp: Solve A(m, n)x(n) = 0 using fast SVD, \ref{_?gsvd2}
	 * @Type: single precision floating point[SFP]
	 * @Param:
		1. pA -- pointer to coeff. matrix A(m,n), which will be overwritten by a bidiagnal matrix.
		2. px -- pointer to k-th singular vector.
		3. m,n -- dimension of A. The pass order of m and n should be exchangded if A is not a square matrix.
		4. k -- index of singular vector will be retrieved by px from Vt(n, n), default is the right one assaciated with the smallest singular value.
	 * @Return: 1 for success, 0 for misconvergence, -4 for dimension inconformity w.r.t. n .*/
	int _sgsvdsv(sptr pA, sptr px, int m, int n, int k = -1);
	/* @Purp: Solve A(m, n)x(n) = 0 using fast SVD, \ref{_?gsvd2}
	 * @Type: double precision floating point[DFP]
	 * @Param:
		1. pA -- pointer to coeff. matrix A(m,n), which will be overwritten by a bidiagnal matrix.
		2. px -- pointer to k-th singular vector.
		3. m,n -- dimension of A. The pass order of m and n should be exchangded if A is not a square matrix.
		4. k -- index of singular vector will be retrieved by px from Vt(n, n), default is the right one assaciated with the smallest singular value.
	 * @Return: 1 for success, 0 for misconvergence, -4 for dimension inconformity w.r.t. n. */
	int _dgsvdsv(dptr pA, dptr px, int m, int n, int k = -1);

	/* @Purp: LU decomposition for square matrix A(n,n)
	 * @Type: double precision floating point[DFP]
	 * @Param:
		1. pA -- pointer to coeff. matrix A(n,n), which will be overwritten by a bidiagnal matrix.
		2. ppiv -- pointer to the row permutation.
		3. n -- dimension of A
	 * @Return: 
			0 for failure since there is a row of zeros in A, 
			+1/-1 for success with an even/odd number of row interchanges*/
	int _dLU(dptr pA, int* ppiv, int n);
	int _sLU(sptr pA, int* ppiv, int n);

	/* @Purp: cholesky decomposition A = L*L^T
	 * @Type: dp_t[inout]
	 * @Param:
		1. pA -- pointer to n-by-n positive sym mat A, which lower triangle part will be overwrriten by L
		2. pp -- pointer to diagnal entry vector;
		3. n   -- rows and cols of A
	   @Return: 1 for success, -i for the ith diagnal entry in A equals to zero.
	 */
	int _dcholy(dptr pA, int n);
	int _scholy(sptr pA, int n);
	int _dcholy2(dptr pA, dptr pp, int n);
	/* @Purp: solve linear equations Ax = b, where A is a positive-definite symmetrix matrix and A should be decomposed by calling _dcholy in advance.
	 * @Type: dp_t[inout]
	 * @Param:
		1. pA -- pointer to n-by-n positive sym mat A;
		2. pbx -- pointer to r.h.s. vector and solution vector;
		3. n   -- rows and cols of A.
	 */
	void _dcholsv(dptr pA, dptr pbx, int n);

	/* @Purp: forward substitution, use lower part
	 * @Type: dp_t[inout]
	 * @Param:
		1. pLU -- pointer to decomposed lower and upper triangle square matrix
		2. pbx -- pointer to rhv[in] and results[out]
		3. n   -- rows and cols for LU matrix
	   @Return: 1 for success, -i for the ith diagnal entry in LU matrix equals to zero.
	 */
	int _dfwsubt(dptr pLU, dptr pbx, int n);
	/* @Purp: backward substitution, use upper part
	 * @Type: dp_t[inout]
	 * @Param:
		1. pLU -- pointer to decomposed lower and upper triangle square matrix
		2. pbx -- pointer to rhv[in] and results[out]
		3. n   -- rows and cols for LU matrix
	 * @Return: 1 for success, -i for the ith diagnal entry in LU matrix equals to zero.
	*/
	int _dbwsubt(dptr pLU, dptr pbx, int n);
}
/* ******
 * Data modeling subroutines
 * ******/
namespace dms{
	/* @Purp: linear fitting. Fit input data as a linear model y = a1 * x + a2.
	   @Type: double precision floating point[DFP]
	   @Para: 
		1. _pData -- source data, saved as(x1, y1, x2, y2, ..., xn, yn);
		2. a -- the best fitted paramters, a1 and a2;
		3. sigma -- the variances in the estimates;
		4. n -- the number of source data points
	   @Retu: the merit value, a good-fit means this function will returns a smaller merit value, and the fitted line is vertical(x + a(2) = 0) if the returned value is a negative.
	 */
	Dfp _dlftab(Dfp* _pData, Dfp* a, Dfp* sigma,
		int n);

	void _dqsft3(dptr _pX, dptr _pY, dptr pcoef, int m, int n);

	void _dcsft2(dptr _pX, dptr _pY, dptr pcoed, int n);
}

/* ******
 * compution geometry subroutines
 * ******/
namespace cgs{
	/* @Purp: find the pedal of line through a point and the input line
	 * @Type: double precision floating point[DFP]
	 * @Para:
		1. _l -- pointer to the input line
		2. p  -- pointer to the fixed point, overwritten by the found pedal when function return
		3. m  -- the dimension of point
	 * @Retu: 
	   3 for success with the pedal lying on the line l but locates outside of the input first endpoint, 
	   2 for success with the pedal lying on the line l but locates outside of the input second endpoint, 
	   1 for success with the pedal lying on the line l,
	   0 for the pedal is coincident with input point p,
	  -1 for failure where the endpoints of input line are coincident.
	 */
	int _dpedal(Dfp* _l, Dfp* p, int m);

	/* @Purp: Rodrigues transformation used to convert a rotation matrix to a rotation vector or vice versa.
	* @Type: single precision floating point[SFP]
	* @Para:
	1. _rv/_rm -- pointer to input rot. vector(3x1) or rot. matrix(3x3)
	2. rm/rv  -- pointer to corresponding output rot. matrix(3x3) or rot. vector(3x1) with returned rot. angle
	*/
	void _srodrgvm(const Sfp* _rv, Sfp* rm);
	Sfp  _srodrgmv(const Sfp* _rm, Sfp* rv);
	/* @Purp: Rodrigues transformation used to convert a rotation matrix to a rotation vector or vice versa.
	* @Type: double precision floating point[DFP]
	* @Para:
	1. _rv/_rm -- pointer to input rot. vector(3x1) or rot. matrix(3x3)
	2. rm/rv  -- pointer to corresponding output rot. matrix(3x3) or rot. vector(3x1) with returned rot. angle
	*/
	void _drodrgvm(const Dfp* _rv, Dfp* rm);
	Dfp  _drodrgmv(const Dfp* _rm, Dfp* rv);
}

/* ***************************************************
 *Desc: special subroutines and functions
 *Name: <domain>_<data-type><name>[s|v|m]
	 @domain -- cvg(convergence number);
	 @data-type -- s(single), d(double);
	 @name -- abbreviation of used algorithm or purpose.
	 @s|v|m -- options to indicate the input are scalar, vector and matrix respectively.
 * **************************************************/
namespace spe{
	/*
	@Purp: calc. cvg number by using normalized square sum
	@Type: double precision floating point[DFP]
	@Para:
	    _x0 -- pointer to vector of reference data
	    _dx -- pointer to vector of object data
	      n -- dimension of vectors _x0 and _dx
	*/
	Dfp cvg_dnssv(Dfp* _x0, Dfp* _dx, int n);

	/*
	 @Purp: get a sub-array from input array A
	 @Type: double precision floating point[DFP]
	 @Para:
		_pA  -- pointer to input array A
		pSub -- pointer to output sub-array
		m, n -- rows and cols of A
		x, y -- current position
		rx,ry-- radiuses of sub-array in row and col directions
		dim  -- number of DFPs in each entry
	*/
	void dsubm(dptr _pA, dptr pSub, int m, int n, int x, int y, int rx, int ry, int dim);
	void ssubm(sptr _pA, sptr pSub, int m, int n, int x, int y, int rx, int ry, int dim);

	/*
	@Purp: get a sub-array from input array A
	@Type: double precision floating point[DFP]
	@Para:
	_pA  -- pointer to input array A
	pSub -- pointer to output sub-array
	m, n -- rows and cols of A
	r0,r1-- indices of start and end rows
	c0,c1-- indices of start and end cols
	*/
	void dsubm_r(dptr _pA, dptr pSub, int m, int n, int r0, int r1, int c0, int c1);
	void ssubm_r(sptr _pA, sptr pSub, int m, int n, int r0, int r1, int c0, int c1);

	/*
	@Purp: copy a block from A to B
	@Type: double precision floating point[DFP]
	@Para:
	_pA  -- pointer to input array A
	m1, n1 -- rows and cols of A
	pB -- pointer to output array B
	m2, n2 -- rows and cols of B
	r0,r1-- indices of start and end rows
	c0,c1-- indices of start and end cols
	*/
	void _dcopy(dptr _pA, int m1, int n1, dptr pB, int m2, int n2, int r0, int r1, int c0, int c1);
	void _scopy(sptr _pA, int m1, int n1, sptr pB, int m2, int n2, int r0, int r1, int c0, int c1);

	void _dcpy(dptr _pA, int m1, int n1, dptr pB, int m2, int n2, int* r1 = nullptr, int* r2 = nullptr);

	/*
	@Purp: padding array A to form array B
	@Type: double precision floating point[DFP]
	@Para:
	_pA -- pointer to input array A
	 pB -- pointer to padded array B
	m, n -- rows and cols of A
	bm, bn -- padding size in row and column directions
	btype -- padding type: 0 for zero value padding, 1 for edge value padding 
	*/
	void dpad(dptr _pA, dptr pB, int m, int n, int bm, int bn, int btype = 0);
	/*
	@Purp: padding array A to form array B
	@Type: uchar to dp_t
	@Para:
	_pA -- pointer to input array A
	pB -- pointer to padded array B
	m, n -- rows and cols of A
	bm, bn -- padding size in row and column directions
	btype -- padding type: 0 for zero value padding, 1 for edge value padding
	*/
	void dpad_u8(ucptr _pA, dptr pB, int m, int n, int bm, int bn, int btype = 0);

	/*
	@Purp: deconvolution by using element-wise division
	@Type: double precision floating point[DFP]
	@Para:
	pA  -- pointer to input and output complex array A
	_pK -- pointer to input FFT kernal array K
	n   -- the number of complex element in A and K
	*/
	void deconv(dptr pA, dptr _pK, int n);

	/*
	@Purp: complex number operation acoording to the input operator "optr": pIO := pIO '/' '*' '-' '+' pIn
	@Type: double precision floating point[DFP]
	@Para:
	pIO  -- pointer to input and output complex array. its  format can be array of _Ty[2] or array of {real, imag}
	optr -- operator, can be one of '+', '-', '*', or '/'
	pIn  -- pointer to right side input complex array
	n   -- the number of complex element in pIO and pIn
	*/
	void dcmplx(dptr pIO, char optr, dptr pIn, int n);
	void scmplx(sptr pIO, char optr, sptr pIn, int n);

	/*
	@Purp: calc. 2d strain
	@Type: double precision floating point[DFP]
	@Para:
		_pF -- pointer to deformation gradient F[2x2]
		pE -- pointer to strain {E11,E22,E12}
		type -- strain type: 0 for Lagrange, 1 for Eulerian, 2 for Infinitesimal, 3 for Engineering
	*/
	void dstrain2(dptr _pF, dptr pE, int type);
	/*
	@Purp: calc. 3d strain
	@Type: double precision floating point[DFP]
	@Para:
	_pF -- pointer to deformation gradient F[3x3]
	pE -- pointer to strain {E11,E22,E33,E12,E23,E31}
	type -- strain type: 0 for Lagrange, 1 for Eulerian, 2 for Infinitesimal, 3 for Engineering
	*/
	void dstrain3(dptr _pF, dptr pE, int type);
}

namespace img{

	/*
	@Purp: the mean of image, return SUM(I)/(m*n)
	@Type: double precision floating point[DFP]
	@Para:
	_pImg -- pointer to image matrix
	m,n -- rows and cols of _pImg
	*/
	Dfp _mean_f32c1(Sfp* _pImg, int m, int n);
	Dfp _mean_f64c1(Dfp* _pImg, int m, int n);
	
	/*
	 * Purp: Gaussian convolution
	 * Type: unsigned char[in], double[out]
	 * Para:
		1. pI -- pointer to padded src(m+2r,n+2r) with uchar type, src must be padded in advance with radius r for each side;
		2. pO -- pointer to res(m,n) with double type;
		3. m,n-- rows and cols of original image;
		4. r -- gauss kernel radius, which size will be (2r+1)x(2r+1);
		5. ¦Ò -- must be a static array, double ¦Ò[2] = {¦Òx, ¦Òy}, default is {2.0, 2.0};
		6. ¦Ì -- double ¦Ì[2] = {¦Ìx, ¦Ìy}, default is {0.0, 0.0}.
	 */
	void _dconvgu(ucptr pI, dptr pO, int m, int n, int r, dptr ¦Ò = nullptr, dptr ¦Ì = nullptr);

	/*
	* Purp: cross correlation
	* Type: unsigned char[in], double[out]
	* Para:
	1. pF -- pointer to reference image F;
	2. pG -- pointer to target image G;
	3. n-- dimension of F and G;
	*/
	dp_t _dcorru(ucptr pF, ucptr pG, int n);

	/*
	* Purp: image integration subroutine
	* Type: unsigned char[in], integer[out]
	* Para:
	1. pI -- pointer to input image;
	2. pO -- pointer to output integral image;
	3. m,n-- dimension of I and O;
	*/
	void _iitguc(ucptr pI, iptr pO, int m, int n);

	/*
	* Purp: copy image. if inc is not nullptr and assigned by values greater then 1, the input image will be reduced acoording to inc parameter
	* Type: unsigned char[in], unsigned char[out]
	* Para:
	1. pI -- pointer to input image;
	2. m1,n1-- dimensions of I;
	3. pO -- pointer to output image;
	4. m2,n2-- dimensions of O;
	5. inc[opt] -- increment step in x and y direction
	*/
	void _ucpy(ucptr pI, int m1, int n1, ucptr pO, int m2, int n2, int* inc = nullptr);

	/*
	* Purp: assign val to IO
	* Type: unsigned char[inout]
	* Para:
	1. pIO -- pointer to input image;
	2. m,n -- dimension of IO;
	3. val[opt] -- if specified means pIO := val
	4. pr[opt] -- pointer to close range {x_begin, x_end, y_begin, y_end}, if specified means pIO := val in the specified range
	*/
	void _uassign(ucptr pIO, int m, int n, unsigned char val = 0, int* pr = nullptr);

	/*
	* Purp: convert I[uchar] into O(double)
	* Type: unsigned char[in], dp(out)
	* Para:
	1. pO -- pointer to output image;
	1. pI -- pointer to input image;
	2. m,n -- dimension of I/O;
	3. normalized[opt] -- default is 1 for dividing O by 255 after conversion, else convert I to O directly 
	*/
	void _dcvtuc(dptr pO, ucptr pI, int m, int n, int normalized = 1);

	/*
	 * Purp: calc. sum of squared difference for an input array with given mean value
	 * Type: uchar [in], dp_t[out]
	 * Param:
		1. pIn -- pointer to input array
		2. m,n -- rows and cols of input array
		3. pAvg -- pointer to mean value[optional], if it is given than the subroutine uses it directly, else it will be calculated.
	 * Return: ssd value
	 */
	dp_t _dssdmu(ucptr pIn, int m, int n, dptr pAvg);
}

namespace utility
{
	/*
	* Purp: calc. sum of squared difference for an input array with given mean value
	* Type: dp_t[inout]
	* Param:
	1. m,n -- rows and cols of input array
	2. pIn -- pointer to input array
	3. pAvg -- mean value of pIn
	4. pMask -- mask for pIn [optional], if pMask is specified, then this function just calculates those elements with true mask value
	* Return: ssd value
	*/
	dp_t _dssdmd(int m, int n, dptr pIn, dp_t Avg, ucptr pMask = nullptr);
}

#ifdef __cplusplus
}
#endif
} 

namespace fblas = fkl::bla;
namespace flapk = fkl::les;
namespace fimgk = fkl::img;
namespace fcges = fkl::cgs;
namespace fspec = fkl::spe;
namespace futil = fkl::utility;