/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

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

#include "../util/_macros.h"

DGE_MATRICE_BEGIN
/**********************************************************************
	This module refers to tag definitions to construct footstone of 
	template specialization and function overloading.
 **********************************************************************/
namespace tag {
	///<brief> tag definitions for categorizing matrix instance </brief>
	struct _Matrix_tag {};
	struct _Managed_matrix_tag : _Matrix_tag {};
	struct _Dynamic_matrix_tag : _Matrix_tag {};
	struct _Device_matrix_tag : _Matrix_tag {};

	///<brief> tag definitions for categorizing matrix iterator </brief>
	struct _Matrix_iterator_tag {};
	struct _Matrix_rwise_iterator_tag {};
	struct _Matrix_cwise_iterator_tag {};

	///<brief> tag definitions for categorizing matrix view </brief>
	struct _Matrix_view_tag {};
	struct _Matrix_row_view_tag : _Matrix_view_tag {};
	struct _Matrix_col_view_tag : _Matrix_view_tag {};
	struct _Matrix_block_view_tag : _Matrix_view_tag {};

	///<brief> tag definitions for categorizing expression </brief>
	struct _Expression_tag {};
	struct _Ewise_add_tag : _Expression_tag {}; // $x_i + y_i$
	struct _Ewise_sub_tag : _Expression_tag {}; // $x_i - y_i$
	struct _Ewise_mul_tag : _Expression_tag {}; // $x_i * y_i$
	struct _Ewise_div_tag : _Expression_tag {}; // $x_i / y_i$
	struct _Ewise_exp_tag : _Expression_tag {}; // $ e^{x_i} $
	struct _Ewise_log_tag : _Expression_tag {}; // $ log(x_i)$
	struct _Ewise_sin_tag : _Expression_tag {}; // $ sin(x_i)$
	struct _Ewise_cos_tag : _Expression_tag {}; // $ cos(x_i)$
	struct _Ewise_sqr_tag : _Expression_tag {}; // $ {x_i}^2 $
	struct _Ewise_sqrt_tag : _Expression_tag {}; // $ \sqrt{x_i} $
	struct _Ewise_max_tag : _Expression_tag {};
	struct _Ewise_min_tag : _Expression_tag {};
	struct _Ewise_abs_tag : _Expression_tag {};
	struct _Ewise_log2_tag : _Expression_tag {};
	struct _Ewise_log10_tag : _Expression_tag {};
	struct _Matrix_mul_tag : _Expression_tag {};// matrix mul.
	struct _Matrix_inv_tag : _Expression_tag {};// matrix inverse
	struct _Matrix_trp_tag : _Expression_tag {};// matrix transpose

	///<brief> variable tag definitions for function </brief>
	struct _Var_tag {};
	struct _Expression_eval_tag : _Var_tag {};
} 
DGE_MATRICE_END
