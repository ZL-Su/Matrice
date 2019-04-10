/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018, Zhilong(Dgelom) Su, all rights reserved.

This program is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#pragma once

#pragma warning (disable: 4067)

#include "_macros.h"

DGE_MATRICE_BEGIN
_CONDITIONS_BEGIN

// \condition expression: _Test_val satisfies _Bool_exp, for exampe _My_val < _Test_val 
#ifndef _COND_VAL(_Bool_exp)
#define _COND_VAL(_Bool_exp) [](const auto& _My_val){return (_My_val _Bool_exp);}
#endif

#ifndef _COND_LT(_Test_val)
#define _COND_LT(_Test_val) _COND_VAL(<_Test_val) //\less than _Test_val
#endif

#ifndef _COND_LQ(_Test_val)
#define _COND_LQ(_Test_val) _COND_VAL(<=_Test_val) //\less than or equal to _Test_val
#endif

#ifndef _COND_EQ(_Test_val)
#define _COND_EQ(_Test_val) _COND_VAL(==_Test_val) //\equal to _Test_val
#endif

#ifndef _COND_GT(_Test_val)
#define _COND_GT(_Test_val) _COND_VAL(>_Test_val) //\greater than _Test_val
#endif

#ifndef _COND_GQ(_Test_val)
#define _COND_GQ(_Test_val) _COND_VAL(>=_Test_val) //\greater than or equal to _Test_val
#endif

#ifndef _COND_EXCEPTION(_Test, _Msg)
#define _COND_EXCEPTION(_Test, _Msg) if(_Test)throw std::exception(_Msg);
#endif

_CONDITIONS_END
DGE_MATRICE_END

