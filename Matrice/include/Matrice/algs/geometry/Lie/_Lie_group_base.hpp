/**************************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
#include "../core/matrix.h"
#include "_Lie_traits.hpp"
#include "_Lie_fwd.hpp"

DGE_MATRICE_BEGIN
_DETAIL_BEGIN
/**
 *\brief Base class for Lie groups, defines the common API.
 */
template<typename _Derived>
class _Lie_group_base {
    using _Myt = _Lie_group_base;
    using _Myderived = _Derived;
    using _Mytraits = traits<_Myderived>;
    using _Myprops = typename _Mytraits::properties;
public:
    static constexpr auto dim = _Myprops::dim;
    static constexpr auto dof = _Myprops::dof;

    using Lie_group_type = typename _Mytraits::group_type;
    using vector_type = typename _Mytraits::vector_type;
    using jacobian_type = typename _Mytraits::jacobian_type;

    /**
     *\brief Access the underlying data by pointer
     */
    MATRICE_GLOBAL_INL const auto data()const noexcept { 
        return _Mydata.data(); 
    }
    MATRICE_GLOBAL_INL auto data() noexcept {
        return _Mydata.data();
    }

    /**
     *\brief Get the adjoint of the group element.
     */
    MATRICE_GLOBAL_INL jacobian_type adj() const {

    }

private:
    MATRICE_GLOBAL_INL const _Myderived& derived() const { 
        return*static_cast<const _Myderived*>(this); 
    }
    MATRICE_GLOBAL_INL _Myderived& derived() {
        return*static_cast<_Myderived*>(this);
    }

protected:
    vector_type _Mydata; //stores underlying data
};

_DETAIL_END
DGE_MATRICE_END