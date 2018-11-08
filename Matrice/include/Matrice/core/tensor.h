#pragma once

#include "../private/tensor/_tensor.hpp"

DGE_MATRICE_BEGIN
template<typename _Ty, std::size_t... _Args>
using tensor = detail::_Tensor_impl<_Ty, _Args...>;
DGE_MATRICE_END
