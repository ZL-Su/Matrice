#pragma once
#include <core/matrix.h>
#include <algorithm>

DGE_MATRICE_BEGIN

template<typename _It>
void _Insertion_sort(const _It _Begin, const _It _End) noexcept {
	const auto _Len = (_End - _Begin);
	auto _First = _Begin;
	for (auto j = 1; j < _Len; ++j) {
		const auto _Key = *(_First + j);
		auto i = j - 1;
		auto _Val = *(_First + i);
		for (; i > 0 && _Val > _Key;) {
			*(_First + i + 1) = _Val;
			i -= 1;
			_Val = *(_First + i);
		}
		*(_First + i + 1) = _Key;
	}
}

DGE_MATRICE_END