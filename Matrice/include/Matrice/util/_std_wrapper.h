/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include <cstdint>
#include <queue>

namespace dgelom {
	using std::uint8_t;
	using std::uint16_t;
	using std::uint32_t;
	using std::uint64_t;
	using std::size_t;
	using diff_t = std::ptrdiff_t;

	template<typename... _Ts>
	using tuple = std::tuple<_Ts...>;
	template<typename... _Ts>
	using priority_queue = std::priority_queue<_Ts...>;

}
