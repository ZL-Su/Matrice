/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2020, Zhilong(Dgelom) Su, all rights reserved.

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
***********************************************************************/
#pragma once
#include <functional>
#include <exception>
#include <string>
#include <vector>
#include "_macros.h"

#if defined(_MSC_VER) && _MSC_VER <= 1900
#define __func__ __FUNCTION__
#endif

DGE_MATRICE_BEGIN _DETAIL_BEGIN
/**
 * \record source code location
 * \Example: _Source_location _Loc{__func__, __FILE__, __LINE__}
 */
struct _Source_location {
	const char* _func = nullptr;
	const char* _file = nullptr;
	long        _line;
};
#define __exceploc__ {__func__, __FILE__, __LINE__}

/**
 * \exception process
 */
struct _Exception_wrapper
{
	using msg_type = std::string;
	using msg_list = std::vector<msg_type>;
	using loc_type = _Source_location;

	class error : public std::exception {
		using _Myt = error;
		using _Mybase = std::exception;
	public:
		using _Mybase::exception;
		error(const loc_type& _Loc, const msg_type& _Msg="None") noexcept
			:_Myloc(_Loc), _Mybase(_Msg.c_str()) {}

		/**
		 * \return exception location
		 */
		MATRICE_HOST_INL loc_type location() const noexcept {
			return (_Myloc);
		}

		/**
		 * \return exception message 
		 */
		MATRICE_HOST_INL msg_type message() const noexcept {
			return std::string(" [MSG]: ") + this->what() 
				+ "\n [Func]: " + _Myloc._func 
				+ "\n [Line]: " + std::to_string(_Myloc._line) 
				+ "\n [File]: " + _Myloc._file;
		}

	private:
		const void* _Mycaller = nullptr;
		msg_type _Mymsg = this->what();
		msg_list _Mymsgstack;
		loc_type _Myloc;
	};

	using handler = std::function<void(const msg_type&, const char*)>;
	static MATRICE_HOST_INL auto warning(loc_type _Loc, msg_type _Msg);
};
_DETAIL_END
using exception = detail::_Exception_wrapper;

#define DGELOM_ERROR(...) \
	throw ::dgelom::exception::error(__exceploc__, \
			::dgelom::exception::msg_type(__VA_ARGS__))

#define DGELOM_CHECK(_Cond, ...) \
	if(!(_Cond)) { \
		DGELOM_ERROR(::dgelom::exception::msg_type(__VA_ARGS__)); \
	}

#define MATRICE_FAIL_TO_SPECIALIZATION \
static_assert(std::false_type::value, "Fail to specialization!");

DGE_MATRICE_END