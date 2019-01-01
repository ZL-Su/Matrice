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
#include<experimental/filesystem>
#include <unordered_map>
#include <vector>
#include "../util/utils.h"
#include "../util/genalgs.h"
#include "io.h"

DGE_MATRICE_BEGIN namespace io {
/**
 * \for filesystem namespace in STD library
 */
namespace fs = std::experimental::filesystem;

_DETAIL_BEGIN

struct folder_tag {};
struct loader_tag {};
struct file_tag {};

template<typename _Tag> struct _Collector {};
template<> struct _Collector<folder_tag> {
	MATRICE_HOST_FINL static auto get(std::string&& path) {
		std::vector<std::string> _Ret;
		if (fs::path _Path(path); fs::exists(_Path)) {
			fs::directory_iterator _End;
			for (decltype(_End) _Begin(_Path); _Begin != _End; ++_Begin) {
				if (fs::is_directory(_Begin->status())) {
					_Ret.emplace_back(_Begin->path().string().substr(path.size()));
				}
			}
			_Ret.shrink_to_fit();
		}
		return std::forward<decltype(_Ret)>(_Ret);
	}
};
template<> struct _Collector<file_tag> {
	MATRICE_HOST_FINL static auto get(std::string&& path) {
		std::vector<std::string> _Ret;
		if (fs::path _Path(path); fs::exists(_Path)) {
			fs::directory_iterator _End;
			for (decltype(_End) _Begin(_Path); _Begin != _End; ++_Begin) {
				if (fs::is_regular_file(_Begin->status())) {
					_Ret.emplace_back(_Begin->path().string().substr(path.size()));
				}
			}
			_Ret.shrink_to_fit();
		}
		return std::forward<decltype(_Ret)>(_Ret);
	}
};

template<typename _Tag = folder_tag> class _Dir_impl {};

template<> class _Dir_impl<folder_tag> {
public:
	using category = folder_tag;
	using value_type = std::string;
	using path_type = fs::path;
	using container = std::vector<value_type>;
	/**
	 * \_Root must has form "./folder_name" or "/folder_name" 
	 */
	_Dir_impl(value_type _Root, path_type&& _Workdir = fs::current_path())
		:_Mypath(_Workdir) {
		if (_Root[0] == '.') _Mypath.concat(_Root.begin() + 1, _Root.end());
		else _Mypath.concat(_Root.begin(), _Root.end());

		_Mysubfolders = _Collector<folder_tag>::get(_Mypath.string());
	}

	/**
	 * \return full path of _Idx-th folder
	 */
	MATRICE_HOST_INL auto operator[] (std::size_t _Idx) const {
		return (_Mypath.string() + folder(_Idx));
	}

	/**
	 * \return _Idx-th folder name under work-path _Mypath
	 */
	MATRICE_HOST_INL value_type folder(std::size_t _Idx) const {
#ifdef _DEBUG
		if (_Idx >= _Mysubfolders.size()) throw
			std::exception("_Idx over range of _Mysubfolders");
#endif
		return (_Mysubfolders[_Idx]);
	}

	/**
	 * \return size (number of subfolders) of current directory
	 */
	MATRICE_HOST_INL const auto size() const {
		return (_Mysubfolders.size());
	}

private:
	path_type _Mypath;
	container _Mysubfolders;
};

template<typename _Tag = loader_tag> class _Data_loader_impl{};

template<> class _Data_loader_impl<loader_tag> {
	using _Mydir_type = _Dir_impl<folder_tag>;
	using _Myt = _Data_loader_impl;

	struct _Loader_iterator {
		_Loader_iterator(const std::add_pointer_t<_Myt> _This)
			:_Mythis(_This), _Mybatchs(_Mythis->batch_size()){
		}

		MATRICE_HOST_FINL auto operator*() const {
			_Mydir_type::container _Ret;
			for (auto _Idx = 0; _Idx < _Mybatchs; ++_Idx) {
				_Ret.emplace_back(_Mythis->directory()[_Idx] + _Mythis->file_names(_Idx)[_Mypos]);
			}
			return std::forward<decltype(_Ret)>(_Ret);
		}

	private:
		std::size_t _Mypos = 0;
		std::size_t _Mybatchs = 0;
		std::add_pointer_t<_Myt> _Mythis = nullptr;
		_Mydir_type::container::iterator _Myitr;
	};
public:
	using category = loader_tag;
	using dir_type = _Mydir_type;
	using iterator = _Loader_iterator;

	_Data_loader_impl(const _Mydir_type& _Dir)
		: _Mydir(_Dir) {
		_Collect_fnames();
	}
	_Data_loader_impl(_Mydir_type&& _Dir)
		: _Mydir(std::forward<_Mydir_type>(_Dir)) {
		_Collect_fnames();
	}

	MATRICE_HOST_INL auto begin() const {

	}
	MATRICE_HOST_INL auto end() const {

	}
	/**
	 * \forward iterate to retrieve data paths
	 */
	MATRICE_HOST_INL auto forward() const {
		_Mypos++;
		return std::make_tuple();
	}
	/**
	 * \reverse iterate to retrieve data paths
	 */
	MATRICE_HOST_INL auto reverse() const {
		_Mypos--;
		return std::make_tuple();
	}

	MATRICE_HOST_FINL dir_type& directory() {
		return (_Mydir);
	}
	MATRICE_HOST_FINL const dir_type& directory() const {
		return (_Mydir);
	}
	MATRICE_HOST_FINL std::size_t batch_size() const {
		return (_Mydir.size());
	}

	/**
	 * \return all file names in currenct work path
	 */
	MATRICE_HOST_FINL std::vector<_Mydir_type::container>& file_names() {
		return (_Mynames);
	}
	MATRICE_HOST_FINL const std::vector<_Mydir_type::container>& file_names() const {
		return (_Mynames);
	}

	/**
	 * \return all file names in _Idx-th subfolder for currenct work path
	 */
	MATRICE_HOST_FINL _Mydir_type::container& file_names(std::size_t _Idx) {
#ifdef _DEBUG
		if (_Idx >= _Mynames.size()) throw
			std::exception("_Idx over range of _Data_loader_impl<>::_Mynames.");
#endif // _DEBUG

		return (_Mynames)[_Idx];
	}
	MATRICE_HOST_FINL const _Mydir_type::container& file_names(std::size_t _Idx) const {
#ifdef _DEBUG
		if (_Idx >= _Mynames.size()) throw
			std::exception("_Idx over range of _Data_loader_impl<>::_Mynames.");
#endif // _DEBUG
		return (_Mynames)[_Idx];
	}

private:
	MATRICE_HOST_INL void _Collect_fnames() {
		_Mynames.resize(_Mydir.size());
		for (const auto _Idx : range(0, _Mynames.size()))
			_Mynames[_Idx] = _Collector<file_tag>::get(_Mydir[_Idx]);
	}
	_Mydir_type _Mydir;
	mutable std::size_t _Mypos = 0;
	std::vector<_Mydir_type::container> _Mynames;
};

_DETAIL_END

using directory = detail::_Dir_impl<detail::folder_tag>;
using data_loader = detail::_Data_loader_impl<detail::loader_tag>;
using folder_collector = detail::_Collector<detail::folder_tag>;
using file_collector = detail::_Collector<detail::file_tag>;

} DGE_MATRICE_END