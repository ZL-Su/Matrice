/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2021, Zhilong(Dgelom) Su, all rights reserved.

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
*********************************************************************/
#pragma once
#include <iostream>
#include <iosfwd>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <functional>
#include <type_traits>
#if _MSC_VER >= 1926
#include <filesystem>
#else
#include <experimental/filesystem>
#endif
#include "core/matrix.h"
#include "core/vector.h"
#include "util/utils.h"
#include "util/genalgs.h"
#include "util/_conditional_macros.h"
#include "util/_exception.h"
#include "private/_range.h"
#include "_dir.hpp"
#include "text_loader.hpp"

DGE_MATRICE_BEGIN
/**
 * \for filesystem namespace in STD library
 */
#if _MSC_VER >= 1926
	namespace fs = std::filesystem;
#else
	namespace fs = std::experimental::filesystem;
#endif

/**
 * \is skip the folder
 */
enum is_skip_folder { N = 0, Y = 1 };

/**
 * \forward declaration for image_instance
 */
struct image_instance;

class IO : MATRICE_STD(ios) {
	using _Mybase = MATRICE_STD(ios);
public:
	using FileMap = std::map<std::string, std::string>;
	using FileInfo = std::vector<FileMap>;
	using _Mybase::in;
	using _Mybase::out;
	using _Mybase::app;
	using _Mybase::ate;
	/**
	 * \extracts folder(s) for a given absolute path.
	 */
	template<size_t _Nsubfolder = 0> 
	struct Dir_  {
		Dir_() 
			: m_path(fs::current_path().string()){ 
			_Init(); 
		}
		Dir_(std::string&& _path) 
			: m_path(_path){ 
			_Init(); 
		}
		Dir_(fs::path&& _path) 
			: m_path(_path.string()){ 
			_Init(); 
		}

		MATRICE_HOST_INL decltype(auto)operator()(size_t i, size_t j) const {
			return std::string(m_path + m_subpaths[i] + m_names[i][j]);
		}
		MATRICE_HOST_INL decltype(auto)operator()(size_t i) const {
			return std::string(m_path + m_names[0][i]);
		}
		MATRICE_HOST_INL decltype(auto)names(size_t i = 0) {
			using names_type = typename decltype(m_names)::value_type;
			return forward<names_type>(m_names)[i];
		}
		MATRICE_HOST_INL decltype(auto)names(size_t i = 0) const {
			using names_type = typename decltype(m_names)::value_type;
			return forward<names_type>(m_names)[i];
		}
		MATRICE_HOST_INL decltype(auto)count(size_t _Idx = 0) const {
			return m_names[_Idx].size(); 
		}

		MATRICE_HOST_FINL decltype(auto)path() noexcept { 
			return forward<decltype(m_path)>(m_path); 
		}
		MATRICE_HOST_FINL decltype(auto)path() const noexcept { 
			return forward<decltype(m_path)>(m_path);
		}
	private:
		MATRICE_HOST_INL void _Init() {
			m_subpaths.clear();
			fs::directory_iterator _End; size_t _Cnt = 0;
			for (decltype(_End) _Begin(m_path); _Begin != _End; ++_Begin) {
				if (fs::is_directory(_Begin->status())) {
					auto _Path = _Begin->path().string().substr(m_path.size());
					m_subpaths.emplace_back(_Path);
					++_Cnt;
				}
			}
			if constexpr (_Nsubfolder == 0) 
				m_subpaths.emplace_back("");
			if (_Nsubfolder != 0 && _Nsubfolder > m_subpaths.size())
				DGELOM_ERROR("No enough subfolders");
			
			_Cnt = 0; 
			m_names.resize(_Nsubfolder == 0 ? 1 : min(_Nsubfolder, m_subpaths.size()));
			for_each(m_names, [&](auto& _List) {
				_List = filenames<Y>(m_path + m_subpaths[_Cnt++]);
			});
		}

		std::string m_path;
		std::vector<std::string> m_subpaths;
		std::vector<std::vector<std::string>> m_names;
	};

	/**
	 * \csv writer
	 */
	class CSV {
	public:
		CSV(std::string&& _Fname, const std::string& _Delm = ",")
			:_Mypath(_Fname), _My_delimeter(_Delm){
		}
		CSV(const fs::path& _Path = {}, const std::string& _Delm = ",")
			:_Mypath(_Path), _My_delimeter(_Delm) {
		}
		CSV(fs::path&& _Path, const std::string& _Delm = ",")
			:_Mypath(_Path), _My_delimeter(_Delm) {
		}

		CSV& reset(const fs::path& _Path) noexcept {
			_Mypath = _Path;
			return (*this);
		}

		MATRICE_HOST_FINL void open(int _Mode) {
			if (_Mode & out || _Mode & app) {
				if (!fs::exists(_Mypath.parent_path())) {
					fs::create_directories(_Mypath.parent_path());
				}
			}
			_My_file.open(_Mypath.string(), _Mode);
		}
		MATRICE_HOST_FINL void close() noexcept {
			_My_file.close(); 
		}

		/// <summary>
		/// \brief Add a dummy line, with range [_First, _Last).
		/// </summary>
		/// <typeparam name="_It">Forward iterator type</typeparam>
		/// <param name="_First">: Begin of the range.</param>
		/// <param name="_Last">: End of the data range.</param>
		/// <returns>Current line number.</returns>
		template<typename _It> 
		MATRICE_HOST_INL size_t append(_It _First, _It _Last) {
			for (; _First != _Last;) {
				_My_file << *_First;
				if (++_First != _Last) 
					_My_file << _My_delimeter;
			}
			_My_file << "\n";

			return (_My_linecnt++);
		}

		/// <summary>
		/// \brief Add a labeled line, with range [_First, _Last).
		/// </summary>
		/// <typeparam name="_It">Forward iterator type</typeparam>
		/// <param name="_Label">: Line label string.</param>
		/// <param name="_First">: Begin of the range.</param>
		/// <param name="_Last">: End of the data range.</param>
		/// <returns>Current line number.</returns>
		template<typename _It>
		MATRICE_HOST_INL size_t append(std::string&& _Label, _It _First, _It _Last) {
			_My_file << _Label << _My_delimeter;
			for (; _First != _Last;) {
				_My_file << *_First;
				if (++_First != _Last)
					_My_file << _My_delimeter;
			}
			_My_file << "\n";

			return (_My_linecnt++);
		}

		MATRICE_HOST_INL size_t append(std::string&& _Label) {
			_My_file << _Label << "\n";
			return (_My_linecnt++);
		}

		/// <summary>
		/// \brief Get count of lines written. 
		/// </summary>
		/// <returns>Number of lines written</returns>
		MATRICE_HOST_INL size_t count() const noexcept {
			return (_My_linecnt);
		}

	private:
		fs::path _Mypath;
		std::string _My_delimeter = ",";
		std::size_t _My_linecnt = 0;
		std::fstream _My_file;
	};

	IO() noexcept {};

	static MATRICE_HOST_FINL auto workspace() noexcept { 
		return fs::current_path(); 
	}

	template<is_skip_folder _Fopt = is_skip_folder::Y>
	static MATRICE_HOST_INL auto filenames(std::string&& dir){
		std::vector<std::string> _Ret;

		fs::path _Path(dir);
		if (!fs::exists(_Path)) return std::move(_Ret);

		fs::directory_iterator _End;
		for (decltype(_End) _It(_Path); _It != _End; ++_It) {
			if (fs::is_regular_file(_It->status()))
				_Ret.push_back(_It->path().string().substr(dir.size()));

			if constexpr (_Fopt == is_skip_folder::N)
				if (fs::is_directory(_It->status())) {
					auto _Inner = filenames<_Fopt>(_It->path().string());
					_Ret.insert(_Ret.end(), _Inner.begin(), _Inner.end());
				}
		}

		return std::move(_Ret);
	}

	/**
	 *\brief Read data from a file.
	 *\param <_Op> The operator to define the data reader.
	 */
	template<typename _Op> 
	static MATRICE_HOST_INL auto read(std::string&& _path, _Op _op) {
		DGELOM_CHECK(fs::exists(fs::path(_path)), 
			"'" + _path + "' does not exist!");

		std::ifstream _Fin(_path);
		DGELOM_CHECK(_Fin.is_open(), 
			"Fail to open the file in " + _path);
			
		return _op(_Fin);
	}
	template<typename _Op>
	static MATRICE_HOST_INL auto read(fs::path&& _path, _Op _op) {
		return read(_path.string(), _op);
	}

	/**
	 *\brief Read data from a file
	 *\param <_N0, _N1> 1-based indices of the begin and end columns to be read.
	 *\param '_Path' refers to the location of target file.
     *\param '_Skips' indicates how many lines to skip over.
	 */
	template<typename _Ty, diff_t _N0, diff_t _N1 = _N0>
	static MATRICE_HOST_INL auto read(std::string&& _Path, size_t _Skips = 0) {
		static_assert(_N0 <= _N1, "_N1 must be greater than or equal to _N0");
		using value_type = _Ty;
		DGELOM_CHECK(fs::exists(fs::path(_Path)),"'"+_Path+"' does not exist!");
		std::ifstream _Fin(_Path);
		DGELOM_CHECK(_Fin.is_open(),"Cannot open file in "+_Path);

		std::string _Line;
		for (const auto _Idx : range(0, _Skips)) {
			std::getline(_Fin, _Line);
		}
		
		if constexpr (_N1 - _N0 > 1) {
			std::array<value_type, -~(_N1 - _N0)> _Myline;
			std::vector<decltype(_Myline)> _Data;
			while (std::getline(_Fin, _Line)) {
				const auto _Res = split<value_type>(_Line, ',');
				if (_Res.size() >= _N1-1) {
					for (auto _Idx = ~- _N0; _Idx < _N1; ++_Idx) {
						_Myline[-~(_Idx - _N0)] = _Res[_Idx];
					}
					_Data.push_back(_Myline);
				}
			}
			return forward<decltype(_Data)>(_Data);
		}
		else {
			std::vector<value_type> _Data;
			while (std::getline(_Fin, _Line)) {
				auto _Res = split<value_type>(_Line, ',');
				_Data.push_back(_Res[_N0-1]);
			}
			return forward<decltype(_Data)>(_Data);
		}
	}
	template<typename _Ty, diff_t _N0, diff_t _N1 = _N0>
	static MATRICE_HOST_INL auto read(fs::path&& _Path, size_t _Skips = 0) {
		return read<_Ty, _N0, _N1>(_Path.string(), _Skips);
	}

	/**
	 *\brief Read data from a file with a given path.
	 *\param '_Path' refers to the location of target file.
     *\param '_Skips' indicates how many lines to skip over.
	 */
	template<typename _Ty>
	static MATRICE_HOST_INL auto read(std::string&& _Path, size_t _Skips = 0) {
		using value_type = _Ty;
		DGELOM_CHECK(fs::exists(fs::path(_Path)), 
			"'" + _Path + "' does not exist!");
		std::ifstream _Fin(_Path);
		DGELOM_CHECK(_Fin.is_open(), 
			"Cannot open file in " + _Path);

		std::string _Line;
		for (const auto _Idx : range(0, _Skips)) {
			std::getline(_Fin, _Line);
		}
		std::vector<std::vector<value_type>> _Data;
		while (std::getline(_Fin, _Line)) {
			if(!_Line.empty())
				_Data.push_back(split<value_type>(_Line, ','));
		}

		Matrix<value_type> _Res(_Data.size(),_Data.front().size());
		for (const auto& _Idx : range(0, _Res.rows())) {
			_Res.rview(_Idx) = _Data[_Idx].data();
		}

		return forward<decltype(_Res)>(_Res);
	}
	template<typename _Ty>
	static MATRICE_HOST_INL auto read(fs::path&& _Path, size_t _Skips = 0) {
		return read<_Ty>(_Path.string(), _Skips);
	}

	/**
	 *\brief Single-stream writter to output data to a file in '_path'.
	 */
	template<typename _Op> 
	static MATRICE_HOST_FINL void write(std::string&& _path, _Op _op) {
		std::ofstream _Fout(_path);
		DGELOM_CHECK(_Fout.is_open(), 
			"Fail to open file: " + _path);
		_Fout.setf(std::ios::fixed, std::ios::floatfield);
		_op(forward<std::ofstream>(_Fout));
		_Fout.close();
	}
	/**
	 * \brief Double-stream output, where '_op' should be a functor that accepts a pair of std::ofstream instances.
	 * \Example:
	    write("C:/data/", {"f1.txt", "f2.txt"}, [&](auto&& fs1, auto&& fs2){ ... })
	 */
	template<typename _Op>
	static MATRICE_HOST_FINL void write(std::string&& _path, std::initializer_list<std::string> _fnames, _Op&& _op) {
		std::cout << "Files are saved in folder: " << _path << std::endl;
		const auto _N = _fnames.begin();
		std::ofstream _O1(_path + _N[0]), _O2(_path + _N[1]);
		DGELOM_CHECK(_O1.is_open(), "Fail to open file: " + _path + _N[0]);
		DGELOM_CHECK(_O2.is_open(), "Fail to open file: " + _path + _N[1]);
		_O1.setf(std::ios::fixed, std::ios::floatfield);
		_O2.setf(std::ios::fixed, std::ios::floatfield);

		_op(forward<std::ofstream>(_O1), forward<std::ofstream>(_O2));

		_O1.close(), _O2.close();
	}
	// \multi-path write
	template<typename _Op> static
	MATRICE_HOST_FINL void write(const std::initializer_list<std::string> _paths, _Op _op)
	{
		for_each(_paths, [&_op](const auto& path) {write(path, _op); });
	}

	// \split a string with the given token
	template<typename _Ty = std::string> static 
	MATRICE_HOST_FINL auto split(const std::string& _str, char _tok) {
		MATRICE_USE_STD(vector); vector<_Ty> _Res;
		const auto _String = std::string(1, _tok) + _str;
		auto _Pos = _String.begin(), _End = _String.end();
		if (_String.back() == _tok) --_End;

		while (_Pos != _End) {
			auto _Pos_last = _Pos+1;
			_Pos = std::find(_Pos_last, _End, _tok);
			MATRICE_USE_STD(distance);
			const auto _Sub = _String.substr(
				distance(_String.begin(), _Pos_last), distance(_Pos_last, _Pos));
			if (!_Sub.empty()) {
				_Res.push_back(stonv<_Ty>(_Sub));
			}
		}

		return (_Res);
	}
	template<typename _Ty = std::string> static
	MATRICE_HOST_FINL auto split(const std::string& _str) {
		char _tok = ',';
		if constexpr (std::is_arithmetic_v<_Ty>) {
			const auto _Tok = std::find_if(_str.begin(), _str.end(),
				[](const auto _Val) {
					return !(std::isdigit(_Val) || _Val == '.' || _Val == '-');
				});
			_tok = *_Tok;
		}

		MATRICE_USE_STD(vector); vector<_Ty> _Res;
		const auto _String = std::string(1, _tok) + _str;
		auto _Pos = _String.begin(), _End = _String.end();
		if (_String.back() == _tok) --_End;

		while (_Pos != _End) {
			auto _Pos_last = _Pos + 1;
			_Pos = std::find(_Pos_last, _End, _tok);
			MATRICE_USE_STD(distance);
			_Res.push_back(stonv<_Ty>(_String.substr(distance(_String.begin(), _Pos_last),
				distance(_Pos_last, _Pos))));
		}

		return (_Res);
	}

};

/// <summary>
/// \brief Concate a file path and a file name.
/// </summary>
/// <param name="path">File path</param>
/// <param name="fname">File name</param>
MATRICE_HOST_INL auto concat(const fs::path& path, const std::string& fname) noexcept {
	auto _Path = path;
	return forward<fs::path>(_Path.append(fname));
}

// \read interface
template<typename _Ty, typename... _Args>
MATRICE_HOST_INL auto read(_Args...args){
	try { return IO::read<_Ty>(args...); } 
	catch (std::exception& _E) { throw _E; } 
};
// \write interface
template<typename... _Args> 
MATRICE_HOST_INL auto write(_Args...args)->decltype(IO::write(args...)){
	try { return IO::write(args...); } 
	catch (std::exception& _E) { throw _E; } 
};
// \definite a path under current folder
template<typename _Str>
MATRICE_HOST_INL auto defpath(const _Str local) noexcept {
	return IO::workspace().append(str(local));
};

/// <summary>
/// \brief Print a given matrix
/// </summary>
/// <typeparam name="_Ty"></typeparam>
/// <param name="m">a matrix with type of dgelom::Matrix_</param>
template<typename _Ty, int _M, int _N>
void print(const Matrix_<_Ty, _M, _N>& m, const std::string& name = "Array") {
	std::cout << "[" << name <<", dim:" << m.rows() << "x" << m.cols() << "\n";
	for (auto _It = m.rwbegin(); _It != m.rwend(); ++_It) {
		std::cout << " ";
		for (auto _Vl : _It)
			std::cout << std::setiosflags(std::ios::left) 
			<< std::setprecision(8) << std::setw(15)
			<< _Vl;
		std::cout << "\n";
	}
	std::cout << "]\n";
}
template<typename _Ty> requires is_scalar_v<_Ty>
void print(const _Ty& m, const std::string& name = "Scalar") {
	std::cout << "[" << name << ", dim:" << 1 << "x" << 1 << "\n";
	std::cout << " " << std::setiosflags(std::ios::left)
			<< std::setprecision(8) << std::setw(15) << _Ty(m);
	std::cout << "\n]\n";
}

/// <summary>
/// \brief CLASS, string_helper
/// </summary>
class string_helper MATRICE_NONHERITABLE {
	using basic_value_type = char;
public:
	using value_type = std::string;
	using const_value_type = std::add_const_t<value_type>;

	/**
	 *\brief Split a given string with a token 
	 *\param '_Str' a string to be splitted;
	 *\param '_Tok' token.
	 *\note: the template is used to specify the returned type.
	 *\example:
	   auto _Str1 = string_helper::value_type("dog,car,plane");
	   auto _Items1 = string_helper::split(_Str1,',');//{dog, car, plane}
	   auto _Str2 = string_helper::value_type("2,3,5");
	   auto _Items2 = string_helper::split<float>(_Str2,',');//{2.f, 3.f, 5.f}
	 */
	template<typename _Ty = value_type> static MATRICE_HOST_INL
	std::vector<_Ty> split(const value_type& _Str, basic_value_type _Tok) {
		std::vector<_Ty> _Res;
		const_value_type _String = value_type(1, _Tok) + _Str;

		auto _Pos = _String.begin();
		while (_Pos != _String.end()) {
			auto _Pos_last = ++_Pos;
			_Pos = std::find(_Pos_last, _String.end(), _Tok);
			_Res.push_back(dgelom::stonv<_Ty>(_String.substr(
				std::distance(_String.begin(), _Pos_last), 
				std::distance(_Pos_last, _Pos)
			)));
		}
		return forward<decltype(_Res)>(_Res);
	}

};

namespace io {
using path_t = MATRICE_STD(filesystem::path);
enum class fmt {
	txt,
	csv,
	xml,
	jsn
};
_DETAIL_BEGIN

struct folder_tag {};
struct loader_tag {
	struct tiff {};
	struct png {};
	struct bmp {};
	struct jpeg {};
	struct xml {};
};
struct file_tag {};

/**
 * \file or folder name collector
 */
template<typename _Tag> struct _Collector {};

/**
 * \build file directory according to a given path
 */
template<typename _Tag = folder_tag> class _Dir_impl {};

/**
 * \load data file from given directory
 */
template<typename _Ty, class _Tag> class _Data_loader_impl{};

/**
 * \file loader: tiff
 */
template<typename _Ty, class _Tag> struct _Loader_impl {};

/**
 * \save data to a file with a given format '_Fmt'.
 */
template<fmt _Fmt> void _Save(const auto& _Data, path_t _Path) {
#ifdef MATRICE_DEBUG
	DGELOM_CHECK(fs::exists(_Path.parent_path()),
		"The given folder path is not exist.");
#endif
	if constexpr (_Fmt == fmt::csv)
		_Path.replace_extension("csv");
	if constexpr (_Fmt == fmt::txt)
		_Path.replace_extension("txt");

	IO::CSV csv(_Path);
	csv.open(IO::out);
	if constexpr (is_matrix_v<remove_all_t<decltype(_Data)>>) {
		for (auto _It = _Data.rwbegin(); _It < _Data.rwend(); ++_It) {
			csv.append(_It.begin(), _It.end());
		}
	}
	csv.close();
}
_DETAIL_END


/**
 *\brief <func> imread: read image file </func>
 *\param [path] the path of the image to be read.
 */
template<typename _Ty = uint8_t, class _Pth = std::string> 
MATRICE_HOST_INL auto imread(const _Pth path);

} DGE_MATRICE_END
#include "inline\_directory.inl"
#include "inline\_image.inl"
#include "inline\_data_loader.inl"

DGE_MATRICE_BEGIN namespace io {
using directory = detail::_Dir_impl<detail::folder_tag>;
using folder_collector = detail::_Collector<detail::folder_tag>;
using file_collector = detail::_Collector<detail::file_tag>;

template<typename _Ty = float_t>
using data_loader = detail::_Data_loader_impl<_Ty, detail::loader_tag>;
using data_loader_uint8 = data_loader<uint8_t>;
using data_loader_uint32 = data_loader<uint32_t>;
using data_loader_f32 = data_loader<float_t>;

template<typename _Ty = uint8_t>
using tiff = detail::_Loader_impl<_Ty, detail::loader_tag::tiff>;

template<size_t _N, typename _Cont>
MATRICE_HOST_FINL auto serial(const _Cont& _L) {
	static_assert(_N>=1, "_N must be no less than 1.");
	DGELOM_CHECK(_N<=_L.size(), "The size _N being serialized over range of _L.");
	return tuple_n<_N - 1>::_(_L.data());
}
template<class _Op>
decltype(auto) make_loader(std::string&& path, _Op&& loader)noexcept {
	using loader_type = data_loader<typename _Op::value_type>;
	return loader_type(directory{ path, path_t() }, loader);
}
template<class _Op>
decltype(auto) make_loader(const path_t& path, _Op&& loader)noexcept {
	using loader_type = data_loader<typename _Op::value_type>;
	return loader_type(directory(path), loader);
}
template<class _Op>
decltype(auto) make_loader(path_t&& path, _Op&& loader)noexcept {
	using loader_type = data_loader<typename _Op::value_type>;
	return loader_type(directory(std::move(path)), loader);
}

template<fmt _Fmt>
MATRICE_HOST_INL void save(auto data, path_t path) noexcept {
	return detail::_Save<_Fmt>(data, path);
}
}

using io::imread;
using io::serial;
using io::make_loader;

DGE_MATRICE_END