/*********************************************************************
	  About License Agreement of this file, see "../lic/license.h"
*********************************************************************/
#pragma once
#include <iosfwd>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <functional>
#include <type_traits>
#include <experimental/filesystem>
#include "../core/matrix.h"
#include "../core/vector.h"
#include "../util/utils.h"
#include "../util/genalgs.h"
#include "../util/_macro_conditions.h"
#include "../util/_exception.h"
#include "../private/_range.h"

DGE_MATRICE_BEGIN
/**
 * \for filesystem namespace in STD library
 */
namespace fs = std::experimental::filesystem;
enum is_skip_folder { N = 0, Y = 1 };

template<typename _Ty = uint8_t>
struct image_instance {
	using value_type = _Ty;
	using image_type = Matrix<value_type>;

	MATRICE_HOST_INL operator Matrix_<value_type,0,0>() noexcept {
		return (m_data);
	}
	MATRICE_HOST_INL image_instance& create(uint32_t nchs) noexcept {
		m_nchs = nchs;
		m_width = m_cols * m_nchs;
		m_data.create(m_rows, m_width, zero<value_type>);

		return (*this);
	}
	template<typename _Op>
	MATRICE_HOST_INL decltype(auto) operator()(_Op&& op) noexcept {
		return (op(m_data));
	}

	uint32_t m_rows = 0, m_cols = 0;
	uint32_t m_nchs = 1, m_width = 0;
	image_type m_data;
};

class IO : std::ios {
	using _Mybase = std::ios;
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
	template<size_t _Nsubfolder = 0> struct Dir_  {
		Dir_() : m_path(fs::current_path().string()) { _Init(); }
		Dir_(const std::string& _path) : m_path(_path) { _Init(); }

		MATRICE_HOST_INL const auto operator()(size_t i, size_t j) const {
			return std::string(m_path + m_subpaths[i] + m_names[i][j]);
		}
		MATRICE_HOST_INL const auto operator()(size_t i) const {
			return std::string(m_path + m_names[0][i]);
		}
		MATRICE_HOST_FINL decltype(auto) names(size_t i = 0) {
			using names_type = typename decltype(m_names)::value_type;
			return forward<names_type>(m_names)[i];
		}
		MATRICE_HOST_FINL const auto& names(size_t i = 0) const { 
			using names_type = typename decltype(m_names)::value_type;
			return forward<names_type>(m_names)[i];
		}
		MATRICE_HOST_FINL const size_t count(size_t _Idx = 0) const { 
			return m_names[_Idx].size(); 
		}

		MATRICE_HOST_FINL decltype(auto) path() noexcept { 
			return forward<decltype(m_path)>(m_path); 
		}
		MATRICE_HOST_FINL decltype(auto) path() const noexcept { 
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
			if constexpr (_Nsubfolder == 0) m_subpaths.emplace_back("");
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
		CSV(const std::string& _Fname, const std::string& _Delm = ",")
			:_My_filename(_Fname), _My_delimeter(_Delm){}

		MATRICE_HOST_FINL void open(int _Mode) {
			_My_file.open(_My_filename, _Mode);
		}
		MATRICE_HOST_FINL void close() noexcept {
			_My_file.close(); 
		}

		template<typename _It> 
		MATRICE_HOST_FINL size_t append(_It _First, _It _Last) {
			for (; _First != _Last;) {
				_My_file << *_First;
				if (++_First != _Last) 
					_My_file << _My_delimeter;
			}
			_My_file << "\n";

			return (_My_linecnt++);
		}

	private:
		std::string _My_filename, _My_delimeter;
		std::size_t _My_linecnt = 0;
		std::fstream _My_file;
	};

	IO() noexcept {};

	static MATRICE_HOST_FINL auto workspace() noexcept { 
		return fs::current_path(); 
	}

	template<is_skip_folder _Fopt = is_skip_folder::Y>
	static MATRICE_HOST_INL auto filenames(const std::string& dir)
	{
		std::vector<std::string> _Ret;

		fs::path _Path(dir);
		if (!fs::exists(_Path)) return std::move(_Ret);

		fs::directory_iterator _End;
		for (fs::directory_iterator _Begin(_Path); _Begin != _End; ++_Begin) {
			if (fs::is_regular_file(_Begin->status()))
				_Ret.emplace_back(_Begin->path().string().substr(dir.size()));

			if constexpr (_Fopt == is_skip_folder::N)
				if (fs::is_directory(_Begin->status())) {
					auto _Inner = filenames<_Fopt>(_Begin->path().string());
					_Ret.insert(_Ret.end(), _Inner.begin(), _Inner.end());
				}
		}

		return std::move(_Ret);
	}

	///<summary>
	//@brief: Template function to cvt. any value to std::string 
	//@author: Zhilong Su - Jan.10.2017 @SEU
	///</summary>
	template<typename _Ty> 
	static MATRICE_HOST_FINL auto strf(_Ty _val) noexcept {
		std::ostringstream strs; strs << _val; return strs.str();
	}
	// \TEMPLATE function to cvt. numeric value to std::string
	template<typename _Ty> 
	static MATRICE_HOST_FINL auto strfn(_Ty _val) noexcept {
#ifdef __CXX11__
		return (std::to_string(_val));
#elif
		return strf(_val);
#endif
	}

	template<typename _Op> 
	static MATRICE_HOST_INL auto read(std::string _path, _Op _op) {
		DGELOM_CHECK(fs::exists(fs::path(_path)), "'" + _path + "' does not exist!");

		std::ifstream _Fin(_path);

		DGELOM_CHECK(_Fin.is_open(), "Fail to open the file in " + _path);
			
		return _op(_Fin);
	}

	/**
	 *\brief read data from a file
	 *\param <_N0, _N1> are the 1-based indices of the begin and end columns to be read. 
	 */
	template<typename _Ty, size_t _N0, size_t _N1 = _N0>
	static MATRICE_HOST_INL auto read(std::string _Path, size_t _Skips = 0) {
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
				auto _Res = split<value_type>(_Line, ',');
				for (auto _Idx = ~- _N0; _Idx < _N1; ++_Idx) {
					_Myline[-~(_Idx - _N0)] = _Res[_Idx];
				}
				_Data.emplace_back(_Myline);
			}
			return forward<decltype(_Data)>(_Data);
		}
		else {
			std::vector<value_type> _Data;
			while (std::getline(_Fin, _Line)) {
				auto _Res = split<value_type>(_Line, ',');
				_Data.emplace_back(_Res[_N0-1]);
			}
			return forward<decltype(_Data)>(_Data);
		}
	}

	/**
	 *\brief read data from a file at _Path
	 *\param [_Path] refers to the location of target file
			   [_Skips] indicates how many lines to skip over
	 */
	template<typename _Ty>
	static MATRICE_HOST_INL auto read(std::string _Path, size_t _Skips = 0) {
		using value_type = _Ty;
		DGELOM_CHECK(fs::exists(fs::path(_Path)), "'" + _Path + "' does not exist!");
		std::ifstream _Fin(_Path);
		DGELOM_CHECK(_Fin.is_open(), "Cannot open file in " + _Path);

		std::string _Line;
		for (const auto _Idx : range(0, _Skips)) {
			std::getline(_Fin, _Line);
		}
		std::vector<std::vector<value_type>> _Data;
		while (std::getline(_Fin, _Line)) {
			_Data.push_back(split<value_type>(_Line, ','));
		}

		Matrix<value_type> _Res(_Data.size(),_Data.front().size());
		for (const auto& _Idx : range(0, _Res.rows())) {
			_Res.rview(_Idx) = _Data[_Idx].data();
		}

		return std::forward<decltype(_Res)>(_Res);
	}

	// \single-stream output
	template<typename _Op> 
	static MATRICE_HOST_FINL void write(std::string _path, _Op _op) {
		std::ofstream _Fout(_path);
		DGELOM_CHECK(_Fout.is_open(), "Fail to open file: " + _path);
		_Fout.setf(std::ios::fixed, std::ios::floatfield);
		_op(forward<std::ofstream>(_Fout));
		_Fout.close();
	}
	/**
	 * \double-stream output, where _op should be a functor accepts a pair of std::ofstream instances.
	 * \Example:
			write("C:/data/", {"f1.txt", "f2.txt"}, [&](auto&& fs1, auto&& fs2){ ... })
	 */
	template<typename _Op>
	static MATRICE_HOST_FINL void write(std::string _path, std::initializer_list<std::string> _fnames, _Op&& _op) {
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
		std::vector<_Ty> _Res;

		const auto _String = std::string(1, _tok) + _str;
		auto _Pos = _String.begin(), _End = _String.end();
		if (_String.back() == _tok) --_End;

		while (_Pos != _End) {
			auto _Pos_last = _Pos+1;
			_Pos = std::find(_Pos_last, _End, _tok);
			_Res.push_back(stonv<_Ty>(_String.substr(std::distance(_String.begin(), _Pos_last), std::distance(_Pos_last, _Pos))));
		}

		return (_Res);
	}

};

// \read interface
template<typename... _Args>
MATRICE_HOST_INL auto read(_Args...args)->decltype(IO::read(args...)){
	try { return IO::read(args...); } 
	catch (std::exception& _E) { throw _E; } 
};
// \write interface
template<typename... _Args> 
MATRICE_HOST_INL auto write(_Args...args)->decltype(IO::write(args...)){
	try { return IO::write(args...); } 
	catch (std::exception& _E) { throw _E; } 
};
// \definite a path under current folder
template<typename T>
MATRICE_HOST_INL decltype(auto) defpath(const T local) {
	return forward<std::string>(IO::workspace().string() + "\\" + IO::strf(local)); 
};

// \Class: std::string helper  
// \Coded by: dgelom su
class string_helper MATRICE_NONHERITABLE {
	using basic_value_type = char;
public:
	using value_type = std::string;
	using const_value_type = std::add_const_t<value_type>;

	template<typename _Ty = value_type> static
	/**
	 *\brief split a given string with a token 
	 *\param [_Str] string being splitted; [_Tok] token.
	 *\note: the template is used to specify the returned type.
	 *\example:
	   auto _Str1 = string_helper::value_type("dog,car,plane");
		auto _Items1 = string_helper::split(_Str1,',');//{dog, car, plane}
		auto _Str2 = string_helper::value_type("2,3,5");
		auto _Items2 = string_helper::split<float>(_Str2,',');//{2.f, 3.f, 5.f}
	 */
	MATRICE_HOST_INL auto split(const value_type& _Str, basic_value_type _Tok) {
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
_DETAIL_END


/**
 *\brief <func> imread: read image file </func>
 *\param [path] the path of the image to be read.
 */
template<typename _Ty = uint8_t, class _Pth = std::string> 
MATRICE_HOST_INL decltype(auto) imread(const _Pth path);

} DGE_MATRICE_END
#include "inline\_directory.inl"
#include "inline\_image.inl"
#include "inline\_data_loader.inl"

DGE_MATRICE_BEGIN namespace io {
using directory = detail::_Dir_impl<detail::folder_tag>;
using folder_collector = detail::_Collector<detail::folder_tag>;
using file_collector = detail::_Collector<detail::file_tag>;
using path_t = directory::path_type;

template<typename _Ty = std::float_t>
using data_loader = detail::_Data_loader_impl<_Ty, detail::loader_tag>;
using data_loader_uint8 = data_loader<uint8_t>;
using data_loader_uint32 = data_loader<uint32_t>;
using data_loader_f32 = data_loader<float_t>;

template<typename _Ty = float>
using tiff = detail::_Loader_impl<_Ty, detail::loader_tag::tiff>;

template<size_t _N, typename _Cont>
MATRICE_HOST_FINL auto serial(const _Cont& _L) {
	static_assert(_N>=1, "_N must be no less than 1.");
	DGELOM_CHECK(_N<=_L.size(), "The size _N being serialized over range of _L.");
	return tuple_n<_N - 1>::_(_L.data());
}
template<typename _Ty, class _Op>
decltype(auto) make_loader(std::string path, _Op&& loader)noexcept {
	using loader_type = io::data_loader<_Ty>;
	return loader_type(io::directory{ path, io::path_t() }, loader);
}
}

using io::imread;
using io::serial;
using io::make_loader;

DGE_MATRICE_END