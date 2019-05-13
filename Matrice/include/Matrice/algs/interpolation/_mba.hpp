/*********************************************************************
This file is part of Matrice, an effcient and elegant C++ library.
Copyright(C) 2018-2019, Zhilong(Dgelom) Su, all rights reserved.

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
#include <map>
#include "_base.h"
#include "../../private/container/_multi_array.hpp"

MATRICE_ALGS_BEGIN
_DETAIL_BEGIN
template<int8_t _N = 2> struct _MBA_base {
	static constexpr auto dim = _N;
	template<typename _Ty> using array_n_type = std::array<_Ty, dim>;
	template<typename _Ty> using functor = std::function<_Ty(const array_n_type<_Ty>&)>;
	template<typename _Ty> struct options_type {
		array_n_type<_Ty> min, max;
		array_n_type<size_t> grid;
	};

	template<typename _Cond, typename _Msg>
	MATRICE_HOST_INL static void _Precond(const _Cond& cond, const _Msg& msg) {
		DGELOM_CHECK(static_cast<bool>(cond), msg);
	}

	// Value of k-th B-Spline basic function at t.
	template<typename _Ty>
	MATRICE_HOST_INL static _Ty _Bicubic_val(size_t k, _Ty t) {
		assert(0 <= t && t < 1); assert(k < 4);

		switch (k) {
		case 0: return (t * (t * (-t + 3) - 3) + 1) / 6;
		case 1: return (t * t * (3 * t - 6) + 4) / 6;
		case 2: return (t * (t * (-3 * t + 3) + 3) + 1) / 6;
		case 3: return t * t * t / 6;
		default: return 0;
		}
	}

	// Checks if p is between lo and hi
	template<typename _Ty>
	MATRICE_HOST_INL static bool _Boxed(const array_n_type<_Ty> &lo, const array_n_type<_Ty> &p, const array_n_type<_Ty> &hi) {
		for (auto i = 0; i < dim; ++i) {
			if (p[i] < lo[i] || p[i] > hi[i]) return false;
		}
		return true;
	}

	template<typename _Ty>
	MATRICE_HOST_INL static _Ty _Safe_divide(_Ty a, _Ty b) {
		return (b == 0.0 ? 0.0 : a / b);
	}

	template<typename _Ty>
	MATRICE_HOST_INL static array_n_type<_Ty> _From(const _Ty* p) {
		array_n_type<_Ty> _Ret;
		std::copy(_Ret.begin(), _Ret.end(), p);
		return (_Ret);
	}
	template<typename _Ty, MATRICE_ENABLE_IF(is_class_v<_Ty>)>
	MATRICE_HOST_INL static array_n_type<typename _Ty::value_type> _From(const _Ty& p) {
		array_n_type<typename _Ty::value_type> _Ret;
		std::copy(_Ret.begin(), _Ret.end(), &p[0]);
		return (_Ret);
	}

};

//template <int N> using point = std::array<default_type, N>;

//template <int N> using index = std::array<size_t, N>;

// N-dimensional dense matrix
//template <class T, int N>
//class multi_array {
//	static_assert(N > 0, "Wrong number of dimensions");
//
//public:
//	multi_array() {}
//	multi_array(index<N> n) { init(n); }
//
//	void resize(index<N> n) { init(n); }
//
//	size_t size() const { return buf.size(); }
//
//	T operator()(index<N> i) const { return buf[idx(i)]; }
//	T& operator()(index<N> i) { return buf[idx(i)]; }
//	T operator[](size_t i) const { return buf[i]; }
//	T& operator[](size_t i) { return buf[i]; }
//	const T* data() const { return buf.data(); }
//	T* data() { return buf.data(); }
//
//private:
//	std::array<int, N> stride;
//	std::vector<T>  buf;
//
//	void init(index<N> n) {
//		size_t s = 1;
//
//		for (int d = N - 1; d >= 0; --d) {
//			stride[d] = s;
//			s *= n[d];
//		}
//
//		buf.resize(s);
//	}
//
//	size_t idx(index<N> i) const {
//		size_t p = 0;
//		for (int d = 0; d < N; ++d)
//			p += stride[d] * i[d];
//		return p;
//	}
//};

// N-dimensional grid iterator (nested loop with variable depth).
template <unsigned _N = 2>
class grid_iterator : _MBA_base<_N> {
	using index_type = typename _MBA_base<_N>::template array_n_type<size_t>;
public:
	explicit grid_iterator(const index_type &dims) noexcept
		: _Dims(dims), _Idx{ 0 }, _Is_end(_Idx == _Dims), _Pos(0) { }

	explicit grid_iterator(size_t dim) noexcept
		: _Idx{ 0 }, _Pos(0), _Is_end(dim == 0) {
		for (auto &v : _Dims)  v = dim;
	}

	MATRICE_HOST_INL size_t operator[](size_t d) const { return _Idx[d]; }
	MATRICE_HOST_INL const index_type& operator*() const { return _Idx; }
	MATRICE_HOST_INL grid_iterator& operator++() {
		_Is_end = std::true_type::value;
		for (size_t d = _N; d--; ) {
			if (++_Idx[d] < _Dims[d]) {
				_Is_end = std::false_type::value;
				break;
			}
			_Idx[d] = 0;
		}
		++_Pos;
		return (*this);
	}
	MATRICE_HOST_INL size_t pos() const { return _Pos; }
	MATRICE_HOST_INL operator bool() const { return !_Is_end; }

private:
	size_t _Pos;
	index_type _Dims, _Idx;
	bool _Is_end = std::true_type::value;
};

template <typename T, size_t N>
std::array<T, N> operator+(std::array<T, N> _Left, const std::array<T, N> &_Right) {
	for (size_t i = 0; i < N; ++i) _Left[i] += _Right[i];
	return _Right;
}

template <typename T, size_t N, typename C>
std::array<T, N> operator-(std::array<T, N> a, C b) {
	for (auto &v : a) v -= b;
	return a;
}

template <typename T, size_t N, typename C>
std::array<T, N> operator*(std::array<T, N> a, C b) {
	for (auto &v : a) v *= b;
	return a;
}

template <typename _Ty, unsigned _N>
class control_lattice : _MBA_base<_N> {
	using _Mybase = _MBA_base<_N>;
public:
	using _Mybase::dim;
	using value_type = _Ty;
	using point_type = _Mybase::template array_n_type<value_type>;
	using index_type = _Mybase::template array_n_type<size_t>;
	using options_type = _Mybase::template options_type<value_type>;
	using functor_type = _Mybase::template functor<value_type>;

	virtual ~control_lattice() {}

	virtual value_type operator()(const point_type &p) const = 0;

	virtual void report(std::ostream&) const = 0;

	template <typename _Cont, class _Iter>
	MATRICE_HOST_INL value_type residual(const _Cont& coords, _Iter vbegin) const {
		auto res = zero<value_type>;

		auto p = coords.begin();
		auto v = vbegin;

		for (; p != coo_end; ++p, ++v) {
			(*v) -= (*this)(*p);
			res = max(res, abs(*v));
		}

		return res;
	}
};

template <typename _Ty, unsigned _N>
class initial_approximation : public control_lattice<_Ty, _N> {
	using _Mybase = control_lattice<_Ty, _N>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	initial_approximation(typename _Mybase::functor_type&& op): _Op(op) {}

	MATRICE_HOST_INL value_type operator()(const point_type& p)const { return _Op(p); }

	MATRICE_HOST_INL void report(std::ostream &os)const { os << "initial approximation"; }
private:
	typename _Mybase::functor_type _Op;
};

template <typename _Ty, unsigned _N>
class control_lattice_dense : public control_lattice<_Ty,_N> {
	using _Mybase = control_lattice<_Ty, _N>;
public:
	using _Mybase::dim;
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	using typename _Mybase::index_type;
	using typename _Mybase::options_type;
	
	template <typename _Cont, class _Iter>
	control_lattice_dense(const _Cont& coords, const _Iter& vbegin, 
		const options_type& opts) noexcept
		: cmin(opts.min), cmax(opts.max), grid(opts.grid)
	{
		for (unsigned i = 0; i < _N; ++i) {
			hinv[i] = (grid[i] - 1) / (cmax[i] - cmin[i]);
			cmin[i] -= 1 / hinv[i];
			grid[i] += 2;
		}

		multi_array<value_type, _N> delta(grid);
		multi_array<value_type, _N> omega(grid);
		phi.resize(grid);

		//std::fill(delta.data(), delta.data() + delta.size(), 0.0);
		//std::fill(omega.data(), omega.data() + omega.size(), 0.0);
		//std::fill(phi.data(), phi.data() + phi.size(), 0.0);

		auto n = coords.size();
		auto m = phi.size();

#pragma omp parallel
		{
			multi_array<value_type, dim> t_delta(grid);
			multi_array<value_type, dim> t_omega(grid);

			//std::fill(t_delta.data(), t_delta.data() + t_delta.size(), 0.0);
			//std::fill(t_omega.data(), t_omega.data() + t_omega.size(), 0.0);

#pragma omp for
			for (ptrdiff_t l = 0; l < n; ++l) {
				auto p = _MBA_base<>::_From(coords[l]);
				auto v = vbegin[l];

				if (!_MBA_base<>::_Boxed(opts.min, p, opts.max)) continue;

				index_type i;
				point_type s;

				for (unsigned d = 0; d < dim; ++d) {
					value_type u = (p[d] - cmin[d]) * hinv[d];
					i[d] = floor(u) - 1;
					s[d] = u - floor(u);
				}

				std::array<value_type, power_nm_v<4,dim>> w;
				value_type sum_w2 = 0.0;

				for (grid_iterator<dim> d(4); d; ++d) {
					value_type prod = 1.0;
					for (unsigned k = 0; k < dim; ++k) 
						prod *= _MBA_base<>::_Bicubic_val(d[k], s[k]);

					w[d.pos()] = prod;
					sum_w2 += prod * prod;
				}

				for (grid_iterator<dim> d(4); d; ++d) {
					value_type w1 = w[d.position()];
					value_type w2 = w1 * w1;
					value_type phi = v * w1 / sum_w2;

					auto j = i + (*d);

					t_delta(j) += w2 * phi;
					t_omega(j) += w2;
				}
			}

			{
				for (ptrdiff_t i = 0; i < m; ++i) {
#pragma omp atomic
					delta[i] += t_delta[i];
#pragma omp atomic
					omega[i] += t_omega[i];
				}
			}
		}

		for (ptrdiff_t i = 0; i < m; ++i) {
			phi[i] = safe_div(delta[i], omega[i]);
		}
	}

	value_type operator()(const point_type &p) const {
		index_type i;
		point_type s;

		for (unsigned d = 0; d < dim; ++d) {
			value_type u = (p[d] - cmin[d]) * hinv[d];
			i[d] = floor(u) - 1;
			s[d] = u - floor(u);
		}

		value_type f = 0;
		for (grid_iterator<dim> d(4); d; ++d) {
			value_type w = 1.0;
			for (unsigned k = 0; k < dim; ++k)
				w *= _MBA_base<>::_Bicubic_val(d[k], s[k]);
			f += w * phi(i + (*d));
		}

		return f;
	}

	void report(std::ostream &os) const {
		std::ios_base::fmtflags ff(os.flags());
		auto fp = os.precision();

		os << "dense  [" << grid[0];
		for (unsigned i = 1; i < dim; ++i)
			os << ", " << grid[i];
		os << "] (" << phi.size() * sizeof(value_type) << " bytes)";

		os.flags(ff);
		os.precision(fp);
	}

	void append_refined(const control_lattice_dense &r) {
		static const std::array<value_type, 5> s = {
			 0.125, 0.500, 0.750, 0.500, 0.125
		};

		for (grid_iterator<dim> i(r.grid); i; ++i) {
			value_type f = r.phi(*i);

			if (f == 0.0) continue;

			for (grid_iterator<dim> d(5); d; ++d) {
				index_type j;
				bool skip = false;
				for (unsigned k = 0; k < dim; ++k) {
					j[k] = 2 * i[k] + d[k] - 3;
					if (j[k] >= grid[k]) {
						skip = true;
						break;
					}
				}

				if (skip) continue;

				value_type c = 1.0;
				for (unsigned k = 0; k < dim; ++k) c *= s[d[k]];

				phi(j) += f * c;
			}
		}
	}

	value_type fill_ratio() const {
		size_t total = phi.size();
		size_t nonzeros = total - std::count(phi.begin(), phi.end(), zero<value_type>);

		return static_cast<value_type>(nonzeros) / total;
	}

private:
	point_type cmin, cmax, hinv;
	index_type grid;

	multi_array<value_type, dim> phi;
};

template<typename _Ty, unsigned _N>
class control_lattice_sparse : public control_lattice<_Ty, _N> {
	using _Mybase = control_lattice<_Ty, _N>;
public:
	using _Mybase::dim;
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	using typename _Mybase::index_type;
	using typename _Mybase::options_type;

	template <typename _Cont, class _Iter>
	control_lattice_sparse(const _Cont& coords, const _Iter& vbegin,
		const options_type& opts) noexcept
		: cmin(opts.min), cmax(opts.max), grid(opts.grid)
	{
		for (unsigned i = 0; i < dim; ++i) {
			hinv[i] = (grid[i] - 1) / (cmax[i] - cmin[i]);
			cmin[i] -= 1 / hinv[i];
			grid[i] += 2;
		}

		std::map<index<dim>, two_doubles> dw;
		auto v = vbegin;
		point_type p;
		for (auto idx = 0; idx < coords.size(); ++idx, ++v) {
			p = _MBA_base<>::_From(coords[idx]);

			if (!_MBA_base<>::_Boxed(opts.min, p, opts.max)) continue;

			index_type i;
			point_type s;

			for (unsigned d = 0; d < dim; ++d) {
				value_type u = (p[d] - cmin[d]) * hinv[d];
				i[d] = floor(u) - 1;
				s[d] = u - floor(u);
			}

			std::array<value_type, power_nm_v<4, dim>> w;
			value_type sum_w2 = 0.0;

			for (grid_iterator<dim> d(4); d; ++d) {
				value_type prod = 1.0;
				for (unsigned k = 0; k < dim; ++k)
					prod *= _MBA_base<>::_Bicubic_val(d[k], s[k]);

				w[d.pos()] = prod;
				sum_w2 += prod * prod;
			}

			for (grid_iterator<dim> d(4); d; ++d) {
				value_type w1 = w[d.pos()];
				value_type w2 = w1 * w1;
				value_type phi = (*v) * w1 / sum_w2;

				const auto ii = i + *d;
				dw[ii][0] += w2 * phi;
				dw[ii][1] += w2;
			}
		}

		phi.insert(
			make_transform_iter(dw.begin(), _Delta_over_omega),
			make_transform_iter(dw.end(),   _Delta_over_omega));
	}

	value_type operator()(const point_type &p) const {
		index_type i;
		point_type s;

		for (unsigned d = 0; d < dim; ++d) {
			value_type u = (p[d] - cmin[d]) * hinv[d];
			i[d] = floor(u) - 1;
			s[d] = u - floor(u);
		}

		value_type f = 0;

		for (grid_iterator<dim> d(4); d; ++d) {
			value_type w = 1.0;
			for (unsigned k = 0; k < dim; ++k) 
				w *= _MBA_base<>::_Bicubic_val(d[k], s[k]);

			f += w * _Get_phi(i + (*d));
		}

		return f;
	}

	void report(std::ostream &os) const {
		std::ios_base::fmtflags ff(os.flags());
		const auto fp = os.precision();

		const size_t grid_size = grid[0];

		os << "sparse [" << grid[0];
		for (unsigned i = 1; i < dim; ++i) {
			os << ", " << grid[i];
			grid_size *= grid[i];
		}

		const size_t bytes = phi.size() * sizeof(std::pair<index_type, value_type>);
		const size_t dense_bytes = grid_size * sizeof(value_type);

		value_type compression = static_cast<value_type>(bytes) / dense_bytes;
		os << "] (" << bytes << " bytes, compression: "
			<< std::fixed << std::setprecision(2) << compression << ")";

		os.flags(ff);
		os.precision(fp);
	}
private:
	point_type cmin, cmax, hinv;
	index_type grid;

	std::map<index_type, value_type> phi;

	typedef std::array<value_type, 2> two_doubles;

	static std::pair<index_type, value_type> _Delta_over_omega(const std::pair<index_type, two_doubles> &dw) {
		return std::make_pair(dw.first, safe_div(dw.second[0], dw.second[1]));
	}

	value_type _Get_phi(const index_type &i) const {
		typename decltype(phi)::const_iterator c = phi.find(i);
		return c == phi.end() ? 0.0 : c->second;
	}
};
_DETAIL_END

template <typename _Ty, unsigned _N>
class linear_approximation {
public:
	static constexpr auto dim = _N;
	using value_type = _Ty;
	using index_type = std::array<size_t, dim>;
	using point_type = std::array<value_type, dim>;
	template <class CooIter, class ValIter>
	linear_approximation(CooIter coo_begin, CooIter coo_end, ValIter val_begin)
	{
		size_t n = std::distance(coo_begin, coo_end);
		for (auto &v : C) v = 0.0;

		if (n <= dim) {
			// Not enough points to get a unique plane
			C[dim] = std::accumulate(val_begin, val_begin + n, 0.0) / n;
			return;
		}

		value_type A[dim + 1][dim + 1] = { 0 };

		CooIter p = coo_begin;
		ValIter v = val_begin;

		value_type sum_val = 0.0;

		// Solve least-squares problem to get approximation with a plane.
		for (; p != coo_end; ++p, ++v, ++n) {
			value_type x[dim + 1];
			for (unsigned i = 0; i < dim; ++i) x[i] = (*p)[i];
			x[dim] = 1.0;

			for (unsigned i = 0; i <= dim; ++i) {
				for (unsigned j = 0; j <= dim; ++j) {
					A[i][j] += x[i] * x[j];
				}
				C[i] += x[i] * (*v);
			}

			sum_val += (*v);
		}

		// Perform LU-factorization of A in-place
		for (unsigned k = 0; k <= dim; ++k) {
			if (A[k][k] == 0) goto singular;

			value_type d = 1 / A[k][k];
			for (unsigned i = k + 1; i <= dim; ++i) {
				A[i][k] *= d;
				for (unsigned j = k + 1; j <= dim; ++j)
					A[i][j] -= A[i][k] * A[k][j];
			}
			A[k][k] = d;
		}

		for (unsigned i = 0; i <= dim; ++i)
			if (A[i][i] == 0.0) goto singular;

		// Lower triangular solve:
		for (unsigned i = 0; i <= dim; ++i) {
			for (unsigned j = 0; j < i; ++j)
				C[i] -= A[i][j] * C[j];
		}

		// Upper triangular solve:
		for (unsigned i = dim + 1; i-- > 0; ) {
			for (unsigned j = i + 1; j <= dim; ++j)
				C[i] -= A[i][j] * C[j];
			C[i] *= A[i][i];
		}

		return;
	singular:
		for (unsigned i = 0; i < dim; ++i) C[i] = 0.0;
		C[dim] = sum_val / n;
	}

	value_type operator()(const point_type &p) const {
		value_type f = C[dim];

		for (unsigned i = 0; i < dim; ++i)
			f += C[i] * p[i];

		return f;
	}
private:
	std::array<value_type, dim+1> C;
};

template<typename _Ty = float, unsigned _N = 2>
class MBA : detail::_MBA_base<_N> {
	using _Mybase = detail::_MBA_base<_N>;
	using lattice = detail::control_lattice<_Ty, _N>;
	using approx_lattice = detail::initial_approximation<_Ty, _N>;
	using dense_lattice = detail::control_lattice_dense<_Ty, _N> ;
	using sparse_lattice = detail::control_lattice_sparse<_Ty, _N>;
public:
	using _Mybase::dim;
	using value_type = _Ty;
	using point_type = _Mybase::template array_n_type<value_type>;
	using index_type = _Mybase::template array_n_type<size_t>;
	using options_type = _Mybase::template options_type<value_type>;
	using functor_type = _Mybase::template functor<value_type>;

	template <class _Cont, class _Iter>
	MBA(const _Cont& coords, const _Iter& vbegin, const options_type& opts,
		unsigned max_levels = 8, value_type tol = 1e-8, value_type min_fill = 0.5,
		functor_type initial = functor_type()) {
		_Init(coords, vbegin, opts, max_levels, tol, min_fill, initial);
	}

	value_type operator()(const point_type &p) const {
		value_type f = 0.0;
		for (auto &psi : cl) f += (*psi)(p);
		return f;
	}

	friend std::ostream& operator<<(std::ostream &os, const MBA &other) {
		size_t level = 0;
		for (auto &psi : other.cl) {
			os << "level " << ++level << ": ";
			psi->report(os);
			os << std::endl;
		}
		return os;
	}

private:
	std::list<std::shared_ptr<lattice>> cl;

	template<class _Cont, class _Iter>
	void _Init(const _Cont& coords, const _Iter& vbegin, const options_type& opts,
		unsigned max_levels, value_type tol, value_type min_fill, functor_type initial)
	{
		for (int i = 0; i < dim; ++i)
			_Mybase::_Precond(grid[i] > 1, "MBA: grid size in each dimension should be more than 1");

		const auto n = coords.size();
		std::vector<value_type> val(vbegin, vbegin + n);

		value_type res, eps = 0.0;
		for (ptrdiff_t i = 0; i < n; ++i)
			eps = max(eps, abs(val[i]));
		eps *= tol;

		if (initial) {
			// Start with the given approximation.
			cl.push_back(std::make_shared<approx_lattice>(initial));
			res = cl.back()->residual(coords, val.begin());
			if (res <= eps) return;
		}

		size_t lev = 1;
		// Create dense head of the hierarchy.
		{
			auto psi = std::make_shared<dense_lattice>(coords, val.begin(), opts);

			res = psi->residual(coords, val.begin());
			const value_type fill = psi->fill_ratio();

			for (; (lev < max_levels) && (res > eps) && (fill > min_fill); ++lev) {
				grid = grid * 2ul - 1ul;

				auto f = std::make_shared<dense_lattice>(coords, val.begin(), opts);

				res = f->residual(coords, val.begin());
				fill = f->fill_ratio();

				f->append_refined(*psi);
				psi.swap(f);
			}

			cl.push_back(psi);
		}

		// Create sparse tail of the hierrchy.
		for (; (lev < max_levels) && (res > eps); ++lev) {
			grid = grid * 2ul - 1ul;

			cl.push_back(std::make_shared<sparse_lattice>(coords, val.begin(), opts));

			res = cl.back()->residual(coords val.begin());
		}
	}
};

/**
 *\brief multilevel bicubic spline interpolation interface
 *\param <_Ty> a scalar template type
 */
template<typename _Ty> 
class _Spline_interpolation<_Ty, _TAG mbicspl_tag> : public
	_Interpolation_base<_Spline_interpolation<_Ty, _TAG mbicspl_tag>> 
{
	static_assert(is_scalar_v<_Ty>, "template type _Ty must be a scalar.");
	using _Myt = _Spline_interpolation<_Ty, _TAG mbicspl_tag>;
	using _Mybase = _Interpolation_base<_Myt>;
public:
	using typename _Mybase::category;
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	using typename _Mybase::matrix_type;
	using _Mybase::_Interpolation_base;
	static constexpr auto dimension = category::dimension;

public:
	 
	template<typename _Cont, typename _Iter>
	_Spline_interpolation(const _Cont& coods, const _Iter& vals) {}

	MATRICE_HOST_FINL void _Coeff_impl();
};

MATRICE_ALGS_END