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
#include <optional>
#include "_base.h"
#include "../../private/_iterator.h"
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
		size_t max_level = 8;
		_Ty tolerance = _Ty(1.0e-8);
		_Ty min_fill = 0.5;
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
		return (b == zero<_Ty> ? zero<_Ty> : a / b);
	}

	template<typename _Ty>
	MATRICE_HOST_INL static array_n_type<_Ty> _From(const _Ty* p) {
		array_n_type<_Ty> _Ret;
		std::copy(p, p + _Ret.size(), _Ret.begin());
		return (_Ret);
	}
};

// N-dimensional grid iterator (nested loop with variable depth).
template <unsigned _N = 2>
class grid_iterator : _MBA_base<_N> {
	using index_type = typename _MBA_base<_N>::template array_n_type<size_t>;
public:
	/**
	 *\brief construct from dims
	 *\param [dims] grid dimensions
	 */
	explicit grid_iterator(const index_type &dims) noexcept
		: _Dims(dims), _Idx{ 0 }, _Is_end(_Idx == _Dims), _Pos(0) { 
	}
	/**
	 *\brief construct from a given dim
	 *\param [dim] dimension
	 */
	explicit grid_iterator(size_t dim) noexcept
		: _Idx{ 0 }, _Pos(0), _Is_end(dim == 0) {
		for (auto &v : _Dims) v = dim;
	}

	MATRICE_HOST_INL size_t operator[](size_t d) const { return _Idx[d]; }
	MATRICE_HOST_INL index_type& operator*() { return _Idx; }
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
	return (_Left);
}

template <typename T, size_t N, typename C>
std::array<T, N> operator-(std::array<T, N> _Left, C b) {
	for (auto &v : _Left) v -= b;
	return (_Left);
}

template <typename T, size_t N, typename C>
std::array<T, N> operator*(std::array<T, N> _Left, C b) {
	for (auto &v : _Left) v *= b;
	return (_Left);
}

template <typename _Ty, unsigned _N>
class control_lattice : public _MBA_base<_N> {
	using _Mybase = _MBA_base<_N>;
public:
	using _Mybase::dim;
	using value_type = _Ty;
	using point_type = typename _Mybase::template array_n_type<value_type>;
	using index_type = typename _Mybase::template array_n_type<size_t>;
	using options_type = typename _Mybase::template options_type<value_type>;
	using functor_type = typename _Mybase::template functor<value_type>;

	virtual ~control_lattice() {}

	virtual value_type operator()(const point_type &p) const = 0;

	virtual void report(std::ostream&) const = 0;

	template <typename _Cont>
	MATRICE_HOST_INL value_type residual(_Cont& data) const {
		auto _Res = zero<value_type>;
		for (ptrdiff_t i = 0; i < data.rows(); ++i) {
			const auto p = _Mybase::_From(data[i]);
			_Res = max(_Res, abs(data[i][dim] -= (*this)(p)));
		}
		return _Res;
	}
};

template <typename _Ty, unsigned _N>
class initial_approximation : public control_lattice<_Ty, _N> {
	using _Mybase = control_lattice<_Ty, _N>;
public:
	using typename _Mybase::value_type;
	using typename _Mybase::point_type;
	initial_approximation(typename _Mybase::functor_type op): _Op(op) {}

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
	
	template <typename _Cont>
	control_lattice_dense(const _Cont& data, const options_type& opts) noexcept
		: cmin(opts.min), cmax(opts.max), grid(opts.grid)
	{
		for (unsigned i = 0; i < dim; ++i) {
			hinv[i] = (grid[i] - 1) / (cmax[i] - cmin[i]);
			cmin[i] -= 1 / hinv[i];
			grid[i] += 2;
		}

		multi_array<value_type, dim> delta(grid);
		multi_array<value_type, dim> omega(grid);
		phi.resize(grid);

		auto n = data.rows();
		auto m = phi.size();
#pragma omp parallel
		{
			multi_array<value_type, dim> t_delta(grid);
			multi_array<value_type, dim> t_omega(grid);

#pragma omp for
			for (ptrdiff_t l = 0; l < n; ++l) {
				auto p = _MBA_base<dim>::_From(data[l]);

				if (!_MBA_base<>::_Boxed(opts.min, p, opts.max)) continue;

				index_type i; point_type s;
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

				const auto v = data[l][dim];
				for (grid_iterator<dim> d(4); d; ++d) {
					value_type w1 = w[d.pos()];
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
		for (auto i = 0; i < m; ++i) phi[i] = safe_div(delta[i], omega[i]);
	}

	MATRICE_HOST_INL value_type operator()(const point_type &p) const {
		index_type i; point_type s;

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

	MATRICE_HOST_INL void append_refined(const control_lattice_dense &r) {
		static const std::array<value_type, 5> s{0.125, 0.500, 0.750, 0.500, 0.125};

		for (grid_iterator<dim> i(r.grid); i; ++i) {
			value_type f = r.phi(*i);

			if (f == zero<value_type>) continue;

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

				auto c = one<value_type>;
				for (unsigned k = 0; k < dim; ++k) c *= s[d[k]];

				phi(j) += f * c;
			}
		}
	}

	MATRICE_HOST_INL value_type fill_ratio() const {
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

	template <typename _Cont>
	control_lattice_sparse(const _Cont& data, const options_type& opts) noexcept
		: cmin(opts.min), cmax(opts.max), grid(opts.grid)
	{
		for (unsigned i = 0; i < dim; ++i) {
			hinv[i] = (grid[i] - 1) / (cmax[i] - cmin[i]);
			cmin[i] -= 1 / hinv[i];
			grid[i] += 2;
		}

		std::map<index_type, two_doubles> dw;
		for (auto idx = 0; idx < data.rows(); ++idx) {
			auto p = _MBA_base<dim>::_From(data[idx]);

			if (!_MBA_base<>::_Boxed(opts.min, p, opts.max)) continue;

			index_type i; point_type s;

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

			const auto v = data[idx][dim];
			for (grid_iterator<dim> d(4); d; ++d) {
				value_type w1 = w[d.pos()];
				value_type w2 = w1 * w1;
				value_type phi = v * w1 / sum_w2;

				auto ii = i + *d;
				dw[ii][0] += w2 * phi;
				dw[ii][1] += w2;
			}
		}

		phi.insert(make_transform_iter(dw.begin(), _Delta_over_omega),
			        make_transform_iter(dw.end(),   _Delta_over_omega));
	}

	MATRICE_HOST_INL value_type operator()(const point_type &p) const {
		index_type i; point_type s;

		for (unsigned d = 0; d < dim; ++d) {
			value_type u = (p[d] - cmin[d]) * hinv[d];
			i[d] = floor(u) - 1;
			s[d] = u - floor(u);
		}

		value_type f = 0;
		for (grid_iterator<dim> d(4); d; ++d) {
			value_type w = 1.0;
			for (unsigned k = 0; k < dim; ++k) 
				w *= _MBA_base<dim>::_Bicubic_val(d[k], s[k]);

			f += w * _Get_phi(i + (*d));
		}

		return f;
	}

	void report(std::ostream &os) const {
		std::ios_base::fmtflags ff(os.flags());
		const auto fp = os.precision();

		size_t grid_size = grid[0];

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

	MATRICE_HOST_INL value_type _Get_phi(const index_type &i) const {
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

/**
 *\brief MBA Multilevel Bicubic-spline Approximation (interpolation)
 *\param <_Ty> data type that MBA uses, 
			<_N> data dimension, which can be 1, 2, or 3. 
 *\note The input scattred data should be organized as
		  MBA<>::container _Data = {
		  x1, y1, f(x1,y1),
		  x2, y2, f(x2,y2),
				... ...
		  xn, yn, f(xn,yn)} 
		  where _N = 2 for this example. The last column of _Data 
		  will be replaced by interpolation residuals.
 */
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
	using container = Matrix<value_type>;
	using point_type = typename _Mybase::template array_n_type<value_type>;
	using index_type = typename _Mybase::template array_n_type<size_t>;
	using options_type = typename _Mybase::template options_type<value_type>;
	using functor_type = typename _Mybase::template functor<value_type>;

	MBA(const container& data, options_type& opts, 
		functor_type&& init_op = functor_type()) {
		_Init(data, opts, init_op);
	}

	MATRICE_HOST_INL value_type operator()(const point_type &p) const {
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

	MATRICE_HOST_INL void _Init(const container& data, options_type& opts, functor_type initial)
	{
		auto& grid = opts.grid;
		for (int i = 0; i < dim; ++i)
			_Mybase::_Precond(grid[i] > 1, 
				"MBA: grid size in each dimension should be more than 1");

		const auto n = data.rows();
		value_type res, eps = 0.0;
		for (auto i = 0; i < n; ++i) eps = max(eps, abs(data[i][dim]));
		eps *= opts.tolerance;

		if (initial) {
			cl.push_back(std::make_shared<approx_lattice>(initial));
			res = cl.back()->residual(data);
			if (res <= eps) return;
		}

		size_t lev = 1;
		{
			auto psi = std::make_shared<dense_lattice>(data, opts);
			res = psi->residual(data);
			value_type fill = psi->fill_ratio();

			for (; (lev < opts.max_level) && (res > eps) && (fill > opts.min_fill); ++lev) {
				for (auto& x : grid) x = x * 2 - 1;

				auto f = std::make_shared<dense_lattice>(data, opts);
				res = f->residual(data);
				fill = f->fill_ratio();
				f->append_refined(*psi);

				psi.swap(f);
			}

			cl.push_back(psi);
		}
		for (; (lev < opts.max_level) && (res > eps); ++lev) {
			for (auto& x : grid) x = x * 2 - 1;

			cl.push_back(std::make_shared<sparse_lattice>(data, opts));

			res = cl.back()->residual(data);
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

DGE_MATRICE_BEGIN
/**
 *\brief 1-dimensional scattered data interpolation.
 *\param <_Ty> data type, must be float or double.
 */
template<typename _Ty>
using multilevel_bicubic_1d = algs::MBA<_Ty, 1>;

/**
 *\brief 2-dimensional scattered data interpolation.
 *\param <_Ty> data type, must be float or double.
 */
template<typename _Ty>
using multilevel_bicubic_2d = algs::MBA<_Ty, 2>;

/**
 *\brief 3-dimensional scattered data interpolation.
 *\param <_Ty> data type, must be float or double.
 */
template<typename _Ty>
using multilevel_bicubic_3d = algs::MBA<_Ty, 3>;
DGE_MATRICE_END