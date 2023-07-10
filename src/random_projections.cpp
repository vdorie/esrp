#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include <Rconfig.h> // for HAVE_ALLOCA_H
#ifdef __GNUC__
// this covers gcc, clang, icc
# undef alloca
# define alloca(x) __builtin_alloca((x))
#elif defined(HAVE_ALLOCA_H)
// needed for native compilers on Solaris and AIX
# include <alloca.h>
#endif
#undef R_NO_REMAP

#include <cmath> // sqrt
#include <cstdint>
#include <cstring>

#include <algorithm> // std::clamp
#include <limits>
#include <chrono> // std::chrono::system_clock

#include <boost/random.hpp>
#include <boost/math/distributions.hpp>

namespace {

enum RandomType {
  RANDOM_TYPE_NORMAL,
  RANDOM_TYPE_DISCRETE
};

struct Options {
  unsigned int d;
  unsigned int output_length;
  unsigned int row_start;
  std::uint32_t seed;
  RandomType random_type;
};


struct NormalGeneratorInfo {
  boost::math::normal_distribution<double> normal_dist;
  double max_int;

  NormalGeneratorInfo() :
    normal_dist(0.0, 1.0), 
    max_int(static_cast<double>(std::numeric_limits<std::uint32_t>::max())) {
  }
};

struct DiscreteGeneratorInfo {
  std::uint32_t p[2];
};


template<typename T>
double generate_normal(T& generator, void* v_info) {
  NormalGeneratorInfo& info(*reinterpret_cast<NormalGeneratorInfo*>(v_info));
  return boost::math::quantile(info.normal_dist, std::clamp(static_cast<double>(generator()) / info.max_int, 0.001, 0.999));
}

template<typename T>
double generate_discrete(T& generator, void* v_info) {
  DiscreteGeneratorInfo& info(*reinterpret_cast<DiscreteGeneratorInfo*>(v_info));
  std::uint32_t value = generator();
  if (value < info.p[0]) return -1.0;
  if (value < info.p[1]) return 0.0;
  return 1.0;
}

}

extern "C" {

SEXP generate_random_projections(SEXP indices_expr, SEXP u_expr, SEXP k_expr, SEXP options_expr) {
  int* indices = INTEGER(indices_expr);
  
  const int** us = NULL;
  unsigned int num_input_vectors = u_expr != R_NilValue ? LENGTH(u_expr) : 0;
  if (num_input_vectors > 0) {
    us = reinterpret_cast<const int**>(alloca(LENGTH(u_expr) * sizeof(const int*)));

    for (unsigned int l = 0; l < num_input_vectors; ++l) {
      us[l] = const_cast<const int*>(INTEGER(VECTOR_ELT(u_expr, l)));
    }
  }
  unsigned int k = static_cast<unsigned int>(INTEGER(k_expr)[0]);
  if (k == 0)
    Rf_error("k must be greater than 0");
  
  Options options;
  options.d = 0;
  options.output_length = 0;
  options.row_start = 0;
  options.seed = 0;
  options.random_type = RANDOM_TYPE_NORMAL;

  bool seed_unset = true;
  
  if (options_expr == R_NilValue)
    Rf_error("option argument cannot be NULL");
  SEXP control_names_expr = Rf_getAttrib(options_expr, R_NamesSymbol);
  if (control_names_expr == R_NilValue)
    Rf_error("option argument must be a named list");
  
  for (unsigned int options_index = 0; options_index < LENGTH(options_expr); ++options_index) {
    if (std::strcmp(CHAR(STRING_ELT(control_names_expr, options_index)), "d") == 0) {
      options.d = static_cast<unsigned int>(INTEGER(VECTOR_ELT(options_expr, options_index))[0]);
    } else if (std::strcmp(CHAR(STRING_ELT(control_names_expr, options_index)), "output_length") == 0) {
      options.output_length = static_cast<unsigned int>(INTEGER(VECTOR_ELT(options_expr, options_index))[0]);
    } else if (std::strcmp(CHAR(STRING_ELT(control_names_expr, options_index)), "row_start") == 0) {
      options.row_start = static_cast<unsigned int>(INTEGER(VECTOR_ELT(options_expr, options_index))[0]);
    } else if (std::strcmp(CHAR(STRING_ELT(control_names_expr, options_index)), "seed") == 0) {
      options.seed = static_cast<std::uint32_t>(INTEGER(VECTOR_ELT(options_expr, options_index))[0]);
      seed_unset = false;
    } else if (std::strcmp(CHAR(STRING_ELT(control_names_expr, options_index)), "random_type") == 0) {
      const char* random_type_name = CHAR(STRING_ELT(VECTOR_ELT(options_expr, options_index), 0));
      if (std::strcmp(random_type_name, "normal") == 0) {
        options.random_type = RANDOM_TYPE_NORMAL;
      } else if (std::strcmp(random_type_name, "discrete") == 0) {
        options.random_type = RANDOM_TYPE_DISCRETE;
      } else {
        Rf_error("unrecognized random type");
      }
    } else {
      Rf_error("unrecognized option name");
    }
  }
  ;
  if (options.output_length == 0)
    options.output_length = k;
  if (options.random_type == RANDOM_TYPE_DISCRETE && options.d == 0)
    Rf_error("d must be set for random type discrete");
  
  unsigned int num_vectors = num_input_vectors == 0 ? 1 : num_input_vectors;
  double** vs = reinterpret_cast<double**>(alloca(num_vectors * sizeof(double*)));
  SEXP result_expr = PROTECT(Rf_allocVector(VECSXP, num_input_vectors == 0 ? 1 : num_input_vectors));
  for (unsigned int l = 0; l < num_vectors; ++l) {
    SEXP v_expr = PROTECT(Rf_allocVector(REALSXP, options.output_length));
    vs[l] = REAL(v_expr);
    SET_VECTOR_ELT(result_expr, l, v_expr);
    UNPROTECT(1);

    for (unsigned int j = 0; j < options.output_length; ++j) {
      vs[l][j] = 0.0;
    }
  }

  void* generator_info;
  double (*generate_random_number)(boost::random::mt19937& generator, void* info);

  if (options.random_type == RANDOM_TYPE_NORMAL) {
    generator_info = reinterpret_cast<void*>(new NormalGeneratorInfo());
    generate_random_number = &generate_normal;
  } else {
    DiscreteGeneratorInfo* d_info = new DiscreteGeneratorInfo();
    double sqrt_d = std::sqrt(static_cast<double>(options.d));
    double uint_max = static_cast<double>(std::numeric_limits<std::uint32_t>::max());
    d_info->p[0] = static_cast<std::uint32_t>(uint_max / (2.0 * sqrt_d));
    d_info->p[1] = static_cast<std::uint32_t>(uint_max * (1.0 - 1.0 / sqrt_d)) + d_info->p[0];

    generator_info = reinterpret_cast<void*>(d_info);
    generate_random_number = &generate_discrete;
  }


  
  unsigned int seed = seed_unset ? std::chrono::system_clock::now().time_since_epoch().count() : options.seed;
  unsigned int row_seed;
  
  boost::random::mt19937 seed_generator(seed);
  boost::math::normal_distribution<double> normal_dist(0.0, 1.0);
  // boost::random::normal_distribution<double> normal_generator;
  
  unsigned int i0 = 0;
  
  if (num_input_vectors == 0) {
    unsigned int i1 = 0;
    unsigned int n1 = static_cast<unsigned int>(LENGTH(indices_expr));
    int* u1 = indices;
    double* v1 = vs[0];

    while (i1 < n1) {
      if (u1[i1] > i0)
        seed_generator.discard(u1[i1] - i0);
      row_seed = seed_generator();
      boost::random::mt19937 row_generator(row_seed);
      
      if (options.row_start > 0)
        row_generator.discard(options.row_start);

      for (unsigned int j = 0; j < options.output_length; ++j) {
        double rn = generate_random_number(row_generator, generator_info);
        v1[j] += (rn - v1[j]) / static_cast<double>(i1 + 1);
      }
      i0 = u1[i1++] + 1;
    }

    for (unsigned int j = 0; j < options.output_length; ++j) {
      v1[j] *= static_cast<double>(i1) / sqrt(static_cast<double>(k));
    }
  } else {
    unsigned int* is = reinterpret_cast<unsigned int*>(alloca(num_vectors * sizeof(unsigned int)));
    unsigned int* ns = reinterpret_cast<unsigned int*>(alloca(num_vectors * sizeof(unsigned int)));
    
    for (unsigned int l = 0; l < num_vectors; ++l) {
      is[l] = 0;
      ns[l] = static_cast<unsigned int>(LENGTH(VECTOR_ELT(u_expr, l)));
    }
    
    bool finished = false;
    while (!finished) {
      int next_index = std::numeric_limits<int>::max();
      // finding the next iteration can be stored in a data structure to be made
      // faster, depending on the number of input vectors; it likely is fastest
      // to brute force for ~<= 5
      for (unsigned int l = 0; l < num_vectors; ++l) {
        if (is[l] < ns[l] && us[l][is[l]] < next_index)
          next_index = us[l][is[l]];
      }
      if (next_index > i0)
        seed_generator.discard(next_index - i0);
      row_seed = seed_generator();
      boost::random::mt19937 row_generator(row_seed);

      if (options.row_start > 0)
        row_generator.discard(options.row_start);

      for (unsigned int j = 0; j < options.output_length; ++j) {
        double rn = generate_random_number(row_generator, generator_info);
        for (unsigned int l = 0; l < num_vectors; ++l) {
          if (us[l][is[l]] == next_index) {
            vs[l][j] += (rn - vs[l][j]) / static_cast<double>(is[l] + 1);
          }
        }
      }
      i0 = next_index + 1;
      
      finished = true;
      for (unsigned int l = 0; l < num_vectors; ++l) {
        if (us[l][is[l]] == next_index)
          is[l]++;
        if (is[l] < ns[l])
          finished = false;
      }
    }
    
    for (unsigned int j = 0; j < options.output_length; ++j) {
      for (unsigned int l = 0; l < num_vectors; ++l) {
        vs[l][j] *= static_cast<double>(is[l]) / sqrt(static_cast<double>(k));
      }
    }
  }

  if (options.random_type == RANDOM_TYPE_NORMAL) {
    delete reinterpret_cast<NormalGeneratorInfo*>(generator_info);
  } else {
    delete reinterpret_cast<DiscreteGeneratorInfo*>(generator_info);
  }

  UNPROTECT(1);
  return result_expr;
}

} // extern "C"

#if __cplusplus >= 202002L
#  include <bit>
#else
#  include <cstring>

namespace std {

#  if __cplusplus >= 201103L
#    include <type_traits>

template <class To, class From>
typename std::enable_if<
  sizeof(To) == sizeof(From) &&
  std::is_trivially_copyable<From>::value &&
  std::is_trivially_copyable<To>::value,
  To>::type
// constexpr support needs compiler magic
bit_cast(const From& src) noexcept
{
  static_assert(std::is_trivially_constructible<To>::value,
    "This implementation additionally requires destination type to be trivially constructible");

  To dst;
  memcpy(&dst, &src, sizeof(To));
  return dst;
}

#  else

// We are only using this to cast function pointers, which are trivially copiable.
// is_trivially_copyable is compiler specific and isn't worth trying to drop in
// an implementation.
template <class To, class From>
To
bit_cast(const From& src)
{
  To dst;
  memcpy(&dst, &src, sizeof(To));
  return dst;
}

#  endif

}

#endif

extern "C" {

#define DEF_FUNC(_N_, _F_, _A_) { _N_, std::bit_cast<DL_FUNC>(&_F_), _A_ }

static R_CallMethodDef R_callMethods[] = {
  DEF_FUNC("esrp_generate_random_projections", generate_random_projections, 4),
  {NULL, NULL, 0}
};

#undef DEF_FUNC

void attribute_visible R_init_esrp(DllInfo *info) {
  R_registerRoutines(info, NULL, R_callMethods, NULL, NULL);
  R_useDynamicSymbols(info, static_cast<Rboolean>(FALSE));
}

} // extern "C"
