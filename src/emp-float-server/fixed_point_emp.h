#pragma once
#include <boost/assert.hpp>
#include <boost/concept_check.hpp>
#include <boost/operators.hpp>
#include <boost/static_assert.hpp>
#include <limits>
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#define __FPML_DEFINED_USE_MATH_DEFINES__
#endif
#include <errno.h>
#include <math.h>

/******************************************************************************/
/*                                                                            */
/* fixed_point_emp */
/*                                                                            */
/******************************************************************************/

template <
    /// The integer part bit count.
    unsigned char I,
    /// The fractional part bit count.
    unsigned char F>
/// A fixed point type.
//!
//! This type is designed to be a plug-in type to replace the floating point
//! types, such as float, double and long double. While it doesn't offer the
//! precision of these types, its operations are all implemented in integer
//! math, and it is therefore hoped that these operations are faster on non-
//! floating-point enabled hardware.
//!
//! The value uses 0/1 bits for the sign, I bits for the integer part and F bits
//! for the fractional part.
//!
//! Here is an example: a signed 8 bit 5:2 fixed_point_emp type would have the
//! following layout:
//!
//! fixed_point_emp<signed char, 5, 2>
//!
//!  sign           integer part \ / fractional part
//!   |                           |
//! +----+----+----+----+----+----+----+----+
//! | S  | I4 | I3 | I2 | I1 | I0 | F0 | F1 |
//! +----+----+----+----+----+----+----+----+
//!
//! where S is the sign-bit, I0 to I4 is the integer part, and F0 to F1 is
//! the fractional part. The range of this type is from -32 to +31.75, the
//! fractional part can encode multiples of 0.25.
//!
//! The class only implements a few operators directly, the others are generated
//! via the ordered_field_operators, unit_steppable and shiftable templates.
//!
//! The ordered_field_operators template takes and generates (==>) the
//! following operators:
//! +=   ==>  + (addable),
//! -=   ==>  - (subtractable),
//! *=   ==>  * (multipliable),
//! /=   ==>  / (dividable),
//! <    ==>  > , >=, <= (less_than_comparable),
//! ==   ==>  != (equality_comparable).
//!
//! The unit_steppable template takes and generates (==>) the following
//! operators:
//! ++   ==>  ++(int), preincrement versus postincrement, (incrementable),
//! --   ==>  --(int), predecrement versus postdecrement, (decrementable).
//!
//! The shiftable template takes and generates the following operators:
//! >>=  ==>  >> (right_shiftable),
//! <<=  ==>  << (left_shiftable).
class fixed_point_emp : public Swappable<fixed_point_emp<I, F>>,
                        public Comparable<fixed_point_emp<I, F>> {
  /// Grant the fixed_point_emp template access to private members. Types with
  /// different template parameters are different types and without this
  /// declaration they do not have access to private members.
  friend class fixed_point_emp;

  template <
      /// The power.
      int P,
      /// Make gcc happy.
      typename T = void>
  /// Calculate 2 to the power of F at compile time.
  //!
  //! The fixed_point_emp class needs 2 to the power of P in several locations
  //! in the code. However, the value depends on compile time constants only and
  //! can therefore be calculated at compile time using this template
  //! trickery. There is no need to call the function pow(2., P) at runtime to
  //! calculate this value.
  //!
  //! The value is calculated by recursively instantiating the power2 template
  //! with successively decrementing P. Finally, 2 to the power of 0 is
  //! terminating the recursion and set to 1.
  struct power2 {
    static const long long value = 2 * power2<P - 1, T>::value;
  };

  template <
      /// Make gcc happy.
      typename P>
  /// Calculate 2 to the power of 0 at compile time.
  //!
  //! The fixed_point_emp class needs 2 to the power of P in several locations
  //! in the code. However, the value depends on compile time constants only and
  //! can therefore be calculated at compile time using this template
  //! trickery. There is no need to call the function pow(2., P) at runtime to
  //! calculate this value.
  //!
  //! The value is calculated by recursively instantiating the power2 template
  //! with successively decrementing P. Finally, 2 to the power of 0 is
  //! terminating the recursion and set to 1.
  struct power2<0, P> {
    static const long long value = 1;
  };

  /// Initializing constructor.
  //!
  //! This constructor takes a value of type Integer and initializes the
  //! internal representation of fixed_point_emp<I, F> with it.
  fixed_point_emp(
      /// The internal representation to use for initialization.
      Integer value,
      /// This value is not important, it's just here to differentiate from
      /// the other constructors that convert its values.
      bool)
      : value_(value) {}

 public:
  /// The base type of this fixed_point_emp class.
  typedef Integer base_type;

  /// The integer part bit count.
  static const unsigned char integer_bit_count = I;

  /// The fractional part bit count.
  static const unsigned char fractional_bit_count = F;

  /// Default constructor.
  //!
  //! Just as with built-in types no initialization is done. The value is
  //! undetermined after executing this constructor.
  fixed_point_emp() {}

  template <
      /// The numeric type. Must be integer.
      typename T>
  /// Converting constructor.
  //!
  //! This constructor takes a numeric value of type T and converts it to
  //! this fixed_point_emp type.
  fixed_point_emp(
      /// The value to convert.
      T value,
      // Who owns the value
      int party = PUBLIC)
      : value_(I + F, value << F, party) {}

  /// Converting constructor.
  //!
  //! This constructor takes a numeric value of type bool and converts it to
  //! this fixed_point_emp type.
  fixed_point_emp(
      /// The value to convert.
      bool value,
      // Who owns the value
      int party = PUBLIC)
      : value_(I + F, value * power2<F>::value, party) {}

  /// Converting constructor.
  //!
  //! This constructor takes a numeric value of type float and converts it to
  //! this fixed_point_emp type.
  //!
  //! The conversion is done by multiplication with 2^F and rounding to the
  //! next integer.
  fixed_point_emp(
      /// The value to convert.
      float value,
      // Who owns the value
      int party = PUBLIC)
      : value_(I + F, value * power2<F>::value + (value >= 0 ? .5 : -.5),
               party) {}

  /// Converting constructor.
  //!
  //! This constructor takes a numeric value of type double and converts it to
  //! this fixed_point_emp type.
  fixed_point_emp(
      /// The value to convert.
      double value,
      // Who owns the value
      int party = PUBLIC)
      : value_(I + F, value * power2<F>::value + (value >= 0 ? .5 : -.5),
               party) {}

  /// Copy constructor.
  fixed_point_emp(fixed_point_emp<I, F> const& rhs) : value_(rhs.value_) {}

  /// Copy constructor.
  fixed_point_emp(Integer const& rhs) : value_(rhs) {}

  template <
      /// The other integer part bit count.
      unsigned char I2,
      /// The other fractional part bit count.
      unsigned char F2>
  /// Converting copy constructor.
  fixed_point_emp(
      /// The right hand side.
      fixed_point_emp<I2, F2> const& rhs)
      : value_(rhs.value_) {
    if (I - I2 > 0)
      value_ = value_ >> (I - I2);
    if (I2 - I > 0)
      value_ = value_ << (I2 - I);
  }

  /// Copy assignment operator.
  fixed_point_emp<I, F>& operator=(
      /// The right hand side.
      fixed_point_emp<I, F> const& rhs) {
    fixed_point_emp<I, F> temp(rhs);
    swap(temp);
    return *this;
  }

  template <
      /// The other integer part bit count.
      unsigned char I2,
      /// The other fractional part bit count.
      unsigned char F2>
  /// Converting copy assignment operator.
  fixed_point_emp<I, F>& operator=(
      /// The right hand side.
      fixed_point_emp<I2, F2> const& rhs) {
    fixed_point_emp<I, F> temp(rhs);
    swap(temp);
    return *this;
  }

  /// Exchanges the elements of two fixed_point_emp objects.
  void swap(
      /// The right hand side.
      fixed_point_emp<I, F>& rhs) {
    std::swap(value_, rhs.value_);
  }

  // Comparable
  Bit geq(const fixed_point_emp<I, F>& rhs) const {
    return value_.geq(rhs.value_);
  }

  Bit equal(const fixed_point_emp<I, F>& rhs) const {
    return value_.equal(rhs.value_);
  }

  // Swappable
  fixed_point_emp<I, F> select(const Bit& sel,
                               const fixed_point_emp<I, F>& rhs) const {
    return value_.select(sel, rhs.value_);
  }

  fixed_point_emp<I, F> operator^(const fixed_point_emp<I, F>& rhs) const {
    return value_ ^ rhs.value_;
  }

  fixed_point_emp<I, F> operator^=(const fixed_point_emp<I, F>& rhs) {
    return value_ ^= rhs.value_;
  }

  /// Negation operator.
  //!
  //! /return true if equal to zero, false otherwise.
  Bit operator!() const { return value_ == Integer(value_.size(), 0, PUBLIC); }

  /// Unary minus operator.
  fixed_point_emp<I, F> operator-() const {
    fixed_point_emp<I, F> result;
    result.value_ = -value_;
    return result;
  }

  /// Increment.
  // fixed_point_emp<I, F> & operator ++()
  //{
  //     value_ = value_ + Integer(value_.size(), power2<F>::value, PUBLIC);
  //     return *this;
  // }

  ///// Decrement.
  // fixed_point_emp<I, F> & operator --()
  //{
  //     value_ = value_ - Integer(value_.size(), power2<F>::value, PUBLIC);
  //     return *this;
  // }

  /// Addition.
  fixed_point_emp<I, F>& operator+=(fixed_point_emp<I, F> const& summand) {
    assert(value_.size() == summand.value_.size());
    size_t size = I + F;

    value_.resize(size * 2, true);
    vector<Bit> summand_bits = summand.value_.bits;
    summand_bits.resize(size * 2, true);
    add_full(value_.bits.data(), nullptr, value_.bits.data(),
             summand_bits.data(), nullptr, size * 2);
    value_.resize(size);

    return *this;
  }
  fixed_point_emp<I, F> operator+(fixed_point_emp<I, F> const& summand) {
    fixed_point_emp<I, F> res(*this);
    res += summand;
    return res;
  }

  /// Subtraction.
  fixed_point_emp<I, F>& operator-=(fixed_point_emp<I, F> const& diminuend) {
    assert(value_.size() == diminuend.value_.size());
    size_t size = I + F;
    value_.resize(size * 2, true);
    vector<Bit> diminuend_bits = diminuend.value_.bits;
    diminuend_bits.resize(size * 2, true);
    sub_full(value_.bits.data(), nullptr, value_.bits.data(),
             diminuend_bits.data(), nullptr, size * 2);
    // Bit MSB = value_.bits[(size*2)-1];
    value_.resize(size);
    // value_.bits[size-1] = MSB;
    return *this;
  }
  fixed_point_emp<I, F> operator-(fixed_point_emp<I, F> const& diminuend) {
    fixed_point_emp<I, F> res(*this);
    res -= diminuend;
    return res;
  }

  /// Multiplication.
  fixed_point_emp<I, F>& operator*=(fixed_point_emp<I, F> const& factor) {
    int size = (I + F);

    // std::cout << value_.reveal<string>() << endl;
    // std::cout << "*" << endl;
    // std::cout << factor.value_.reveal<string>() << endl;
    //
    //         //// do mult by hand instead of resize() then using * operator
    //         //// to take advantage of operators only being size bits (not
    //         2*size). Bit * sum = new Bit[size*2]; Bit * temp = new Bit[size];
    //         for(int i = 0; i < size*2; ++i)
    //             sum[i]=false;
    //         for(int i=0;i<size;++i) {
    //             for (int k = 0; k < size; ++k)
    //                 temp[k] = value_.bits[k] & factor.value_.bits[i];
    //
    //             add_full(sum+i, nullptr, sum+i, temp, nullptr, size);
    //
    //             sum[size+i] = sum[size+i-1]; // sign extend sum by one
    //             position
    //
    // vector<Bit> vsum(sum, sum + (2*size));
    // Integer isum(vsum);
    // std::cout << isum.reveal<string>() << endl;
    //         }
    //         memcpy(value_.bits.data(), sum+F, sizeof(Bit)*size);
    //         delete[] sum;
    //         delete[] temp;
    //         return *this;

    Integer fc(factor.value_);
    value_.resize(size * 2, true);
    fc.resize(size * 2, true);
    // value_ = (value_ * fc) >> F;
    value_ = value_ * fc;
    // right shift with sign extend in place
    for (int i = F; i < size; ++i)
      value_.bits[i - F] = value_.bits[i];
    for (int i = size - F; i < size; ++i)
      value_.bits[i] = value_.bits[size - 1];
    value_.resize(size);
    return *this;
  }
  fixed_point_emp<I, F> operator*(fixed_point_emp<I, F> const& factor) {
    fixed_point_emp<I, F> res(*this);
    res *= factor;
    return res;
  }

  /// Division.
  fixed_point_emp<I, F>& operator/=(
      fixed_point_emp<I, F> const& divisor)  // divisor cannot be const
  {
    Integer dc(divisor.value_);  // non const copy
    value_.resize((I + F) * 2, true);
    dc.resize((I + F) * 2, true);
    // fixed_point_emp<I, F> res((value_ << F) / dc);
    value_ = (value_ << F) / dc;
    value_.resize(I + F, true);
    return *this;
  }
  fixed_point_emp<I, F> operator/(
      fixed_point_emp<I, F> const& divisor)  // divisor cannot be const
  {
    fixed_point_emp<I, F> res(*this);
    res /= divisor;
    return res;
  }

  /// Shift right.
  fixed_point_emp<I, F>& operator>>=(
      /// Count of positions to shift.
      size_t shift) {
    value_ = value_ >> shift;
    return *this;
  }
  fixed_point_emp<I, F>& operator>>(
      /// Count of positions to shift.
      size_t shift) {
    fixed_point_emp res(*this);
    res >>= shift;
    return res;
  }

  /// Shift left.
  fixed_point_emp<I, F>& operator<<=(
      /// Count of positions to shift.
      size_t shift) {
    value_ = value_ << shift;
    return *this;
  }
  fixed_point_emp<I, F>& operator<<(
      /// Count of positions to shift.
      size_t shift) {
    fixed_point_emp res(*this);
    res <<= shift;
    return res;
  }

  double reveal(int party = PUBLIC) const {
    if (I + F == 32) {
      int32_t ct = value_.reveal<int32_t>(party);
      return ((double)ct) / power2<F>::value;
    } else if (I + F == 64) {
      int64_t ct = value_.reveal<int64_t>(party);
      return ((double)ct) / power2<F>::value;
    }
  }

  /**************************************************************************/
  /*                                                                        */
  /* fabs                                                                   */
  /*                                                                        */
  /**************************************************************************/

  /// Calculates the absolute value.
  //!
  //! The fabs function computes the absolute value of its argument.
  //!
  //! /return The absolute value of the argument.
  friend fixed_point_emp<I, F> fabs(
      /// The argument to the function.
      fixed_point_emp<I, F> x) {
    return x.value_.abs();
  }

  /**************************************************************************/
  /*                                                                        */
  /* fmod                                                                   */
  /*                                                                        */
  /**************************************************************************/

  /// Calculates the remainder.
  //!
  //! The fmod function computes the fixed point remainder of x/y.
  //!
  //! /return The fixed point remainder of x/y.
  friend fixed_point_emp<I, F> fmod(
      /// The argument to the function.
      fixed_point_emp<I, F> x,
      /// The argument to the function.
      fixed_point_emp<I, F> y) {
    fixed_point_emp<I, F> result;
    result.value_ = x.value_ % y.value_;
    return result;
  }

  /**************************************************************************/
  /*                                                                        */
  /* modf                                                                   */
  /*                                                                        */
  /**************************************************************************/

  ///// Split in integer and fraction parts.
  ////!
  ////! The modf function breaks the argument into integer and fraction parts,
  ////! each of which has the same sign as the argument. It stores the integer
  ////! part in the object pointed to by ptr.
  ////!
  ////! /return The signed fractional part of x/y.
  // friend fixed_point_emp<I, F> modf(
  //     /// The argument to the function.
  //     fixed_point_emp<I, F> x,
  //     /// The pointer to the integer part.
  //     fixed_point_emp<I, F> * ptr)
  //{
  //     fixed_point_emp<I, F> integer;
  //     integer.value_ = x.value_ & ~(power2<F>::value-1);
  //     *ptr = x < fixed_point_emp<I, F>(0) ?
  //         integer + fixed_point_emp<I, F>(1) : integer;

  //    fixed_point_emp<I, F> fraction;
  //    fraction.value_ = x.value_ & (power2<F>::value-1);

  //    return x < fixed_point_emp<I, F>(0) ? -fraction : fraction;
  //}

  /**************************************************************************/
  /*                                                                        */
  /* exp                                                                    */
  /*                                                                        */
  /**************************************************************************/

  ///// Calculates the exponential.
  ////!
  ////! The function computes the exponential function of x. The algorithm uses
  ////! the identity e^(a+b) = e^a * e^b.
  ////!
  ////! /return The exponential of the argument.
  // friend fixed_point_emp<I, F> exp(
  //     /// The argument to the exp function.
  //     fixed_point_emp<I, F> x)
  //{
  //     // a[x] = exp( (1/2) ^ x ), x: [0 ... 31]
  //     fixed_point_emp<I, F> a[] = {
  //         1.64872127070012814684865078781,
  //         1.28402541668774148407342056806,
  //         1.13314845306682631682900722781,
  //         1.06449445891785942956339059464,
  //         1.03174340749910267093874781528,
  //         1.01574770858668574745853507208,
  //         1.00784309720644797769345355976,
  //         1.00391388933834757344360960390,
  //         1.00195503359100281204651889805,
  //         1.00097703949241653524284529261,
  //         1.00048840047869447312617362381,
  //         1.00024417042974785493700523392,
  //         1.00012207776338377107650351967,
  //         1.00006103701893304542177912060,
  //         1.00003051804379102429545128481,
  //         1.00001525890547841394814004262,
  //         1.00000762942363515447174318433,
  //         1.00000381470454159186605078771,
  //         1.00000190735045180306002872525,
  //         1.00000095367477115374544678825,
  //         1.00000047683727188998079165439,
  //         1.00000023841860752327418915867,
  //         1.00000011920929665620888994533,
  //         1.00000005960464655174749969329,
  //         1.00000002980232283178452676169,
  //         1.00000001490116130486995926397,
  //         1.00000000745058062467940380956,
  //         1.00000000372529030540080797502,
  //         1.00000000186264515096568050830,
  //         1.00000000093132257504915938475,
  //         1.00000000046566128741615947508 };

  //    fixed_point_emp<I, F> e(2.718281828459045);

  //    fixed_point_emp<I, F> y(1);
  //    for (int i=F-1; i>=0; --i)
  //    {
  //        if (!(x.value_ & 1<<i))
  //            y *= a[F-i-1];
  //    }

  //    int x_int = (int)(floor(x));
  //    if (x_int<0)
  //    {
  //        for (int i=1; i<=-x_int; ++i)
  //            y /= e;
  //    }
  //    else
  //    {
  //        for (int i=1; i<=x_int; ++i)
  //            y *= e;
  //    }

  //    return y;
  //}

  /**************************************************************************/
  /*                                                                        */
  /* cos                                                                    */
  /*                                                                        */
  /**************************************************************************/

  /// Calculates the cosine.
  //!
  //! The algorithm uses a MacLaurin series expansion.
  //!
  //! First the argument is reduced to be within the range -Pi .. +Pi. Then
  //! the MacLaurin series is expanded. The argument reduction is problematic
  //! since Pi cannot be represented exactly. The more rounds are reduced the
  //! less significant is the argument (every reduction round makes a slight
  //! error), to the extent that the reduced argument and consequently the
  //! result are meaningless.
  //!
  //! The argument reduction uses one division. The series expansion uses 3
  //! additions and 4 multiplications.
  //!
  //! /return The cosine of the argument.
  friend fixed_point_emp<I, F> cos(
      /// The argument to the cos function.
      fixed_point_emp<I, F> x) {
    fixed_point_emp<I, F> x_ = fmod(x, fixed_point_emp<I, F>(M_PI * 2));

    // if (x_ > fixed_point_emp<I, F>(M_PI))
    //     x_ -= fixed_point_emp<I, F>(M_PI * 2);
    x_ = x_.If(x_ > fixed_point_emp<I, F>(M_PI),
               x_ - fixed_point_emp<I, F>(M_PI * 2));

    fixed_point_emp<I, F> xx = x_ * x_;

    fixed_point_emp<I, F> y =
        -xx * fixed_point_emp<I, F>(1. / (2 * 3 * 4 * 5 * 6));
    y += fixed_point_emp<I, F>(1. / (2 * 3 * 4));
    y *= xx;
    y -= fixed_point_emp<I, F>(1. / (2));
    y *= xx;
    y += fixed_point_emp<I, F>(1);

    return y;
  }

  /**************************************************************************/
  /*                                                                        */
  /* sin                                                                    */
  /*                                                                        */
  /**************************************************************************/

  /// Calculates the sine.
  //!
  //! The algorithm uses a MacLaurin series expansion.
  //!
  //! First the argument is reduced to be within the range -Pi .. +Pi. Then
  //! the MacLaurin series is expanded. The argument reduction is problematic
  //! since Pi cannot be represented exactly. The more rounds are reduced the
  //! less significant is the argument (every reduction round makes a slight
  //! error), to the extent that the reduced argument and consequently the
  //! result are meaningless.
  //!
  //! The argument reduction uses one division. The series expansion uses 3
  //! additions and 5 multiplications.
  //!
  //! /return The sine of the argument.
  friend fixed_point_emp<I, F> sin(
      /// The argument to the sin function.
      fixed_point_emp<I, F> x) {
    fixed_point_emp<I, F> x_ = fmod(x, fixed_point_emp<I, F>(M_PI * 2));
    // if (x_ > fixed_point_emp<I, F>(M_PI))
    //     x_ -= fixed_point_emp<I, F>(M_PI * 2);
    x_ = x_.If(x_ > fixed_point_emp<I, F>(M_PI),
               x_ - fixed_point_emp<I, F>(M_PI * 2));

    fixed_point_emp<I, F> xx = x_ * x_;

    fixed_point_emp<I, F> y =
        -xx * fixed_point_emp<I, F>(1. / (2 * 3 * 4 * 5 * 6 * 7));
    y += fixed_point_emp<I, F>(1. / (2 * 3 * 4 * 5));
    y *= xx;
    y -= fixed_point_emp<I, F>(1. / (2 * 3));
    y *= xx;
    y += fixed_point_emp<I, F>(1);
    y *= x_;

    return y;
  }

  /**************************************************************************/
  /*                                                                        */
  /* sqrt                                                                   */
  /*                                                                        */
  /**************************************************************************/

  /// Calculates the square root.
  //!
  //! The sqrt function computes the nonnegative square root of its argument.
  //!
  //! Calculates an approximation of the square root using an integer
  //! algorithm. The algorithm is described in Wikipedia:
  //! http://en.wikipedia.org/wiki/Methods_of_computing_square_roots
  //!
  //! The algorithm seems to have originated in a book on programming abaci by
  //! Mr C. Woo.
  //!
  //! /return The square root of the argument. If the argument is negative,
  //! the function returns 0.
  friend fixed_point_emp<I, F> sqrt(fixed_point_emp<I, F> x) {
    int size = (I + F);
    Integer op(x.value_);
    op.resize(size * 2, true);
    op = op << I - 1;

    Integer res(size * 2, 0, PUBLIC);

    int64_t o = 1;
    o <<= ((size * 2) - 1 - 1);  // additional -1 for sign bit
    Integer one(size * 2, o, PUBLIC);
    // std::cout << one.reveal<string>() << endl;

    // while (one > op)
    //     one >>= 2;
    for (int i = 0; i < size; i++)
      one = one.If(one > op, one >> 2);

    // std::cout << one.reveal<string>() << endl;
    // std::cout << op.reveal<string>() << endl;

    // while (one != 0)
    for (int i = 0; i < size * 2; i++)  // is 2* required?
    {
      // if (op >= res + one)
      //{
      //     op = op - (res + one);
      //     res = res + (one << 1);
      // }
      Bit ogtrpo = (op >= (res + one));
      op = op.If(ogtrpo, op - (res + one));
      res = res.If(ogtrpo, res + (one << 1));
      // res >>= 1;
      res = res >> 1;
      // one >>= 2;
      one = one >> 2;
      std::cout << res.reveal<string>() << endl;
    }
    std::cout << res.reveal<string>() << endl;

    fixed_point_emp<I, F> root(res);
    root.value_.resize(I + F, true);
    return root;
  }

  // private:
  /// The value in fixed point format.
  Integer value_;
};

template <
    /// The input stream type.
    typename S,
    /// The integer part bit count of the fixed_point_emp type.
    unsigned char I,
    /// The fractional part bit count of the fixed_point_emp type.
    unsigned char F>
/// Stream input operator.
//!
//! A value is first input to type double and then the read value is converted
//! to type fixed_point_emp before it is returned.
//!
//! /return A reference to this input stream.
S& operator>>(
    /// The input stream.
    S& s,
    /// A reference to the value to be read.
    fixed_point_emp<I, F>& v) {
  double value = 0.;
  s >> value;
  if (s)
    v = value;
  return s;
}

template <
    /// The output stream type.
    typename S,
    /// The integer part bit count of the fixed_point_emp type.
    unsigned char I,
    /// The fractional part bit count of the fixed_point_emp type.
    unsigned char F>
/// Stream output operator.
//!
//! The fixed_point_emp value is first converted to type double and then the
//! output operator for type double is called.
//!
//! /return A reference to this output stream.
S& operator<<(
    /// The output stream.
    S& s,
    /// A const reference to the value to be written.
    fixed_point_emp<I, F> const& v) {
  double value;
  value = v.reveal();
  s << value;
  return s;
}

#ifdef __FPML_DEFINED_USE_MATH_DEFINES__
#undef _USE_MATH_DEFINES
#undef __FPML_DEFINED_USE_MATH_DEFINES__
#endif
