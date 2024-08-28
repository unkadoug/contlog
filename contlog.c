#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "contlog.h"
#define REP_NBITS (8*sizeof(CONTLOG_BASE))
#define SGNBIT_POS (REP_NBITS - 1)
#define MAGN_BITS (REP_NBITS + 1 - CONTLOG_SIGNED - CONTLOG_UNBOUNDED)
#define MINVAL (CONTLOG_SIGNED ? ((CONTLOG_BASE)1 << SGNBIT_POS) : 0)

#define ffs(X) _Generic((X),			\
			char: ffs,		\
			unsigned char: ffs,	\
			short: ffs,		\
			unsigned short: ffs,	\
			int: ffs,		\
			unsigned: ffs,		\
			long: ffsl,		\
			unsigned long: ffsl,	\
			long long: ffsll,	\
			unsigned long long: ffsll\
	  )(X)

static int
short_fls(contlog_t x) {
     int f = fls(x);
     return (f <= REP_NBITS? f : REP_NBITS);
}

#define fls(X) _Generic((X),			\
			char: short_fls,	\
			unsigned char: fls,	\
			short: short_fls,	\
			unsigned short: fls,	\
			int: fls,		\
			unsigned: fls,		\
			long: flsl,		\
			unsigned long: flsl,	\
			long long: flsll,	\
			unsigned long long: flsll\
	  )(X)

/*
 * For n >= d, compute the floor(log2(n/d)) with bit operations.
 */
static int
lgratio(fracpart_t n, fracpart_t d)
{
     if (n < d)
	  return (0);
     int lg = fls(n) - fls(d);
     n -= d << lg;
     n >>= 8 * sizeof(n) - 1;
     return (lg + n);
}

/*
 * Swap a pair of consecutive values, negating one of them.
 */
static void
contlog_swap_negate(int max_shift, int lo, contlog_t arg, fracpart_t pair[])
{
     pair[lo] += pair[lo^1];
     pair[lo^1] -= pair[lo];
     pair[lo] += pair[lo^1];
     if (arg != 0) {
	  int shift = ffs(arg) - 1;
	  pair[lo] -= pair[lo^1] >> shift;
	  if (MAGN_BITS - max_shift != shift + 1)
	       pair[lo^1] /= 2;
     }
}

static int
min(int a, int b)
{
     return (a < b ? a : b);
}

/* Complete the contlog binary encoding of the fraction in pair[]. */
static contlog_t
contlog_encode_exact(int max_shift, int lo, contlog_t arg, fracpart_t pair[])
{
     if (-pair[lo] >= 0) {
	  /* result <= 0; flip to positive value */
	  contlog_swap_negate(max_shift, lo, arg, pair);
	  lo ^= 1;
     }
     if (pair[lo] > pair[lo^1]) {
	  /* result > 1 */
#if (!CONTLOG_SIGNED && !CONTLOG_UNBOUNDED)
	  if (pair[lo] == pair[lo^1])
	       return (0);	/* avoids too-much bit shift */
#endif
	  if (arg == 0) {
		  if (MINVAL || pair[lo^1] != 0)
			  arg = 1;		/* result = 1 / result */
		  lo ^= 1;
	  } else {
		  /*
		   * 1 < result < 3, from encode_bounds.  Flip to a value in
		   * (1/4, 1/2), shift, and process one more bit.
		   */
		  pair[lo] -= pair[lo^1];
		  pair[lo^1] *= 2;
		  arg = 2 * arg + (lo? -1: 1); /* result = (result - 1) / 2 */
		  --max_shift;
	  }
     }

     for (;;) {
	  int shift = max_shift;
	  if (pair[lo] != 0)
	       shift = min(lgratio(pair[lo^1], pair[lo]), shift);
	  pair[lo] <<= shift;
	  arg <<= shift;	/* result <<= shift */
	  if ((max_shift -= shift) == 0)
	       break;
	  pair[lo^1] -= pair[lo];
	  arg = 2 * arg + (lo? -1: 1); /* result = (1 - result) / result */
	  lo ^= 1;
	  --max_shift;
     }
     return (arg);
}

/*
 * Given a numerator/denominator pair, compute its binary representation.
 */
contlog_t
contlog_encode_frac(fracpart_t pair[])
{

     if (-pair[1] > 0) {
	  pair[0] = -pair[0];
	  pair[1] = -pair[1];
     }

     return (contlog_encode_exact(MAGN_BITS, 0, 0, pair));
}

/*
 * Represent the partially-constructed contlog result of a computation.
 */
struct contlog_encode_state {
     contlog_t arg;		/* bits so far */
     int max_shift;		/* # left shifts of arg remaining */
     int lo;			/* position of numerator of lower bound */
};

static void
contlog_encode_state_init(struct contlog_encode_state *ces, fracpart_t quad[])
{
     fracpart_t mask = 0;
     for (int i = 0; i < 4; ++i)
	  mask |= quad[i] ^ (quad[i] >> 1);
     int shift = REP_NBITS- 1 - fls(mask);
     for (int i = 0; i < 4; ++i)
	  quad[i] <<= shift;
     ces->max_shift = MAGN_BITS;
     ces->lo = 0;
     ces->arg = 0;
     if (quad[0] >= quad[1]) {
	  /* result >= 1 */
	  if (MINVAL || quad[1] != 0)
	       ces->arg = 1;		/* result = 1 / result */
	  ces->lo = 3;
     }
}

/*
 * Extract as many bits as possible from quad and pack them into ces.  Return
 * true when enough bits are packed.  Return false when the lower and upper
 * bounds in quad are too far apart to extract more bits.
 */
static int
contlog_encode_bounds(struct contlog_encode_state *ces, fracpart_t quad[])
{
     int max_shift = ces->max_shift;
     int lo = ces->lo;
     contlog_t arg = ces->arg;
     if (-quad[lo^2] >= 0) {
	  /* result <= 0, flip to positive value */
	  contlog_swap_negate(max_shift, lo&1, arg, &quad[0]);
	  contlog_swap_negate(max_shift, lo&1, arg, &quad[2]);
	  lo ^= 3;
     }

     /* Extract bits into arg until either arg is filled, or lower
      * (quad[lo]/quad[lo^1]) and upper (quad[lo^2]/quad[lo^3]) bound ratios
      * have been pushed too far apart.
      */
     for (;;) {
	  /* Shift 'result' up without pushing its bounds outside [-1, 1] */
	  int shift = max_shift;
	  if (-quad[lo] > 0)
	       shift = min(lgratio(quad[lo^1], -quad[lo]), shift);
	  if (quad[lo^2] > 0)
	       shift = min(lgratio(quad[lo^3], quad[lo^2]), shift);
	  quad[lo] <<= shift;
	  quad[lo^2] <<= shift;
	  arg <<= shift;	/* result <<= shift */

	  if ((max_shift -= shift) == 0 || /* finished */
	      quad[lo] <= quad[lo^1] / 4)  /* result could be <= 1/4. */
	       break;

	  /*
	   * Result is > 1/4, so writing a '1' in arg is okay, because at worst
	   * it will have to be flipped to a '01'.
	   */
	  quad[lo^1] -= quad[lo];
	  quad[lo^3] -= quad[lo^2];
	  arg = 2 * arg + ((lo&1)? -1: 1); /* result = (1 - result) / result */
	  lo ^= 3;
	  --max_shift;
     }
     ces->lo = lo;
     ces->max_shift = max_shift;
     ces->arg = arg;
     return (max_shift == 0);
}

/*
 * Compute the overflow (0 or 1) that resulted when two values were added to
 * produce a sum.
 */
static int
add_overflow(fracpart_t add1, fracpart_t add2, fracpart_t sum)
{
     return (1 & (((add1 & add2) |
		   (add1 & ~sum) |
		   (add2 & ~sum)) >> SGNBIT_POS));
}

/*
 * Compute (x1 >> overflow) + x2, capturing numerical overflow.
 */
static void
oversum(fracpart_t sum[], int overflow,
	fracpart_t x1, fracpart_t x2)
{
     x1 >>= overflow;
     sum[1] = x1 + x2;
     sum[0] = (x1 >> SGNBIT_POS) + (x2 >> SGNBIT_POS) +
	  add_overflow(x1, x2, sum[1]);
}

/*
 * Find how much the largest of the values stored in pairs in sum[] overflowed
 * into the higher-order member of the pair.  Shift all the values right by that
 * amount and copy them back into result.  Return the overflow amount.
 */
static int
contlog_pack(fracpart_t result[], fracpart_t sum[])
{
     fracpart_t mask = 0;
     for (int i = 0; i < 2; ++i) {
	  fracpart_t hi = sum[2*i], lo = sum[2*i+1];
	  if (hi < 0) {
	       lo = -lo;
	       hi = ~hi + (lo == 0);
	  }
	  mask |= 2 * hi + (lo < 0);
     }
     int overflow = fls(mask);
     mask = ((fracpart_t)1 << 1 << (SGNBIT_POS - overflow)) - 1;
     for (int i = 0; i < 2; ++i)
	  result[i] = (sum[2*i] << 1 << (SGNBIT_POS - overflow)) |
	       ((sum[2*i+1] >> overflow) & mask);
     return (overflow);
}

/*
 * Given a quad (a, b, c, d) and an operand x, compute (a+cx)/(b+dx).
 */
static contlog_t
contlog_arith(contlog_t operand, fracpart_t quad[])
{
     int j = 2;
#if (CONTLOG_SIGNED)
     operand <<= 1;
#endif
#if (CONTLOG_UNBOUNDED)
     if (operand >> SGNBIT_POS) {
	  j ^= 2;
	  operand = -operand;
     }
     operand <<= 1;
#endif

     struct contlog_encode_state ces;
     contlog_encode_state_init(&ces, quad);
     int overflow = 0;
     while (operand != 0 || j == 0) {
	  /* Find the leftmost set bit position of operand and shift the bit
	   * out. */
	  if (operand != 0) {
	       int shift = REP_NBITS - fls(operand);
	       operand <<= shift;
	       operand = -2 * operand;
	       if (operand == 0 && j == 0)
		    ++shift;
	       overflow += shift;
	  }

	  /* Update quad to shrink range containing the result */
	  fracpart_t sum[4];
	  oversum(&sum[0], overflow, quad[j^0], quad[j^2]);
	  oversum(&sum[2], overflow, quad[j^1], quad[j^3]);
	  overflow = contlog_pack(&quad[j], sum);
	  if (contlog_encode_bounds(&ces, quad))
	       return (ces.arg);
	  j ^= 2;
     }
     return (contlog_encode_exact(ces.max_shift, ces.lo&1, ces.arg, quad));
}

/*
 * Translate operand into fraction frac[] = {numer, denom} that lies at the
 * 'center' of the interval of values represented by operand.
 */
static int
contlog_decode(contlog_t operand, ufracpart_t frac[])
{
     int neg = 0;
#if (CONTLOG_SIGNED)
     if (operand >> SGNBIT_POS) {
	  neg = 1;
	  operand = -operand;
     }
     operand <<= 1;
#endif
     int improper = 0;
     int zero = operand == 0;
#if (CONTLOG_UNBOUNDED)
     if (operand >> SGNBIT_POS) {
	  improper = 1;
	  operand = -operand;
     } else if (operand == 0)
	  improper = neg;
     operand <<= 1;
#endif
     improper ^= zero;
     fracpart_t pair[] = {0, 1};
     for (int lo = 1, lobit = operand ? ffs(operand) - 1: REP_NBITS;
	  !zero && lobit < REP_NBITS + 1; ) {
	  operand += lo ? (operand | -operand): (operand & -operand);
	  int nextbit = operand ? ffs(operand) - 1: REP_NBITS + lo;
	  lo ^= 1;
	  pair[lo] += pair[lo^1];
	  pair[lo] <<= nextbit - lobit - 1;
	  lobit = nextbit;
     }

     fracpart_t mask = pair[0] | pair[1];
     int shift = ffs(mask) - 1;
     frac[improper^1] = pair[0] >> shift;
     frac[improper^0] = pair[1] >> shift;
     return (neg);
}

/* Compute x+y */
contlog_t
contlog_add(contlog_t op0, contlog_t op1)
{
     int neg = op0 < 0;
     if (neg) {
	     op0 = -op0;
	     op1 = -op1;
     }
     ufracpart_t frac[2];
     if (contlog_decode(op1, frac))
	  frac[0] = -frac[0];
     fracpart_t quad[] = {frac[0], frac[1], frac[1], 0};
     contlog_t val = contlog_arith(op0, quad);
     return (neg ? -val : val);
}

/* Compute x-y */
contlog_t
contlog_sub(contlog_t op0, contlog_t op1)
{
     int neg = op0 < 0;
     if (neg) {
	     op0 = -op0;
	     op1 = -op1;
     }
     ufracpart_t frac[2];
     if (!contlog_decode(op1, frac))
	  frac[0] = -frac[0];
     fracpart_t quad[] = {frac[0], frac[1], frac[1], 0};
     contlog_t val = contlog_arith(op0, quad);
     return (neg ? -val : val);
}

/* Compute x*y */
contlog_t
contlog_mult(contlog_t op0, contlog_t op1)
{
     if (op0 < 0) {
	  op0 = -op0;
	  op1 = -op1;
     }
     ufracpart_t frac[2];
     int neg = contlog_decode(op1, frac);
     fracpart_t quad[] = {0, frac[1], frac[0], 0};
     contlog_t val = contlog_arith(op0, quad);
     return (neg ? -val : val);
}

/* Compute x*(1-y) */
contlog_t
contlog_compmult(contlog_t op0, contlog_t op1)
{
     if (op0 < 0) {
	  op0 = -op0;
	  op1 = -op1;
     }
     ufracpart_t frac[2];
     int neg = contlog_decode(op1, frac);
     if (neg)
	  frac[0] = -frac[0];
     frac[0] = frac[1] - frac[0];
     neg = frac[0] < 0;
     if (neg)
	  frac[0] = -frac[0];
     fracpart_t quad[] = {0, frac[1], frac[0], 0};
     contlog_t val = contlog_arith(op0, quad);
     return (neg ? -val : val);
}


/* Compute x/y */
contlog_t
contlog_div(contlog_t op0, contlog_t op1)
{
     if (op0 < 0) {
	  op0 = -op0;
	  op1 = -op1;
     }
     ufracpart_t frac[2];
     int neg = contlog_decode(op1, frac);
     fracpart_t quad[] = {0, frac[0], frac[1], 0};
     contlog_t val = contlog_arith(op0, quad);
     return (neg ? -val : val);
}

/* Compute (x+y)/(1-x*y) */
contlog_t
contlog_atnsum(contlog_t op0, contlog_t op1)
{
     int neg = op0 < 0;
     if (neg) {
	  op0 = -op0;
	  op1 = -op1;
     }
     if (CONTLOG_UNBOUNDED && op1 >= MINVAL - op0) {
	  op0 = MINVAL - op0;
	  op1 = MINVAL - op1;
	  neg ^= 1;
     }
     ufracpart_t frac[2];
     if (contlog_decode(op1, frac))
	  frac[0] = -frac[0];
     fracpart_t quad[] = {frac[0], frac[1], frac[1], -frac[0]};
     contlog_t val = contlog_arith(op0, quad);
     return (neg ? -val : val);
}

/* Compute (x*y)/(x+y) */
contlog_t
contlog_parallel(contlog_t op0, contlog_t op1)
{
     int neg = op0 < 0;
     if (neg) {
	  op0 = -op0;
	  op1 = -op1;
     }
     if (op1 < 0 && -op1 < op0) {
	  op1 += op0;		/* op0: op0+op1 */
	  op0 -= op1;		/* -op1; op0+op1 */
	  op1 = -op0 - op1;	/* -op1; -op0 */
	  neg ^= 1;
     }
     ufracpart_t frac[2];
     if (contlog_decode(op1, frac))
	  frac[1] = -frac[1];
     fracpart_t quad[] = {0, frac[0], frac[0], frac[1]};
     contlog_t val = contlog_arith(op0, quad);
     return (neg ? -val : val);
}

/*
 * Compute the high-order part of the product of a and b, where 'a*b' computes
 * the low-order part.
 */
static fracpart_t
hiprod(fracpart_t a, fracpart_t b)
{
     const unsigned int halfbits = REP_NBITS / 2;
     const fracpart_t halfmask = ((fracpart_t)1 << halfbits) - 1;
     fracpart_t a_lo = a & halfmask;
     fracpart_t a_hi = a >> halfbits;
     fracpart_t b_lo = b & halfmask;
     fracpart_t b_hi = b >> halfbits;
     fracpart_t hilo = (((a_lo * b_lo) >> halfbits) & halfmask) + a_hi * b_lo;
     fracpart_t lohi = (hilo & halfmask) + a_lo * b_hi;
     return ((hilo >> halfbits) + (lohi >> halfbits) + a_hi * b_hi);
}

/*
 * Compute (x1*y1 >> overflow) + x2*y2, capturing numerical overflow.
 */
static void
dotprod(fracpart_t sum[], int overflow,
	fracpart_t x1, fracpart_t x2,
	fracpart_t y1, fracpart_t y2)
{
     fracpart_t prod1 = x1 * y1;
     fracpart_t prod2 = x2 * y2;
     sum[0] = hiprod(x1, y1);
     if (overflow > 0) {
	  int backflow = REP_NBITS - overflow;
	  prod1 >>= overflow;
	  prod1 &= ((fracpart_t)1 << backflow) - 1;
	  prod1 |= sum[0] << backflow;
	  sum[0] >>= overflow;
     }
     sum[1] = prod1 + prod2;
     sum[0] += hiprod(x2, y2) + add_overflow(prod1, prod2, sum[1]);
}

/*
 * Replace x = column j with ((a*x)>>overflow) + b*y,
 * Return new overflow.
 */
static int
contlog_axpby(int overflow, int j, int d, fracpart_t quad[], fracpart_t a, fracpart_t b)
{
     fracpart_t sum[4];
     dotprod(&sum[0], overflow, quad[j^0], quad[j^2], a, b);
     dotprod(&sum[2], overflow, quad[j^1], quad[j^3], a, b);
     return (contlog_pack(&quad[d], sum));
}

/* Compute sqrt(numer/denom), where numer <= denom. */
static contlog_t
contlog_sqrt_frac(ufracpart_t numer, ufracpart_t denom)
{
     /*
      * Scale the argument up by a power of 4, to the range [1/4, 1), and the
      * result down by a power of 2, to avoid the higher iteration count that
      * comes with larger values.
      */
     int shift = lgratio(denom, numer) / 2;
     int ffsd = min(2 * shift, ffs(denom) - 1);
     numer <<= 2 * shift - ffsd;
     denom >>= ffsd;
     fracpart_t quad[] = {2*numer, numer+denom, 1, 1};
     struct contlog_encode_state ces;
     contlog_encode_state_init(&ces, quad);
     ces.max_shift -= shift;
  
     int overflow = 0;
     int j = 2;
     do {
	  /* Update quad to shrink range containing the result */
	  overflow = contlog_axpby(overflow, j, j, quad, denom-numer, 2);
	  overflow += contlog_axpby(overflow, j^2, j, quad, 0, numer);
	  j ^= 2;
     } while (!contlog_encode_bounds(&ces, quad));
     return (ces.arg);
}

/* Compute f(x) = sqrt(x) */
contlog_t
contlog_sqrt(contlog_t operand)
{
     ufracpart_t frac[2];
     if (contlog_decode(operand, frac))
	  return (MINVAL);
     int improper = CONTLOG_UNBOUNDED && frac[0] >= frac[1];
     contlog_t arg = contlog_sqrt_frac(frac[improper], frac[!improper]);
     if (improper)
	  arg = MINVAL - arg;
     return (arg);
}

/* Compute f(x) = x / sqrt(1 + x*x); ie, sin(arctan(x)) when 0 <= x < 1 */
static contlog_t
contlog_sinarctan_unit(contlog_t operand)
{
     ufracpart_t frac[2];
     (void)contlog_decode(operand, frac);
     fracpart_t numer = frac[0];
     fracpart_t denom = frac[1];
     fracpart_t quad[] = {0, 1, numer, denom};
     struct contlog_encode_state ces;
     contlog_encode_state_init(&ces, quad);
     int overflow = 0;
     int j = 0;
     do {
	  /* Update quad to shrink range containing the result */
	  overflow = contlog_axpby(overflow, j, j, quad, numer, 0);
	  overflow += contlog_axpby(overflow, j^2, j, quad, 2*denom, numer);
	  j ^= 2;
     } while (!contlog_encode_bounds(&ces, quad));
     return (ces.arg);
}

/* Compute f(x) = x / sqrt(1 + x*x); ie, sin(arctan(x)) */
contlog_t
contlog_sinarctan(contlog_t operand)
{
     int neg = operand < 0;
     if (neg) {
	  operand = -operand;
	  if (operand < 0)
	       return (-(operand >> 1));
     }
     contlog_t val;
     if (CONTLOG_UNBOUNDED && operand > MINVAL - operand) {
	  val = contlog_sinarctan_unit(MINVAL-operand);
	  val = contlog_mult(operand, val);
     } else
	  val = contlog_sinarctan_unit(operand);
     return (neg ? -val: val);
}

/* Compute sqrt(x*x + y*y) */
contlog_t
contlog_hypot(contlog_t op0, contlog_t op1)
{
     if (op0 < 0)
	  op0 = -op0;
     if (op1 < 0)
	  op1 = -op1;
     contlog_t div;
     if (op0 == op1)
	  div = contlog_sqrt_frac(1, 2);
     else {
	  if (op0 > op1) {
	       op1 = op0 - op1;
	       op0 -= op1;
	       op1 += op0;
	  }
	  if (op0 == 0)
	       return (op1);
	  div = contlog_sinarctan_unit(contlog_div(op0, op1));
     }
     return (contlog_div(op0, div));
}

/* Compute log(1 + numer/denom) */
static contlog_t
contlog_log1p_frac(fracpart_t numer, fracpart_t denom)
{
     fracpart_t int_numer = numer;
     fracpart_t odd_denom = 3 * denom;
     fracpart_t quad[] = {0, 1, numer, denom};
     int overflow = 0;
     int j = 0;
     struct contlog_encode_state ces;
     contlog_encode_state_init(&ces, quad);
     do {
	  /* Update quad to shrink range containing the result */
	  fracpart_t s1 = int_numer, s2 = 2;
	  if (j != 0) {
	       int_numer += numer;
	       s2 = odd_denom;
	       odd_denom += 2 * denom;
	  }
	  overflow = contlog_axpby(overflow, j, j, quad, s1, s2);
	  j ^= 2;
     } while (!contlog_encode_bounds(&ces, quad));
     return (ces.arg);
}

/* Compute log(1 + x) */
contlog_t
contlog_log1p(contlog_t operand)
{
     ufracpart_t frac[2];
     int neg = contlog_decode(operand, frac);
     ufracpart_t numer = frac[0];
     ufracpart_t denom = frac[1];
     /* log1p(-numer/denom) == -log1p(numer/(denom-numer)) */
     if (neg) {
	  if (denom <= numer)
	       return (MINVAL);
	  denom -= numer;
     }

     int shift = 0;
     if (numer >= denom) {
	  /*
	   * The log1p continued fraction converges slowly for large argments,
	   * which leads to overflows and inaccuracy.  To speed up for large
	   * arguments, use the identity: log1p(n/d) == shift * log1p(1) +
	   * log1p((n+d-D) / D), where D = d<<shift.  First, replace numer,denom
	   * with shift-adjusted values.
	   */
	  numer += denom;
	  shift = lgratio(numer, denom);
	  denom <<= shift;
	  int nzbits = ffs(numer | denom) - 1;
	  numer >>= nzbits;
	  denom >>= nzbits;
	  numer -= denom;
     }
     contlog_t arg = contlog_log1p_frac(numer, denom);
     if (shift != 0) {
	  /* Compute actual log1p from shift-adjusted log1p. */
	  contlog_t log2 = contlog_log1p_frac(1, 1);
	  (void)contlog_decode(log2, frac);
	  fracpart_t quad[] = {frac[0] * shift, frac[1], frac[1], 0};
	  arg = contlog_arith(arg, quad);
     }
     return (neg ? -arg: arg);
}

/* Compute 1/e**x. */
contlog_t
contlog_expm(contlog_t operand)
{
     ufracpart_t frac[2];
     int neg = contlog_decode(operand, frac);
     fracpart_t numer = frac[0];
     fracpart_t denom = frac[1];
     if (numer / (SGNBIT_POS + 1) >= denom)
	  return (neg? 0 : MINVAL);
     fracpart_t odd_denom = denom;
     fracpart_t quad[] = {0, 1, 1, 1};
     int overflow = 0;
     int j = 0;
     struct contlog_encode_state ces;
     contlog_encode_state_init(&ces, quad);
     do {
	  /* Update quad to shrink range containing the result */
	  fracpart_t s1 = numer, s2 = odd_denom;
	  if (j != 0) {
	       s1 = -s1;
	       s2 = 2;
	       odd_denom += 2 * denom;
	       ces.lo ^= 2;
	  }
	  overflow = contlog_axpby(overflow, j, j, quad, s1, s2);
	  j ^= 2;
     } while (j != 0 || !contlog_encode_bounds(&ces, quad));
     if (neg)
	  ces.arg = MINVAL - ces.arg;
     return (ces.arg);
}

/* Compute cos(sqrt(x)) or sin(sqrt(x))/sqrt(x). */
contlog_t
contlog_cssqrt(contlog_t operand, int n)
{
     ufracpart_t frac[2];
     int neg = contlog_decode(operand, frac);
     if (neg)
	  return (MINVAL);
     fracpart_t numer = frac[0];
     fracpart_t denom = frac[1];
     fracpart_t nn1_d = denom;
     fracpart_t quad[] = {0, 1, denom, denom};
     struct contlog_encode_state ces;
     contlog_encode_state_init(&ces, quad);
     int overflow = 0;
     int j = 0;
     do {
	  /* Update quad to shrink range containing the result */
	  overflow = contlog_axpby(overflow, j, j, quad, nn1_d, 0);
	  nn1_d = n * (n+1) * denom;
	  n += 2;
	  overflow += contlog_axpby(overflow, j^2, j, quad, nn1_d-numer, numer);
	  j ^= 2;
     } while (!contlog_encode_bounds(&ces, quad));
     return (ces.arg);
}

/* Compute cos(sqrt(x)). */
contlog_t
contlog_cosqrt(contlog_t operand)
{
     return (contlog_cssqrt(operand, 1));
}

/* Compute sin(sqrt(x))/sqrt(x). */
contlog_t
contlog_sisqrt(contlog_t operand)
{
     return (contlog_cssqrt(operand, 2));
}

/*
 * Given a pair of ratios that bound a range (open or closed), return the ratio
 * with smallest elements in that range.  Lower bound is bound[0]/bound[1];
 * upper bound is bound[2]/bound[3].  'frac' is a 4-element array; the result is
 * stored in its first two elements.
 */
static void
contlog_find_simplest(ufracpart_t bound[], int open, ufracpart_t frac[])
{
     int lo;
     fracpart_t gap = 0, val;
     for (lo = 0; gap == 0 && (bound[lo] != 0 || open); lo ^= 3) {
	  val = gap = 1;
	  if (bound[lo^2] != 0) {
	       val = bound[lo^3] / bound[lo^2];
	       bound[lo^3] %= bound[lo^2];
	  }
	  if (bound[lo] != 0) {
	       gap = bound[lo^1] / bound[lo] - val;
	       bound[lo^1] %= bound[lo];
	       if (bound[lo^1] == 0 && open)
		    --gap;	/* exclude range boundary */
	  }
	  if (gap != 0)
	       ++val;
	  frac[lo] += val * frac[lo^2];
	  frac[lo^1] += val * frac[lo^3];
     }
     lo ^= 3;
     lo &= 2;
     frac[0] = frac[lo];
     frac[1] = frac[lo^1];
}

/*
 * Translate operand into fraction frac[] = {numer, denom} that lies within the
 * interval of values represented by operand and has least numer+denom, for
 * presentation, not for calculation.  For even operand, frac may lie on the
 * boundary of the interval.
 */
int
contlog_decode_frac(contlog_t operand, ufracpart_t pair[])
{
     int nbits = REP_NBITS;
     int j = 0;
     int neg = 0;
#if (CONTLOG_SIGNED)
     if (operand >> SGNBIT_POS) {
	  neg = 1;
	  operand = -operand;
     }
     operand <<= 1;
     --nbits;
#endif
     int improper = 0;
     int zero = operand == 0;
     int big = 0;
#if (CONTLOG_UNBOUNDED)
     if (operand >> SGNBIT_POS) {
	  improper = 1;
	  operand = -operand;
     } else if (operand == 0)
	  improper = neg;
     operand <<= 1;
     --nbits;
#else
     zero &= !neg;
#if (!CONTLOG_SIGNED)
     /*
      * If the value is in [1/2, 1), negate operand to compute not numer/denom,
      * but (denom-numer)/(2*numer), because that has smaller intermediate
      * values and reduces overflow risk.  This eliminates overflow in the
      * all-ones case.  When completed, transfrom the result back to
      * numer/denom.
      */
     if (operand >> SGNBIT_POS) {
	  big = 1;
	  improper = 1;
	  operand = -operand;
     }
#endif
#endif

     /*
      * Find lower and upper bounds on the range of fractions represented by
      * 'operand'.
      */
     ufracpart_t quad[] = {0, 1, 1, 0};
     while (operand != 0) {
	  /* Find the leftmost set bit position of operand and shift the bit
	   * out. Negate the rest. */
	  int shift = REP_NBITS - fls(operand);
	  operand <<= shift;
	  operand = -2 * operand;
	  nbits -= shift + 1;
	  /* Update quad to shrink range containing the result */
	  j ^= 2;
	  quad[j^0] += (quad[j^2] <<= shift);
	  quad[j^1] += (quad[j^3] <<= shift);
     }

     /* Compute the simplest fraction in the interval, unless compution leads to
      * overflow; in that case assume that the exact value is simplest.
      */
     ufracpart_t frac[] = {1, 0, 0, 1};
     if (zero) {
	  frac[0] = 0;
	  frac[1] = 1;
     } else {
	  /*
	   * Start the computation of contlog_find_simplest before computing
	   * 'mid', the exact ratio, in order to reduce values so that
	   * subsequent calculation doesn't overflow.
	   */
	  for (int lo = 0; quad[lo^2] != 0; lo ^= 3) {
	       ufracpart_t val = quad[lo^3] / quad[lo^2];
	       quad[lo^1] -= val * quad[lo^0];
	       quad[lo^3] -= val * quad[lo^2];
	       frac[lo^0] += val * frac[lo^2];
	       frac[lo^1] += val * frac[lo^3];
	       if (quad[lo^0] == 0 || val != quad[lo^1] / quad[lo^0])
		    break;
	  }
	  /*
	   * Ignoring that everything was just transformed according to 'frac',
	   * compute 'mid', the exact ratio expressed by operand, and update
	   * 'quad' to the values of the lower and upper bounds of the set of
	   * values that map to operand.  Then find the simplest ratio in the
	   * interval.
	   */
	  ufracpart_t mid[] = {quad[j^0] + quad[j^2], quad[j^1] + quad[j^3]};
	  mid[0] <<= nbits;
	  mid[1] <<= nbits;
	  quad[j^0] += mid[0];
	  quad[j^1] += mid[1];
	  quad[j^2] += quad[j^2] + mid[0];
	  quad[j^3] += quad[j^3] + mid[1];
	  contlog_find_simplest(quad, nbits == 0, frac);
     }
     if (big) {
	  /* We have (denom-numer)/(2*numer).  Transform back to numer/denom. */
	  if (frac[1] % 2 == 0)
	       frac[0] += (frac[1] /= 2);
	  else
	       frac[0] += frac[0] + frac[1];
     }
     pair[0] = frac[improper];
     pair[1] = frac[!improper];
     return (neg);
}
