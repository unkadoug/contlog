#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "contlog.h"
#define SIGNED(T) (-(T)1 < 0)
#define MAXBITS(T) (8*sizeof(T)-(SIGNED(T)?1:0))
#define SGNBIT_POS(T) (8*sizeof(T) - 1)
#define MINVAL(T) (SIGNED(T) ? (T)1 << SGNBIT_POS(T) : 0)

#define ffs(X) _Generic((X), \
		char: ffs, \
		unsigned char: ffs, \
		short: ffs, \
		unsigned short: ffs, \
		int: ffs, \
		unsigned: ffs, \
		long: ffsl,	\
		unsigned long: ffsl,	\
		long long: ffsll \
		)(X)
#define fls(X) _Generic((X), \
		char: fls, \
		unsigned char: fls, \
		short: fls, \
		unsigned short: fls, \
		int: fls, \
		unsigned: fls, \
		long: flsl, \
		unsigned long: flsl, \
		long long: flsll \
		)(X)

static inline int
lobit(contlog_t operand)
{
  int invpos = MAXBITS(contlog_t);
  return (operand ? ffs(operand) - 1 : invpos);
}

static void
contlog_op_to_frac(contlog_t operand, fracpart_t frac[], int lo)
{
  int w = lobit(operand);
  while (operand != 0) {
    lo ^= 1;
    frac[lo] += frac[lo^1];
    operand &= operand - 1;
    int w1 = lobit(operand);
    frac[lo] <<= w1 - w - 1;
    w = w1;
  }
}

/*
 * Translate operand into fraction frac[] = {numer, denom}.
 */
static int
contlog_decode(contlog_t operand, fracpart_t frac[])
{
  int neg = 0;
  if (SIGNED(contlog_t)) {
    contlog_t hibit = (contlog_t)1 << SGNBIT_POS(contlog_t);
    neg = (operand & hibit) != 0;
    operand ^= operand << 1;
    operand &= ~hibit;
  }
  else
    operand ^= operand << 1;
  frac[neg] = 0;
  frac[!neg] = 1;
  contlog_op_to_frac(operand, frac, !neg);
  int shift = ffs(frac[0] | frac[1]) - 1;
  frac[0] >>= shift;
  frac[1] >>= shift;
  return (neg);
}

static void
contlog_to_frac_ubound(contlog_t operand, fracpart_t frac[])
{
  operand ^= (operand << 1) | 1;
  if (SIGNED(contlog_t))
    operand &= ~((contlog_t)1 << SGNBIT_POS(contlog_t));
  frac[0] = operand & -operand;
  frac[1] = 1;
  contlog_op_to_frac(operand, frac, 0);
}

void
contlog_decode_frac(contlog_t operand, fracpart_t pair[])
{
  contlog_t hibit = (contlog_t)1 << SGNBIT_POS(contlog_t);
  int neg = SIGNED(contlog_t) && (operand & hibit) != 0;
  if (neg)
    operand = -operand;
  int improper = MINVAL(contlog_t) - operand <= operand;
  if (improper)
    operand = MINVAL(contlog_t) - operand;
  int lo = 0;
  fracpart_t frac[] = {1, 0, 0, 1};

  if (operand != 0) {
    contlog_t bound[4];
    contlog_to_frac_ubound(operand-1, (fracpart_t *)&bound[0]);
    contlog_to_frac_ubound(operand, (fracpart_t *)&bound[2]);
    
    for (;; lo ^= 3) {
      fracpart_t gap = 1;
      fracpart_t val = bound[lo^3] / bound[lo^2];
      if (bound[lo] != 0 && val == bound[lo^1] / bound[lo]) {
	bound[lo^1] %= bound[lo];
	bound[lo^3] %= bound[lo^2];
	gap = 0;
      }
      val += gap;
      frac[lo] += val * frac[lo^2];
      frac[lo^1] += val * frac[lo^3];
      if (gap != 0)
	break;
    }
  }
  lo &= 2;
  if (improper)
    lo ^= 1;
  pair[0] = frac[lo];
  pair[1] = frac[lo^1];
  if (neg)
    pair[0] = -pair[0];
}

static int
lgratio(fracpart_t n, fracpart_t d)
{
  if (d == 0)
    return (8 * sizeof(contlog_t));

  int lg = fls(n) - fls(d);
  if (lg >= 0) {
    if (n < d << lg)
      --lg;
  } else {
    abort();
    if (n << -lg < d)
      --lg;
  }
  return (lg);
}

static int
min(int a, int b)
{
  return (a < b ? a : b);
}

static contlog_t
contlog_encode_exact(int nbits, int lo, contlog_t arg, fracpart_t pair[])
{
  if (-pair[lo] > 0) {
    /* result < 0 */
    int shift = ffs(arg) - 1;
    pair[lo] += pair[lo^1];
    pair[lo^1] -= pair[lo];
    pair[lo] += pair[lo^1];
    if (shift >= 0) {
      pair[lo] -= pair[lo^1] >> shift;
      if (MAXBITS(contlog_t) != nbits + shift + 1)
	pair[lo^1] /= 2;
    }
    lo ^= 1;
  }
  if (pair[lo] >= pair[lo^1]) {
    /* result >= 1 */
    if (SIGNED(contlog_t) || pair[lo^1] != 0)
      arg = 2 * (arg - lo) + 1;
    lo ^= 1;
  }

  do {
    int shift = min(lgratio(pair[lo^1], pair[lo]), nbits);
    pair[lo] <<= shift;
    arg <<= shift;
    if ((nbits -= shift) == 0)
      break;
    pair[lo^1] -= pair[lo];
    arg = 2 * (arg - lo) + 1;
    lo ^= 1;
  } while (--nbits != 0);
  return (arg);
}

contlog_t
contlog_encode_frac(fracpart_t pair[])
{

  if (-pair[1] > 0) {
    pair[0] = -pair[0];
    pair[1] = -pair[1];
  }

  return (contlog_encode_exact(MAXBITS(contlog_t), 0, 0, pair));
}

struct contlog_encode_state {
  contlog_t arg;
  int nbits;
  int lo;
};

static void
contlog_encode_state_init(struct contlog_encode_state *ces, fracpart_t quad[])
{
  int nbits = 8 * sizeof(contlog_t);
  fracpart_t mask = 0;
  for (int i = 0; i < 4; ++i)
    mask |= quad[i] ^ (quad[i] >> 1);
  int shift = nbits - fls(mask) - 1;
  for (int i = 0; i < 4; ++i)
    quad[i] <<= shift;
  ces->nbits = MAXBITS(contlog_t);
  ces->lo = 0;
  ces->arg = 0;
}

static int
contlog_encode_bounds(struct contlog_encode_state *ces, fracpart_t quad[])
{
  int nbits = ces->nbits;
  int lo = ces->lo;
  contlog_t arg = ces->arg;
  if (-quad[lo^2] > 0) {
    /* result < 0 */
    int shift = ffs(arg) - 1;
    for (int i = lo&1; i < 4; i += 2) {
      quad[i] += quad[i^1];
      quad[i^1] -= quad[i];
      quad[i] += quad[i^1];
      if (shift >= 0) {
	quad[i] -= quad[i^1] >> shift;
	if (MAXBITS(contlog_t) != nbits + shift + 1)
	  quad[i^1] /= 2;
      }
    }
    lo ^= 3;
  }
  if (quad[lo] >= quad[lo^1]) {
    /* result >= 1 */
    if (SIGNED(contlog_t) || quad[lo^1] != 0)
      arg = 2 * (arg - (lo&1)) + 1;
    lo ^= 3;
  }

  /* Extract bits into arg until either arg is filled,
   * or lower and upper bound ratios have been pushed too far apart.
   */
  do {
    if (quad[lo^2] <= quad[lo^3] / 2) {
      int shift = nbits;
      if (-quad[lo] >= 0) {
        if (-quad[lo] <= quad[lo^1] / 2)
          shift = min(lgratio(quad[lo^1], -quad[lo]), shift);
	else
	  break;
      }
      /* -1/2 <= result <= 1/2 */
      shift = min(lgratio(quad[lo^3], quad[lo^2]), shift);
      quad[lo] <<= shift;
      quad[lo^2] <<= shift;
      arg <<= shift;
      if ((nbits -= shift) == 0)
        break;
    }
    if (quad[lo] <= 0 ||
	quad[lo] <= (quad[lo^1] - quad[lo]) / 2 ||
	quad[lo^3] <= quad[lo^2] / 2 ||
	(quad[lo^3] - quad[lo^2] <= quad[lo^2] / 2 &&
	 quad[lo] <= quad[lo^1] / 2))
      break;
    /* 1/2 < result < 2, or 1/3 < result < 2/3 */
    quad[lo^1] -= quad[lo];
    quad[lo^3] -= quad[lo^2];
    arg = 2 * (arg - (lo&1)) + 1;
    lo ^= 3;
  } while (--nbits != 0);
  ces->lo = lo;
  ces->nbits = nbits;
  ces->arg = arg;
  return (nbits == 0);
}

static int
pack(int n, fracpart_t result[], fracpart_t sum[])
{
  fracpart_t mask = 0;
  for (int i = 0; i < n; ++i) {
    fracpart_t s = sum[2*i];
    s ^= (s << 1) | ((sum[2*i+1] >> SGNBIT_POS(fracpart_t)) & 1);
    mask |= s;
  }
  int overflow = fls(mask);
  mask = ((fracpart_t)1 << 1 << (SGNBIT_POS(fracpart_t) - overflow)) - 1;
  for (int i = 0; i < n; ++i)
    result[i] = (sum[2*i] << 1 << (SGNBIT_POS(fracpart_t) - overflow)) |
	    ((sum[2*i+1] >> overflow) & mask);
  return (overflow);
}

static int
add_overflow(fracpart_t add1, fracpart_t add2, fracpart_t sum)
{
  return (1 & (((add1 & add2) |
		(add1 & ~sum) |
		(add2 & ~sum)) >> SGNBIT_POS(fracpart_t)));
}

static void
oversum(fracpart_t sum[], int overflow,
	fracpart_t x1, fracpart_t x2)
{
  int nbits = SGNBIT_POS(fracpart_t);
  x1 >>= overflow;
  sum[0] = (x1 >> nbits) + (x2 >> nbits) + add_overflow(x1, x2, x1 + x2);
  sum[1] = x1 + x2;
}

static void
debug_print(contlog_t operand, fracpart_t quad[], int j)
{
#ifdef DEBUG
  int sh = 4*sizeof(contlog_t);
  __uintmax_t val = operand;
  val -= val >> sh >> sh << sh << sh;
  fprintf(stderr, " %0*jx", (int)(2*sizeof(operand)), val);
  for (int i = 0; i < 4; ++i) {
    __intmax_t val = quad[i^1];
    val -= val >> sh >> sh << sh << sh;
    fprintf(stderr, " %10jx", val);
  }
  fprintf(stderr, "%12g%12g",
	  (double)quad[j^0]/quad[j^1],
	  (double)quad[j^2]/quad[j^3]);
  fprintf(stderr, "\n");
#endif
}

static contlog_t
contlog_arith(contlog_t operand, fracpart_t quad[])
{
  contlog_t hibit = (contlog_t)1 << SGNBIT_POS(contlog_t);
  if (SIGNED(contlog_t)) {
    if (operand & hibit) {
      for (int i = 0; i < 2; ++i) {
	quad[i] -= quad[i^2];
	quad[i^2] += quad[i];
	quad[i] -= quad[i^2];
      }
    }
    operand <<= 1;
  }
  int j = (operand & hibit)? 2 : 0;
  operand ^= operand << 1;

  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  int overflow = 0;
  while (operand != 0) {
    debug_print(operand, quad, j);
    j ^= 2;

    /* Find the leftmost set bit position of operand and shift the bit out. */
    int shift = SGNBIT_POS(contlog_t) - (fls(operand) - 1);
    operand <<= shift;
    operand <<= 1;
    overflow += shift;

    /* Update quad to shrink range containing the result */
    fracpart_t sum[4];
    oversum(&sum[0], overflow, quad[j^0], quad[j^2]);
    oversum(&sum[2], overflow, quad[j^1], quad[j^3]);
    overflow = pack(2, &quad[j], sum);
    debug_print(operand, quad, j);
    if (contlog_encode_bounds(&ces, quad))
      break;
  }
  if (operand == 0)
    debug_print(operand, quad, j);
    ces.arg = contlog_encode_exact(ces.nbits, ces.lo&1, ces.arg, &quad[j]);
  return (ces.arg);
}

static void
contlog_upshift(fracpart_t frac[])
{
  int shift = SGNBIT_POS(fracpart_t) - fls(frac[0] | frac[1]);
  frac[0] <<= shift;
  frac[1] <<= shift;
}

contlog_t
contlog_add(contlog_t op0, contlog_t op1)
{
  fracpart_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  if (neg)
    frac[0] = -frac[0];
  fracpart_t quad[] = {frac[0], frac[1], frac[1], 0};
  return (contlog_arith(op0, quad));
}

contlog_t
contlog_sub(contlog_t op0, contlog_t op1)
{
  fracpart_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  if (!neg)
    frac[0] = -frac[0];
  fracpart_t quad[] = {frac[0], frac[1], frac[1], 0};
  return (contlog_arith(op0, quad));
}

contlog_t
contlog_mult(contlog_t op0, contlog_t op1)
{
  fracpart_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  fracpart_t quad[] = {0, frac[1], frac[0], 0};
  contlog_t val = contlog_arith(op0, quad);
  return (neg ? -val : val);
}

contlog_t
contlog_div(contlog_t op0, contlog_t op1)
{
  fracpart_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  fracpart_t quad[] = {0, frac[0], frac[1], 0};
  contlog_t val = contlog_arith(op0, quad);
  return (neg ? -val : val);
}

contlog_t
contlog_atnsum(contlog_t op0, contlog_t op1)
{
  fracpart_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  if (neg)
    frac[0] = -frac[0];
  fracpart_t quad[] = {frac[0], frac[1], frac[1], -frac[0]};
  return (contlog_arith(op0, quad));
}

contlog_t
contlog_harmsum(contlog_t op0, contlog_t op1)
{
  fracpart_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  if (neg)
    frac[0] = -frac[0];
  fracpart_t quad[] = {0, frac[0], frac[0], frac[1]};
  return (contlog_arith(op0, quad));
}

contlog_t
contlog_hypot(contlog_t op0, contlog_t op1)
{
  if (op0 < 0)
    op0 = -op0;
  else if (op0 == 0)
    return (op1);
  if (op1 < 0)
    op1 = -op1;
  else if (op1 == 0)
    return (op0);
  if (op1 > op0) {
    contlog_t tmp = op0;
    op0 = op1;
    op1 = tmp;
  }
  return (contlog_div(op1, contlog_recip_hypot1(contlog_div(op0, op1))));
}

static fracpart_t
hiprod(fracpart_t a, fracpart_t b)
{
  const unsigned int halfbits = 4 * sizeof(fracpart_t);
  const fracpart_t halfmask = ((fracpart_t)1 << halfbits) - 1;
  fracpart_t a_lo = a & halfmask;
  fracpart_t a_hi = a >> halfbits;
  fracpart_t b_lo = b & halfmask;
  fracpart_t b_hi = b >> halfbits;
  fracpart_t hilo = (((a_lo * b_lo) >> halfbits) & halfmask) + a_hi * b_lo;
  fracpart_t lohi = (hilo & halfmask) + a_lo * b_hi;
  return ((hilo >> halfbits) + (lohi >> halfbits) + a_hi * b_hi);
}

static void
dotprod(fracpart_t sum[],
	fracpart_t x1, fracpart_t x2,
	fracpart_t y1, fracpart_t y2)
{
  fracpart_t prod1 = x1 * y1;
  fracpart_t prod2 = x2 * y2;
  sum[0] = hiprod(x1, y1) + hiprod(x2, y2) +
	  add_overflow(prod1, prod2, prod1 + prod2);
  sum[1] = prod1 + prod2;
}

static void
dotprod2(fracpart_t sum[], int overflow,
	fracpart_t x1, fracpart_t x2,
	fracpart_t y1, fracpart_t y2)
{
  fracpart_t prod1 = x1 * y1;
  fracpart_t prod2 = x2 * y2;
  sum[0] = hiprod(x1, y1);
  if (overflow > 0) {
    int backflow = 8 * sizeof(fracpart_t) - overflow;
    prod1 >>= overflow;
    prod1 &= ((fracpart_t)1 << backflow) - 1;
    prod1 |= sum[0] << backflow;
    sum[0] >>= overflow;
  }
  sum[0] += hiprod(x2, y2) + add_overflow(prod1, prod2, prod1 + prod2);
  sum[1] = prod1 + prod2;
}

/* Compute f(x) = sqrt(x) */
contlog_t
contlog_sqrt(contlog_t operand)
{
  fracpart_t frac[2];
  if (contlog_decode(operand, frac))
    return (MINVAL(contlog_t));
  int improper = frac[0] >= frac[1];
  fracpart_t numer = frac[improper];
  fracpart_t denom = frac[!improper];

  /*
   * Scale the argument down by a power of 4, to the range [1, 4), and the
   * result up by a power of 2, to avoid the higher iteration count that comes
   * with larger values.
   */
  int shift = lgratio(denom, numer) / 2;
  numer <<= 2 * shift;
  fracpart_t quad[] = {2*numer, numer+denom, 1, 1};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  ces.nbits -= shift;
  
  int overflow = 0;
  int j = 0;
  for (;;) {
    /* Update quad to shrink range containing the result */
    fracpart_t sum[4];
    j ^= 2;
    dotprod2(&sum[0], overflow,
	     quad[j^0], quad[j^2], denom-numer, 2);
    dotprod2(&sum[2], overflow,
	     quad[j^1], quad[j^3], denom-numer, 2);
    overflow = pack(2, &quad[j], sum);
    if (contlog_encode_bounds(&ces, quad))
      break;
    dotprod2(&sum[0], 0,
	     quad[j^0], quad[j^2], numer, 0);
    dotprod2(&sum[2], 0,
	     quad[j^1], quad[j^3], numer, 0);
    overflow += pack(2, &quad[j], sum);
  }
  if (improper)
    ces.arg = MINVAL(contlog_t) - ces.arg;
  return (ces.arg);
}

/* Compute f(x) = 1 / sqrt(1 + x*x) */
contlog_t
contlog_recip_hypot1(contlog_t operand)
{
  if (operand < 0) {
    operand = -operand;
    if (operand < 0)
      return (-(operand >> 1));
  }
  fracpart_t frac[2];
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
    fracpart_t sum[4];
    dotprod2(&sum[0], overflow,
	     quad[j^0], quad[j^2], numer*numer, 2*denom);
    dotprod2(&sum[2], overflow,
	     quad[j^1], quad[j^3], numer*numer, 2*denom);
    overflow = pack(2, &quad[j], sum);
    j ^= 2;
  } while (!contlog_encode_bounds(&ces, quad));
  return (ces.arg);
}

static contlog_t
contlog_log1p_frac(fracpart_t numer, fracpart_t denom)
{
  fracpart_t n_11234 = numer;
  fracpart_t n_12345 = numer;
  fracpart_t d_13579 = denom;
  fracpart_t quad[] = {0, 1, 1, 0};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  do {
    /* Update quad to shrink range containing the result */
    fracpart_t sum[8];
    for (int i = 0; i < 2; ++i) {
      int j = i ^ 2;
      dotprod(&sum[2*i], quad[i], quad[j], n_12345 + 2*d_13579, 2*n_11234);
      dotprod(&sum[2*j], quad[i], quad[j], d_13579, n_11234);
    }
    (void)pack(4, quad, sum);
    n_11234 = n_12345;
    n_12345 += numer;
    d_13579 += 2 * denom;
  } while (!contlog_encode_bounds(&ces, quad));
  return (ces.arg);
}

contlog_t
contlog_log1p(contlog_t operand)
{
  fracpart_t frac[2];
  int neg = contlog_decode(operand, frac);
  int numer = frac[0];
  int denom = frac[1];
  /* log1p(-numer/denom) == -log1p(numer/(denom-numer)) */
  if (neg) {
    if (denom <= numer)
      return (MINVAL(contlog_t));
    denom -= numer;
  }
  int shift = 0;
  if (numer/4 >= denom) {
    /*
     * The log1p continued fraction converges slowly for large argments, which
     * leads to overflows and inaccuracy.  To speed up for arguments >= 4, use
     * the identity:
     * log1p(n/d) == shift * log1p(1) + log1p((n+d-D) / D),
     * where D = d<<shift.
     * First, replace numer,denom with shift-adjusted values.
     */
    numer += denom;
    shift = lgratio(numer, denom);
    denom <<= shift;
    numer -= denom;
  }
  contlog_t arg = contlog_log1p_frac(numer, denom);
  if (shift != 0) {
    /* Compute actual log1p from shift-adjusted log1p. */
    contlog_t log2 = contlog_log1p_frac(1, 1);
    (void)contlog_decode(log2, frac);
    frac[0] *= shift;
    contlog_upshift(frac);
    fracpart_t quad[] = {frac[0], frac[1], frac[1], 0};
    arg = contlog_arith(arg, quad);
  }
  if (neg)
    arg = -arg;
  return (arg);
}

/* Compute e**x.
 */
contlog_t
contlog_exp(contlog_t operand)
{
  fracpart_t frac[2];
  int neg = contlog_decode(operand, frac);
  fracpart_t numer = frac[0];
  fracpart_t denom = frac[1];
  if (numer / (SGNBIT_POS(contlog_t) + 1) >= denom)
    return (neg? 0 : MINVAL(contlog_t));
  fracpart_t d_13579 = denom;
  fracpart_t quad[] = {0, 1, 1, 1};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  do {
    /* Update quad to shrink range containing the result */
    fracpart_t sum[8];
    for (int i = 0; i < 2; ++i) {
      int j = i ^ 2;
      dotprod(&sum[2*i], quad[i], quad[j], numer, d_13579);
      dotprod(&sum[2*j], quad[i], quad[j], 2*numer, 2*d_13579 - numer);
    }
    (void)pack(4, quad, sum);
    d_13579 += 2 * denom;
    ces.lo ^= 2;
  } while (!contlog_encode_bounds(&ces, quad));
  if (!neg)
    ces.arg = MINVAL(contlog_t) - ces.arg;
  return (ces.arg);
}

/* Compute cos(sqrt(x)) or sin(sqrt(x))/sqrt(x).
 */
contlog_t
contlog_cssqrt(contlog_t operand, int n)
{
  fracpart_t frac[2];
  int neg = contlog_decode(operand, frac);
  fracpart_t numer = neg ? -frac[0] : frac[0];
  fracpart_t denom = frac[1];
  fracpart_t L = 1;
  fracpart_t quad[] = {0, 1, denom, denom};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  int overflow = 0;
  int j = 0;
  do {
    /* Update quad to shrink range containing the result */
    fracpart_t sum[4];
    fracpart_t H = n++;
    H *= n++;
    dotprod2(&sum[0], overflow, quad[j^0], quad[j^2], L*numer, H);
    dotprod2(&sum[2], overflow, quad[j^1], quad[j^3], L*numer, H);
    overflow = pack(2, &quad[j], sum);
    L = H;
    j ^= 2;
    dotprod2(&sum[0], overflow, quad[j^0], quad[j^2], -numer, denom);
    dotprod2(&sum[2], overflow, quad[j^1], quad[j^3], -numer, denom);
    overflow += pack(2, &quad[j^2], sum);
  } while (!contlog_encode_bounds(&ces, quad));
  return (ces.arg);
}

/* Compute cos(sqrt(x)).
 */
contlog_t
contlog_cosqrt(contlog_t operand)
{
  return (contlog_cssqrt(operand, 1));
}

/* Compute sin(sqrt(x))/sqrt(x).
 */
contlog_t
contlog_sisqrt(contlog_t operand)
{
  return (contlog_cssqrt(operand, 2));
}
