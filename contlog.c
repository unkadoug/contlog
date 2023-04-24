#include <stdlib.h>

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "contlog.h"

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

#define UNS(val) (sizeof(val) == sizeof(long long)?		\
		(unsigned long long) val :			\
		sizeof(val) == sizeof(long)?			\
		(unsigned long) val :				\
		sizeof(val) == sizeof(int)?			\
		(unsigned) val:					\
		sizeof(val) == sizeof(short)?			\
		  (unsigned short) val: (unsigned char)val)	\

static inline int
lobit(contlog_t operand)
{
  unsigned int invpos = 8*sizeof(contlog_t) - (SIGNED(contlog_t) ? 1 : 0);
  return (operand ? ffs(operand) - 1 : invpos);
}

/*
 * Translate operand into fraction frac[] = {denom, numer}.
 */
int
contlog_decode(contlog_t operand, contlog_t frac[])
{
  contlog_t hibit = (contlog_t)1 << SGNBIT_POS(contlog_t);
  int neg = 0;
  if (SIGNED(contlog_t)) {
    neg = (operand & hibit) != 0;
    operand ^= operand << 1;
    operand &= ~hibit;
  }
  else
    operand ^= operand << 1;
  frac[neg] = 1;
  frac[!neg] = 0;
  unsigned int shift;
  unsigned int numer = neg;
  unsigned int w = lobit(operand);
  while (operand != 0) {
    operand &= operand - 1;
    numer ^= 1;
    frac[numer] += frac[numer^1];
    shift = lobit(operand) - w;
    w += shift--;
    frac[numer] <<= shift;
  }
  shift = lobit(frac[0] | frac[1]);
  frac[0] >>= shift;
  frac[1] >>= shift;
  return neg;
}

static void
contlog_to_frac_ubound(contlog_t operand, contlog_t frac[])
{
  operand ^= (operand << 1) | 1;
  unsigned int numer = 1;
  unsigned int invpos = SGNBIT_POS(contlog_t) + (SIGNED(contlog_t) ? 0 : 1);
  unsigned int w = lobit(operand);
  frac[numer] = 1 << w;
  frac[numer^1] = 1;
  while (w < invpos) {
    operand &= operand - 1;
    numer ^= 1;
    frac[numer] += frac[numer^1];
    unsigned int shift = lobit(operand) - w;
    w += shift--;
    frac[numer] <<= shift;
  }
}

void
contlog_to_frac(contlog_t operand, contlog_t *n, contlog_t *d)
{
  contlog_t hibit = (contlog_t)1 << SGNBIT_POS(contlog_t);
  int neg = (operand & hibit) != 0;
  if (neg)
    operand = -operand;
  int numer = 1;
  contlog_t frac[2][2] = {{1, 0}, {0, 1}};

  if (operand == hibit)
    numer = 0;
  else if (operand != 0) {
    contlog_t bound[2][2];
    contlog_to_frac_ubound(operand-1, &bound[0][0]);
    contlog_to_frac_ubound(operand, &bound[1][0]);
    
    for (;;) {
      contlog_t val[2];
      for (int b = 0; b < 2; ++b) {
	val[b] = bound[b][numer] / bound[b][numer^1];
	bound[b][numer] %= bound[b][numer^1];
      }

      if (0 == operand % 2 && bound[numer^1][numer] == 0)
	/* The lower bound of the closed interval is found.  Make
	   sure the val fields are different to escape the loop.
	 */
	--val[numer^1];
      switch (val[numer] - val[numer^1]) {
      default:
	val[numer] = val[numer^1] + 1;
	break;
      case 0:
	if (bound[numer^1][numer] != 0)
	  break;
	/* The lower bound of the open interval is unavailable, so
	 * add a small extra cf term without passing the upper bound.
	 */
	for (int b = 0; b < 2; ++b)
	  frac[numer][b] += val[numer] * frac[numer^1][b];
	numer ^= 1;
	val[numer^1] = bound[numer^1][numer] / bound[numer^1][numer^1];
	val[numer] = val[numer^1] + 1;
	break;
      case 1:
	if (bound[numer][numer] != 0 ||
	    0 == operand % 2)
	  break;
	/* The upper bound of the open interval is unavailable, so
	 * add one or two small extra cf terms to the lower bound.
	 */
	for (int b = 0; b < 2; ++b)
	  frac[numer][b] += val[numer^1] * frac[numer^1][b];
	numer ^= 1;
	if (bound[numer][numer] >= 2 * bound[numer][numer^1])
	  val[numer^1] = 1;
	else {
	  for (int b = 0; b < 2; ++b)
	    frac[numer][b] += frac[numer^1][b];
	  numer ^= 1;
	  val[numer^1] = bound[numer^1][numer] / 
	    (bound[numer^1][numer^1] - bound[numer^1][numer]);
	}
	val[numer] = val[numer^1] + 1;
	break;
      }

      for (int b = 0; b < 2; ++b)
	frac[numer][b] += val[numer] * frac[numer^1][b];
      if (val[numer^1] < val[numer])
	break;
      numer ^= 1;
    }
  }
    
  *n = frac[numer][0];
  *d = frac[numer][1];
  
  if (neg)
    *n = -*n;
  }

static int
lgratio(contlog_t n, contlog_t d)
{
  if (d == 0)
    return 8 * sizeof(contlog_t);

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

contlog_t
frac_to_contlog(contlog_t n, contlog_t d)
{
  int lo = 0;
  int nbits = 8 * sizeof(contlog_t);

  if (SIGNED(contlog_t)) {
    if (n >> SGNBIT_POS(contlog_t)) {
      n = -n;
      lo ^= 1;
    }
    if (d >> SGNBIT_POS(contlog_t)) {
      d = -d;
      lo ^= 1;
    }
    else if (d == 0)
      lo = 1;
    --nbits;
  }

  contlog_t pair[2];
  pair[lo] = n;
  pair[lo^1] = d;
  contlog_t arg = n >= d;
  lo ^= !arg;
  for (;;) {
    int shift = lgratio(pair[lo], pair[lo^1]);
    shift = min(shift, nbits);
    pair[lo^1] <<= shift;
    nbits -= shift;
    arg <<= shift;
    if (nbits == 0)
      break;
    pair[lo] -= pair[lo^1];
    arg = 2 * arg + 2 * lo - 1;
    --nbits;
    lo ^= 1;
  }
  return (arg);
}

static contlog_t
hiprod(contlog_t a, contlog_t b)
{
  const unsigned int halfbits = 4 * sizeof(contlog_t);
  const contlog_t halfmask = ((contlog_t)1 << halfbits) - 1;
  contlog_t a_lo = a & halfmask;
  contlog_t a_hi = a >> halfbits;
  contlog_t b_lo = b & halfmask;
  contlog_t b_hi = b >> halfbits;
  contlog_t hilo = (((a_lo * b_lo) >> halfbits) & halfmask) + a_hi * b_lo;
  contlog_t lohi = (hilo & halfmask) + a_lo * b_hi;
  return ((hilo >> halfbits) + (lohi >> halfbits) + a_hi * b_hi);
}

static int
add_overflow(contlog_t add1, contlog_t add2, contlog_t sum)
{
  return (1 & (((add1 & add2) |
		(add1 & ~sum) |
		(add2 & ~sum)) >> SGNBIT_POS(contlog_t)));
}

static void
dotprod(contlog_t *sum,
	contlog_t x1, contlog_t x2,
	contlog_t y1, contlog_t y2)
{
  contlog_t prod1 = x1 * y1;
  contlog_t prod2 = x2 * y2;
  sum[0] = hiprod(x1, y1) + hiprod(x2, y2);
  sum[1] = prod1 + prod2;
  sum[0] += add_overflow(prod1, prod2, sum[1]);
}

static void
dotprod2(contlog_t *sum, int overflow,
	contlog_t x1, contlog_t x2,
	contlog_t y1, contlog_t y2)
{
  contlog_t prod1 = x1 * y1;
  contlog_t prod2 = x2 * y2;
  sum[0] = hiprod(x1, y1);
  if (overflow > 0) {
    int backflow = 8 * sizeof(contlog_t) - overflow;
    prod1 >>= overflow;
    prod1 &= ((contlog_t)1 << backflow) - 1;
    prod1 |= sum[0] << backflow;
    sum[0] >>= overflow;
  }
  sum[0] += hiprod(x2, y2);
  sum[1] = prod1 + prod2;
  sum[0] += add_overflow(prod1, prod2, sum[1]);
}

static void
oversum(contlog_t *sum, int overflow,
	contlog_t x1, contlog_t x2)
{
  int nbits = 8 * sizeof(contlog_t);
  sum[0] = x1 >> (nbits - 1);
  if (overflow > 0) {
    x1 >>= overflow;
    x1 &= ((contlog_t)1 << (nbits - overflow)) - 1;
    x1 |= sum[0] << (nbits - overflow);
  }
  sum[0] += (x2 >> (nbits - 1)) + add_overflow(x1, x2, x1 + x2);
  sum[1] = x1 + x2;
}

struct contlog_encode_state {
  contlog_t arg;
  int nbits;
  int lo;
};

static void
contlog_encode_state_init(struct contlog_encode_state *ces, contlog_t quad[])
{
  int nbits = 8 * sizeof(contlog_t);
  if (SIGNED(contlog_t))
    --nbits;
  contlog_t mask = 0;
  for (int i = 0; i < 4; ++i)
    mask |= quad[i] ^ (quad[i] >> 1);
  int shift = nbits - fls(mask);
  for (int i = 0; i < 4; ++i)
    quad[i] <<= shift;
  ces->nbits = nbits;
  ces->lo = 0;
  ces->arg = 0;
}

static void
contlog_encode_exact(struct contlog_encode_state *ces, contlog_t pair[])
{
  int nbits = ces->nbits;
  int lo = ces->lo & 1;
  contlog_t arg = ces->arg;
  if (-pair[lo] > 0) {
    /* result < 0 */
    int shift = ffs(arg) - 1;
    pair[lo] += pair[lo^1];
    pair[lo^1] -= pair[lo];
    pair[lo] += pair[lo^1];
    if (shift >= 0) {
      pair[lo] -= pair[lo^1] >> shift;
      if (SGNBIT_POS(contlog_t) != nbits + shift + 1)
	pair[lo^1] /= 2;
    }
    lo ^= 1;
  }
  if (pair[lo] > pair[lo^1]) {
    /* result > 1 */
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
  ces->lo = lo;
  ces->nbits = nbits;
  ces->arg = arg;
}

static int
contlog_encode_bounds(struct contlog_encode_state *ces, contlog_t quad[])
{
  contlog_t hibit = (contlog_t)1 << SGNBIT_POS(contlog_t);
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
	if (SGNBIT_POS(contlog_t) != nbits + shift + 1)
	  quad[i^1] /= 2;
      }
    }
    lo ^= 3;
  }
  if (quad[lo] > quad[lo^1]) {
    /* result > 1 */
    arg = 2 * (arg - (lo&1)) + 1;
    lo ^= 3;
  }

  /* Extract bits into arg until either arg is filled,
   * or lower and upper bound ratios have been pushed too far apart.
   */
  do {
    if (quad[lo] < 0 && quad[lo^1] / 2 < -quad[lo])
      break;
    /* -1/2 <= result */
    if (quad[lo^2] <= quad[lo^3] / 2) {
      /* -1/2 <= result <= 1/2 */
      int shift = min(lgratio(quad[lo^3], quad[lo^2]), nbits);
      if (quad[lo] < 0)
        shift = min(lgratio(quad[lo^1], -quad[lo]), shift);
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
pack(int n, contlog_t result[], contlog_t sum[])
{
  contlog_t mask = 0;
  for (int i = 0; i < n; ++i)
    mask |= sum[2*i] ^ ((sum[2*i] << 1) |
			((sum[2*i+1] >> SGNBIT_POS(contlog_t)) & 1));
  int overflow = fls(mask);
  mask = ((contlog_t)1 << 1 << (SGNBIT_POS(contlog_t) - overflow)) - 1;
  for (int i = 0; i < n; ++i)
    result[i] = (sum[2*i] << 1 << (SGNBIT_POS(contlog_t) - overflow)) |
	    ((sum[2*i+1] >> overflow) & mask);
  return (overflow);
}
static void
debug_print(contlog_t operand, contlog_t quad[], int j)
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

contlog_t
contlog_arith(contlog_t operand, contlog_t quad[])
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
  int j = 0;		/* posn in quad of coeff of 1 in numerator */
  if (operand & hibit)
    j ^= 2;
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

    /* Update quad to shrink range containing the result */
    contlog_t sum[4];
    oversum(&sum[0], overflow, quad[j^0] >> shift, quad[j^2]);
    oversum(&sum[2], overflow, quad[j^1] >> shift, quad[j^3]);
    overflow = pack(2, &quad[j], sum);
    (void)contlog_encode_bounds(&ces, quad);
  }
  if (operand == 0)
    contlog_encode_exact(&ces, &quad[j]);
  return (ces.arg);
}

/* Compute f(x) = sqrt(x) */
contlog_t
contlog_sqrt(contlog_t operand)
{
  contlog_t frac[2];
  if (contlog_decode(operand, frac))
    return (MINVAL(contlog_t));
  int ge_1 = frac[0] <= frac[1];
  contlog_t numer = frac[ge_1];
  contlog_t denom = frac[!ge_1];

  /*
   * Scale the argument down by a power of 4, to the range [1, 4), and the
   * result up by a power of 2, to avoid the higher iteration count that comes
   * with larger values.
   */
  int shift = lgratio(numer, denom) / 2;
  denom <<= 2 * shift;
  contlog_t quad[] = {2*denom, numer+denom, 1, 1};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  ces.nbits -= shift;
  
  int overflow = 0;
  int j = 0;
  for (;;) {
    /* Update quad to shrink range containing the result */
    contlog_t sum[4];
    j ^= 2;
    dotprod2(&sum[0], overflow,
	     quad[j^0], quad[j^2], numer-denom, 2);
    dotprod2(&sum[2], overflow,
	     quad[j^1], quad[j^3], numer-denom, 2);
    overflow = pack(2, &quad[j], sum);
    if (contlog_encode_bounds(&ces, quad))
      break;
    dotprod2(&sum[0], 0,
	     quad[j^0], quad[j^2], denom, 0);
    dotprod2(&sum[2], 0,
	     quad[j^1], quad[j^3], denom, 0);
    overflow += pack(2, &quad[j], sum);
  }
  if (ge_1)
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
  contlog_t frac[2];
  (void)contlog_decode(operand, frac);
  int proper = frac[1] <= frac[0];
  contlog_t numer = frac[!proper];
  contlog_t denom = frac[proper];
  contlog_t quad[] = {0, 1, denom, numer};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  int overflow = 0;
  int j = 0;
  do {
    /* Update quad to shrink range containing the result */
    contlog_t sum[4];
    dotprod2(&sum[0], overflow,
	     quad[j^0], quad[j^2], denom*denom, 2*numer);
    dotprod2(&sum[2], overflow,
	     quad[j^1], quad[j^3], denom*denom, 2*numer);
    overflow = pack(2, &quad[j], sum);
    j ^= 2;
  } while (!contlog_encode_bounds(&ces, quad));
  if (proper)
    ces.arg = contlog_div(ces.arg, operand);
  return (ces.arg);
}

static contlog_t
contlog_log1p_frac(contlog_t numer, contlog_t denom)
{
  contlog_t n_11234 = numer;
  contlog_t n_12345 = numer;
  contlog_t d_13579 = denom;
  contlog_t quad[] = {0, 1, 1, 0};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  do {
    /* Update quad to shrink range containing the result */
    contlog_t sum[8];
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
  contlog_t frac[2];
  int neg = contlog_decode(operand, frac);
  int numer = frac[1];
  int denom = frac[0];
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
    frac[1] *= shift;
    contlog_upshift(frac);
    contlog_t quad[] = {frac[1], frac[0], frac[0], 0};
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
  contlog_t frac[2];
  int neg = contlog_decode(operand, frac);
  if (frac[1] / (SGNBIT_POS(contlog_t) + 1) >= frac[0])
    return (neg? 0 : MINVAL(contlog_t));
  contlog_t numer = frac[1];
  contlog_t denom = frac[0];
  contlog_t d_13579 = denom;
  contlog_t quad[] = {0, 1, 1, 1};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  do {
    /* Update quad to shrink range containing the result */
    contlog_t sum[8];
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
  contlog_t frac[2];
  int neg = contlog_decode(operand, frac);
  contlog_t numer = neg ? -frac[1] : frac[1];
  contlog_t denom = frac[0];
  contlog_t L = 1;
  contlog_t quad[] = {0, 1, denom, denom};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  int overflow = 0;
  int j = 0;
  do {
    /* Update quad to shrink range containing the result */
    contlog_t sum[4];
    contlog_t H = n++;
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
  return contlog_cssqrt(operand, 1);
}

/* Compute sin(sqrt(x))/sqrt(x).
 */
contlog_t
contlog_sisqrt(contlog_t operand)
{
  return contlog_cssqrt(operand, 2);
}
