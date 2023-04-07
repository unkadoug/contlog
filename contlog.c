#include <stdlib.h>

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "contlog.h"

#define ffs(X) _Generic((X), \
		int: ffs, \
		long: ffsl, \
		long long: ffsll \
		)(X)
#define fls(X) _Generic((X), \
		int: fls, \
		long: flsl, \
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

#define DEBUG
static void
debug_print(contlog_t operand, contlog_t quad[], int nDims)
{
#ifdef DEBUG
  int sh = 4*sizeof(contlog_t);
  __uintmax_t val = operand;
  val -= val >> sh >> sh << sh << sh;
  fprintf(stderr, " %0*jx", (int)(2*sizeof(operand)), val);
  for (int b = 0; b < (1 << nDims); ++b) {
    __intmax_t val = quad[b];
    val -= val >> sh >> sh << sh << sh;
    fprintf(stderr, " %10jx", val);
  }
  fprintf(stderr, "\n");
#endif
}


static inline int
lobit(contlog_t operand)
{
  unsigned int invpos = 8*sizeof(contlog_t) - (SIGNED(contlog_t) ? 1 : 0);
  return (operand ? ffs(operand) - 1 : invpos);
}

/*
 * Translate operand into fraction frac[] = {denom, numer}.
 */
static int
contlog_decode(contlog_t operand, contlog_t frac[])
{
  int neg;
  if (SIGNED(contlog_t)) {
    neg = ((operand >> SGNBIT_POS(contlog_t)) & 1);
    operand ^= operand << 1;
    operand &= ~((contlog_t)1 << SGNBIT_POS(contlog_t));
  }
  else {
    neg = 0;
    operand ^= operand << 1;
  }
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
contlog_upshift(contlog_t frac[])
{
  int shift = SGNBIT_POS(contlog_t) - fls(frac[0] | frac[1]);
  frac[0] <<= shift;
  frac[1] <<= shift;
}

void
contlog_load_arg(contlog_t operand, contlog_t frac[])
{
  int neg = contlog_decode(operand, frac);
  contlog_upshift(frac);
  if (neg)
    frac[1] = -frac[1];
}

contlog_t
contlog_fold(contlog_t operand, contlog_t quad[])
{
  contlog_t hibit = MINVAL(contlog_t);
  int j = 0;		/* posn in quad of coeff of 1 in numerator */
  if (SIGNED(contlog_t)) {
    if (operand & hibit) {
      for (int i = 0; i < 2; ++i)
	quad[i^2] = -quad[i^2];
      j ^= 2;
    }
    operand <<= 1;
  }
  if (operand & hibit)
    j ^= 2;
  operand ^= operand << 1;

  while (operand != 0) {
    /* Swap the with-d face of the quad with the without-d
     * face, and let idx_nopd identify the position of the constant
     * coefficient of the numerator in the new without-d face.
     */
    debug_print(operand, quad, 2);

    /* Find the leftmost set bit position of operand and shift the bit
       out. */
    int shift = SGNBIT_POS(contlog_t) - (fls(operand) - 1);
    operand <<= shift + 1;

    /* Shift the coeffs without operand and add to them the corresponding
     * coeffs with operand.  If necessary, divide everything by 2 first to
     * avoid overflow.
     */
    int overflow = 0;
    for (int i = 0; !overflow && i < 2; ++i) {
      contlog_t addend = quad[i^j^2];
      addend += (contlog_t)1 << shift >> 1;
      addend >>= shift;
      overflow = sum_overflows(addend, quad[i^j]);
    }
    shift += overflow;
    for (int i = 0; i < 2; ++i) {
      quad[i^j^2] += (contlog_t)1 << shift >> 1;
      quad[i^j^2] >>= shift;
      quad[i^j] >>= overflow;
      quad[i^j^2] += quad[i^j];
    }
    j ^= 2;
  }
  debug_print(operand, quad, 2);

  return (frac_to_contlog(quad[j+1], quad[j]));
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
  contlog_t hibit = MINVAL(contlog_t);
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
  }

  contlog_t pair[2];
  pair[lo] = n;
  pair[lo^1] = d;
  int nbits = SGNBIT_POS(contlog_t);
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

static void
dotprod(contlog_t *sum,
	contlog_t x1, contlog_t x2,
	contlog_t y1, contlog_t y2)
{
  contlog_t prod1 = x1 * y1;
  contlog_t prod2 = x2 * y2;
  sum[0] = hiprod(x1, y1) + hiprod(x2, y2);
  sum[1] = prod1 + prod2;
  sum[0] += 1 & (((prod1 & prod2) |
		 (prod1 & ~sum[1]) |
		 (prod2 & ~sum[1])) >> SGNBIT_POS(contlog_t));
  long ans = (long)x1 * y1 + (long)x2 * y2;
  if (sum[1] != (contlog_t)ans ||
      sum[0] != (contlog_t)(ans >> 32))
	  printf("Stop!\n");
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
  sum[0] += 1 & (((prod1 & prod2) |
		 (prod1 & ~sum[1]) |
		 (prod2 & ~sum[1])) >> SGNBIT_POS(contlog_t));
#if 1
  long ans = (((long)x1 * y1) >> overflow)  + (long)x2 * y2;
  if (sum[1] != (contlog_t)ans ||
      sum[0] != (contlog_t)(ans >> 32))
	  printf("Stop!\n");
#endif
}

struct contlog_encode_state {
  contlog_t arg;
  int nbits;
  int lo;
};

void
contlog_encode_state_init(struct contlog_encode_state *ces, contlog_t quad[])
{
  contlog_t mask = 0;
  for (int i = 0; i < 4; ++i)
    mask |= quad[i];
  int shift = SGNBIT_POS(contlog_t) - fls(mask);
  for (int i = 0; i < 4; ++i)
    quad[i] <<= shift;
  ces->nbits = SGNBIT_POS(contlog_t);
  ces->lo = 0;
  ces->arg = 0;
}

int
contlog_encode(struct contlog_encode_state *ces, contlog_t quad[])
{
  int nbits = ces->nbits;
  int lo = ces->lo;
  contlog_t arg = ces->arg;
  if (quad[lo] >= 0 && quad[lo^2] <= 0) {
	  printf("splat!\n");
  }
  if (quad[lo^2] <= 0) {
    /* result < 0 */
    int shift = ffs(arg) - 1;
    for (int i = lo&1; i < 4; i += 2) {
      quad[i] += quad[i^1];
      quad[i^1] -= quad[i];
      quad[i] += quad[i^1];
      if (SGNBIT_POS(contlog_t) != nbits) {
	quad[i] -= quad[i^1] >> shift;
	if (SGNBIT_POS(contlog_t) != nbits + shift + 1)
	  quad[i^1] /= 2;
      }
    }
    lo ^= 3;
  }
  if (quad[lo] >= quad[lo^1]) {
    /* result >= 1 */
    arg = 1;
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
  int shift = 0;//lgratio(numer, denom) / 2;
  denom <<= 2 * shift;
  contlog_t quad[] = {2*denom, numer+denom, 1, 1};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  ces.nbits -= shift;
  int overflow = 0;
  int j = 0;
  do {
    /* Update quad to shrink range containing the result */
    contlog_t sum[4];
    j ^= 2;
    for (int i = 0; i < 2; ++i)
      dotprod2(&sum[2*i], overflow,
	       quad[i^j], quad[i^j^2], (numer-denom)*denom, 2*denom);
    overflow = pack(2, &quad[j], sum);
  } while (!contlog_encode(&ces, quad));
  if (ge_1)
    ces.arg = MINVAL(contlog_t) - ces.arg;
  return (ces.arg);
}

/* Compute f(x) = 1 / sqrt(1 + x*x) */
contlog_t
contlog_recip_hypot1(contlog_t operand)
{
  if (operand < 0)
    operand = -operand;
  contlog_t frac[2];
  (void)contlog_decode(operand, frac);
  contlog_t numer = frac[1];
  contlog_t denom = frac[0];
  contlog_t quad[] = {0, 1, denom, numer};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  int overflow = 0;
  int j = 0;
  do {
    /* Update quad to shrink range containing the result */
    contlog_t sum[4];
    for (int i = 0; i < 2; ++i)
      dotprod2(&sum[2*i], overflow,
	       quad[i^j], quad[i^j^2], denom*denom, 2*numer);
    overflow = pack(2, &quad[j], sum);
    j ^= 2;
  } while (!contlog_encode(&ces, quad));
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
  } while (!contlog_encode(&ces, quad));
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
    contlog_t quad[] = {frac[0], frac[1], 0, frac[0]};
    arg = contlog_fold(arg, quad);
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
  } while (!contlog_encode(&ces, quad));
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
  contlog_t H = 1, L;
  contlog_t quad[] = {0, 1, 1, 1};
  struct contlog_encode_state ces;
  contlog_encode_state_init(&ces, quad);
  int overflow = 0;
  int j = 0;
  do {
    /* Update quad to shrink range containing the result */
    contlog_t sum[4];
    L = H;
    H = denom * n++;
    H *= n++;
    for (int i = 0; i < 2; ++i)
      dotprod2(&sum[2*i], overflow, quad[i^j], quad[i^j^2], numer*L, H-numer);
    overflow = pack(2, &quad[j], sum);
    j ^= 2;
  } while (!contlog_encode(&ces, quad));
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
