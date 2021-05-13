#include <stdlib.h>

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "contlog.h"

#define FLS(val) (sizeof(val) == sizeof(long long) ? flsll(val) :	\
		  sizeof(val) == sizeof(long) ? flsl(val) :		\
		  sizeof(val) == sizeof(int) ? fls(val) :		\
		  fls(val & ((1 << 8*sizeof(val)) - 1)))
#define FFS(val) (sizeof(val) == sizeof(long long) ? ffsll(val) :	\
		  sizeof(val) == sizeof(long) ? ffsl(val) :		\
		  ffs(val))

static void
debug_print(contlog_t operand, contlog_t box[], int nDims)
{
  int sh = 4*sizeof(contlog_t);
  __uintmax_t val = operand;
  val -= val >> sh >> sh << sh << sh;
  fprintf(stderr, " %0*jx", (int)(2*sizeof(operand)), val);
  for (int b = 0; b < (1 << nDims); ++b) {
    __intmax_t val = box[b];
    val -= val >> sh >> sh << sh << sh;
    fprintf(stderr, " %10jx", val);
  }
  fprintf(stderr, "\n");
}


static inline int
lobit(contlog_t operand)
{
  unsigned int invpos = 8*sizeof(contlog_t) - (SIGNED(contlog_t) ? 1 : 0);
  return (operand ? FFS(operand) - 1 : invpos);
}

/*
 * Translate operand into fraction frac[] = {denom, numer}.
 */
static int
contlog_decode(contlog_t operand, contlog_t frac[])
{
  const unsigned int maxbits = 8*sizeof(operand);
  int neg;
  if (SIGNED(contlog_t)) {
    neg = ((operand >> (maxbits-1)) & 1);
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
  unsigned int invpos = maxbits - (SIGNED(contlog_t) ? 1 : 0);
  unsigned int w = lobit(operand);
  while (operand != 0) {
    operand ^= operand & -operand;
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

void
contlog_load_arg(contlog_t operand, contlog_t frac[])
{
  const unsigned int maxbits = 8*sizeof(operand);
  int neg = contlog_decode(operand, frac);
  int shift = maxbits - FLS(frac[0] | frac[1]) - 1;
  if (neg)
    frac[1] = -frac[1];
  frac[0] <<= shift;
  frac[1] <<= shift;
}

contlog_t *
contlog_fold(contlog_t operand, contlog_t box[], int nDims)
{
  const unsigned int maxbits = 8*sizeof(contlog_t);
  contlog_t hibit = (contlog_t)1 << (maxbits - 1);
  unsigned int idx_n1;		/* posn in box of coeff of 1 in numerator */
  int b, d = nDims - 1;
  int bit_d = 1 << d;
  if (SIGNED(contlog_t)) {
    if (operand & hibit) {
      for (b = 0; b < bit_d; ++b)
	box[b^bit_d] = -box[b^bit_d];
    }
    operand ^= operand << 1;
    idx_n1 = (operand & hibit) ? bit_d : 0;
    operand <<= 1;
  }
  else {
    idx_n1 = (operand & hibit) ? bit_d : 0;
    operand ^= operand << 1;
  }

  while (operand != 0) {
    /* Swap the with-d face of the box with the without-d
     * face, and let idx_nopd identify the position of the constant
     * coefficient of the numerator in the new without-d face.
     */
    debug_print(operand, box, nDims);
    unsigned int idx_nopd = idx_n1;
    idx_n1 ^= bit_d;

    /* Find the leftmost set bit position of operand and shift the bit
       out. */
    int shift = maxbits - FLS(operand);
    operand <<= shift + 1;

    /* Shift the coeffs without operand and add to them the corresponding
     * coeffs with operand.  If necessary, divide everything by 2 first to
     * avoid overflow.
     */
    int overflow = 0;
    for (b = 0; !overflow && b < bit_d; ++b) {
      contlog_t addend = box[idx_n1^b];
      addend += (contlog_t)1 << shift >> 1;
      addend >>= shift;
      overflow = sum_overflows(addend, box[idx_nopd^b]);
    }
    shift += overflow;
    for (b = 0; b < bit_d; ++b) {
      box[idx_n1^b] += (contlog_t)1 << shift >> 1;
      box[idx_n1^b] >>= shift;
      box[idx_nopd^b] >>= overflow;
      box[idx_n1^b] += box[idx_nopd^b];
    }
  }
  debug_print(operand, box, nDims);

  return &box[idx_n1];
}

static void
contlog_to_frac_ubound(contlog_t operand, contlog_t frac[])
{
  const unsigned int maxbits = 8*sizeof(operand);
  operand ^= (operand << 1) | 1;
  unsigned int numer = 1;
  unsigned int invpos = maxbits - (SIGNED(contlog_t) ? 1 : 0);
  unsigned int w = lobit(operand);
  frac[numer] = 1 << w;
  frac[numer^1] = 1;
  while (w < invpos) {
    operand ^= operand & -operand;
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
  const unsigned int maxbits = 8*sizeof(operand);
  contlog_t hibit = (contlog_t)1 << (maxbits - 1);
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
    return 8 * sizeof(d);

  int lg = FLS(n) - FLS(d);
  if (n >= d << lg)
    ++lg;
  return (lg);
}

contlog_t
frac_to_contlog(contlog_t n, contlog_t d)
{
  contlog_t operand = 0;
  int w = 8 * sizeof(contlog_t);
  unsigned int neg = 0;

  if (SIGNED(contlog_t)) {
    if (d >> SGNBIT_POS(contlog_t)) {
      d = -d;
      neg = (n > 0);
    }
    else if (d != 0) {
      neg = (n >> SGNBIT_POS(contlog_t)) != 0;
    }
    else
      neg = 1;
    if (n >> SGNBIT_POS(contlog_t))
      n = -n;
    operand |= neg << --w;
  }

  contlog_t frac[2];
  frac[neg] = d;
  frac[neg^1] = n;
  unsigned int numer = frac[0] < frac[1];
  int shift;
  while (frac[numer^1] != 0 &&
	 (shift = lgratio(frac[numer], frac[numer^1])) <= w) {
    w -= shift;
    if (numer)
      operand |= (((contlog_t)1 << shift) - 1) << w;
    numer ^= 1;
    if (--shift)
      frac[numer] <<= shift;
    frac[numer^1] -= frac[numer];
  }
  if (numer)
    operand |= (contlog_t)1 << w;
  return operand;
}

/*-----------------------------------------------------------------*\
 * Compute the integer square root and remainder
 */
typedef struct {
  contlog_t quot;
  contlog_t rem;
} contlog_div_t;

static inline void
isqrt_help(contlog_div_t *res, contlog_t arg)
{
  int place = res->quot == 0 ?
    ((FLS(arg) + 1) / 2) :
    4 * sizeof(contlog_t);
  while (place-- > 0) {
    res->quot <<= 1;
    res->rem <<= 2;
    res->rem += (arg >> 2*place) & 3;
    if ((res->rem - 1) / 2 >= res->quot)
      res->rem -= 1 + 2 * res->quot++;
  }
}

static inline contlog_div_t
isqrt_prod(contlog_t a, contlog_t b)
{
  const unsigned int maxbits = 8*sizeof(contlog_t);
  const unsigned int halfbits = maxbits / 2;
  const contlog_t halfmask = ((contlog_t)1 << halfbits) - 1;
  contlog_t hibit = (contlog_t)1 << (maxbits - 1);

  contlog_t a_hi = a >> halfbits;
  contlog_t b_hi = b >> halfbits;
  a -= a_hi << halfbits;
  b -= b_hi << halfbits;

  contlog_t ab_hi = a_hi * b_hi;
  contlog_t ab_mid = a_hi * b + b_hi * a;
  contlog_t ab_lo = a * b;

  contlog_t carry = (ab_lo >> halfbits) & halfmask;
  ab_hi += ((carry + ab_mid) >> halfbits) & halfmask;
  ab_lo += ab_mid << halfbits;
  contlog_div_t res = {0, 0};
  isqrt_help(&res, ab_hi);
  isqrt_help(&res, ab_lo);
  return res;
}

contlog_t
contlog_sqrt(contlog_t operand)
{
  const unsigned int maxbits = 8*sizeof(contlog_t);
  contlog_t hibit = (contlog_t)1 << (maxbits - 1);
  if (operand < 0)
    return (hibit);

  contlog_t frac[2];
  (void)contlog_decode(operand, frac);
  unsigned int numer = frac[0] < frac[1];
  int shift = (lgratio(frac[numer], frac[numer^1]) - 1) / 2;
  frac[numer^1] <<= 2 * shift;

  int w = 8 * sizeof(contlog_t) - 1;
  operand = 0;
  w -= shift;
  if (numer)
    operand |= (((contlog_t)1 << shift) - 1) << w;

  shift = maxbits - FLS(frac[0] | frac[1]) - 1;
  if (shift >= 0) {
    frac[0] <<= shift;
    frac[1] <<= shift;
  }
  else {
    frac[0] >>= -shift;
    frac[0] &= ((contlog_t)1 << (maxbits + shift)) - 1;
    frac[1] >>= -shift;
    frac[1] &= ((contlog_t)1 << (maxbits + shift)) - 1;
  }
  contlog_div_t gmean = isqrt_prod(frac[0], frac[1]);
  contlog_t mix = 0;
  while (frac[numer^1] != 0 &&
	 (shift = lgratio(gmean.quot + mix, frac[numer^1])) <= w) {
    w -= shift;
    if (numer)
      operand |= (((contlog_t)1 << shift) - 1) << w;
    numer ^= 1;
    if (--shift) {
      frac[numer] <<= shift;
      frac[numer^1] >>= shift;
      frac[numer^1] &= ((contlog_t)1 << (maxbits - shift)) - 1;
    }
    frac[numer^1] -= (frac[numer] - mix) - mix;
    mix = frac[numer] - mix;
  }
  if (numer)
    operand |= (contlog_t)1 << w;
  return operand;
}

struct contlog_extractor {
  int w;
  int numer;
  int operand;
};

void
contlog_extractor_init(struct contlog_extractor *xtractor, int bits)
{
  xtractor->w = bits;
  xtractor->numer = 0;
  xtractor->operand = 0;
}

int
contlog_extract(struct contlog_extractor *xtractor, unsigned box[])
{
  int numer = xtractor->numer;
  int w = xtractor->w;
  int operand = xtractor->operand;

  /* Extract info from box to set operand bits */
  if (box[numer^2] <= box[numer^3])
    numer ^= 3;
  int shift;
  while ((shift = lgratio(box[numer], box[numer^1])) <= w && shift > 0) {
    box[numer^1] <<= shift - 1;
    box[numer^3] <<= shift - 1;
    int ival_spans_2 = box[numer^2] / 2 >= box[numer^3];
    if (ival_spans_2)
      --shift;
    w -= shift;
    if ((numer&1) == 0)
      operand  |= (((contlog_t)1 << shift) - 1) << w;
    if (ival_spans_2) {
      shift = 1;
      break;
    }
    numer ^= 3;
    assert(box[numer^1] >= box[numer]);
    assert(box[numer^3] >= box[numer^2]);
    box[numer^1] -= box[numer];
    box[numer^3] -= box[numer^2];
  }
  if (shift > w && (numer&1) == 0)
    operand |= (contlog_t)1 << w;
  xtractor->numer = numer;
  xtractor->w = w;
  xtractor->operand = operand;
  return (shift > w);
}

contlog_t
contlog_log1p(contlog_t operand)
{
  const unsigned int maxbits = 8 * sizeof(operand);
  contlog_t frac[2];
  int neg = contlog_decode(operand, frac);

  unsigned int x = frac[1];
  unsigned int y = frac[0];
  /* log1p(-x/y) == -log1p(x/(y-x)) */
  if (neg) {
    if (y <= x)
      return (contlog_t)1 << (maxbits - 1);
    y -= x;
  }
  unsigned box[] = {0, y, x, y};
  contlog_t ix = 0;
  int overflow = 0;
  struct contlog_extractor xtract;
  contlog_extractor_init(&xtract, maxbits - 1);
  for (int i = 2; !contlog_extract(&xtract, box); ++i) {
    /* Update box to shrink range containing the result */
    int j = (i&1) ? 2 : 0;
    if (j == 0)
      ix += x;
    unsigned long long quot = (j == 0) ? 2 * y : i;
    unsigned long long ixL = ix;
    unsigned long long box0 = ixL * box[j^0];
    unsigned long long box1 = ixL * box[j^1];
    if (overflow > 0) {
      box0 = (box0 + (1 << (overflow - 1))) >> overflow;
      box1 = (box1 + (1 << (overflow - 1))) >> overflow;
    }
    box0 += quot * box[j^2];
    box1 += quot * box[j^3];
    overflow = flsll(box0 | box1) - maxbits;
    if (overflow > 0) {
      box0 = (box0 + (1ULL << (overflow - 1))) >> overflow;
      box1 = (box1 + (1ULL << (overflow - 1))) >> overflow;
    }
    box[j^0] = box0;
    box[j^1] = box1;
    printf("%u/%u %u/%u\n", box[j^0], box[j^1], box[j^2], box[j^3]);
    printf("%g %g\n", (unsigned)box[j^0]/(double)box[j^1],
	   (unsigned)box[j^2]/(double)(unsigned)box[j^3]);
  }
  operand = xtract.operand;
    
  if (neg)
    operand = -operand;
  return operand;
}

/* Compute e**x by computing e**x/(1 + e**x) and shifting the result
 * one position left.
 */
contlog_t
contlog_exp(contlog_t operand)
{
  const unsigned int maxbits = 8 * sizeof(operand);
  contlog_t frac[2];
  int neg = contlog_decode(operand, frac);

  unsigned int x = frac[1];
  unsigned int y = frac[0];
  
  if (y == 0)
    return (contlog_t)1 << (maxbits - 1);

  unsigned long long ay = 6 * y;
  unsigned box[] = {2 * x + 2 * y, 2 * x + 4 * y, x + 2 * y, 4 * y};

  struct contlog_extractor xtract;
  contlog_extractor_init(&xtract, maxbits);

  while (!contlog_extract(&xtract, box)) {
    /* Update box to shrink range containing the result */
    unsigned long long box0 = box[0];
    unsigned long long box1 = box[1];
    unsigned long long box2 = box[2];
    unsigned long long box3 = box[3];
    box0 = ay * box2 + x * box0;
    box1 = ay * box3 + x * box1;
    ay += 4 * y;
    assert(box0 > x * box2);
    assert(box1 > x * box3);
    box2 = box0 - x * box2;
    box3 = box1 - x * box3;
    int overflow = flsll(box0 | box1) - maxbits;
    if (overflow > 0) {
      box0 = (box0 + (1ULL << (overflow - 1))) >> overflow;
      box1 = (box1 + (1ULL << (overflow - 1))) >> overflow;
      box2 = (box2 + (1ULL << (overflow - 1))) >> overflow;
      box3 = (box3 + (1ULL << (overflow - 1))) >> overflow;
    }
    box[0] = box0;
    box[1] = box1;
    box[2] = box2;
    box[3] = box3;
    xtract.numer ^= 2;
  }
  operand = xtract.operand;
  if (neg)
    operand = -operand ^ (1UL << (maxbits - 1));
  return operand;
}
