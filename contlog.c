#include <stdlib.h>

#include <string.h>
#include <stdio.h>
#include "contlog.h"


static void
debug_print(contlog_t operand, contlog_t box[], int nDims)
{
  int sh = 4*sizeof(contlog_t);
  __uintmax_t val = operand;
  val -= val >> sh >> sh << sh << sh;
  fprintf(stderr, " %0*jx", 2*sizeof(operand), val);
  for (int b = 0; b < (1 << nDims); ++b) {
    __intmax_t val = box[b];
    val -= val >> sh >> sh << sh << sh;
    fprintf(stderr, " %10jx", val);
  }
  fprintf(stderr, "\n");
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
      if (shift > 0) {
	addend += 1 << (shift - 1);
	addend >>= shift;
      }
      overflow = sum_overflows(addend, box[idx_nopd^b]);
    }
    shift += overflow;
    for (b = 0; b < bit_d; ++b) {
      if (shift > 0) {
	box[idx_n1^b] += 1 << (shift - 1);
	box[idx_n1^b] >>= shift;
	box[idx_nopd^b] >>= overflow;
      }
      box[idx_n1^b] += box[idx_nopd^b];
    }
  }

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
    operand ^= (contlog_t)1 << w;
    numer ^= 1;
    frac[numer] += frac[numer^1];
    unsigned int next_w = lobit(operand);
    frac[numer] <<= next_w - w - 1;
    w = next_w;
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
  int shift;
  if (operand >> (maxbits - 2)) {
    shift = maxbits - 2 - FLS(operand ^ (hibit - 1));
    operand <<= shift;
    operand ^= hibit;
  }
  else {
    shift = FLS(operand) - (maxbits - 2);
    operand <<= -shift;
  }

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
    if (shift >= 0) {
      int shift2 = FFS(frac[numer][1]) - 1;
      if (shift2 > shift)
	shift2 = shift;
      shift -= shift2;
      frac[numer][0] <<= shift;
      frac[numer][1] >>= shift2;
    }
    else {
      int shift2 = FFS(frac[numer][0]) - 1;
      if (shift2 > -shift)
	shift2 = -shift;
      shift += shift2;
      frac[numer][0] >>= shift2;
      frac[numer][1] <<= -shift;
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

  while (frac[numer^1] != 0 && w > 0) {
    int shift = lgratio(frac[numer], frac[numer^1]);
    if (shift > w)
      break;
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
  while (frac[numer^1] != 0 && w > 0) {
    int shift = lgratio(gmean.quot + mix, frac[numer^1]);
    if (shift > w)
      break;
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
