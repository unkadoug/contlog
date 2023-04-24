#include <string.h>

typedef int contlog_t;
#define SIGNED(T) (-(T)1 < 0)
#define SGNBIT_POS(T) (8*sizeof(T) - 1)
#define MINVAL(T) (SIGNED(T) ? (T)1 << SGNBIT_POS(T) : 0)
#define MAXVAL(T) (~MINVAL(T))


/* use continued logrithms, as described by Gosper,
 * but with logs a, b, c, d, ... represented as bits
 * 11^a 00^b 11^c 00^d ...
 * so that x < y for rationals x and y iff rep(x) < rep(y),
 * where rep(x) is the integer that stores the representation of x.
 */

contlog_t contlog_arith(contlog_t operand, contlog_t quad[]);
int contlog_decode(contlog_t operand, contlog_t frac[]);
void contlog_to_frac(contlog_t operand, contlog_t *n, contlog_t *d);
contlog_t frac_to_contlog(contlog_t n, contlog_t d);
contlog_t contlog_sqrt(contlog_t operand);
contlog_t contlog_recip_hypot1(contlog_t operand);
contlog_t contlog_log1p(contlog_t operand);
contlog_t contlog_exp(contlog_t operand);
contlog_t contlog_cosqrt(contlog_t operand);
contlog_t contlog_sisqrt(contlog_t operand);

static void
contlog_upshift(contlog_t frac[])
{
  int shift = SGNBIT_POS(contlog_t) - fls(frac[0] | frac[1]);
  frac[0] <<= shift;
  frac[1] <<= shift;
}

static contlog_t
contlog_add(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  if (neg)
    frac[1] = -frac[1];
  contlog_t quad[] = {frac[1], frac[0], frac[0], 0};
  return (contlog_arith(op0, quad));
}

static contlog_t
contlog_sub(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  if (!neg)
    frac[1] = -frac[1];
  contlog_t quad[] = {frac[1], frac[0], frac[0], 0};
  return (contlog_arith(op0, quad));
}

static contlog_t
contlog_mult(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  contlog_t quad[] = {0, frac[0], frac[1], 0};
  contlog_t val = contlog_arith(op0, quad);
  return (neg ? -val : val);
}

static contlog_t
contlog_div(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  contlog_t quad[] = {0, frac[1], frac[0], 0};
  contlog_t val = contlog_arith(op0, quad);
  return (neg ? -val : val);
}

static contlog_t
contlog_atnsum(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  if (neg)
    frac[1] = -frac[1];
  contlog_t quad[] = {frac[1], frac[0], frac[0], -frac[1]};
  return (contlog_arith(op0, quad));
}

static contlog_t
contlog_harmsum(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  int neg = contlog_decode(op1, frac);
  contlog_upshift(frac);
  if (neg)
    frac[1] = -frac[1];
  contlog_t quad[] = {0, frac[1], frac[1], frac[0]};
  return (contlog_arith(op0, quad));
}

static contlog_t
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
