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

static inline int
sum_overflows(contlog_t a, contlog_t b)
{
  return (SIGNED(contlog_t) ? (a >> SGNBIT_POS(contlog_t) == 0) : 1)?
    (b > MAXVAL(contlog_t)-a) : (b < MINVAL(contlog_t)-a);
}

contlog_t contlog_fold(contlog_t operand, contlog_t quad[]);
void contlog_load_arg(contlog_t operand, contlog_t frac[]);
void contlog_to_frac(contlog_t operand, contlog_t *n, contlog_t *d);
contlog_t frac_to_contlog(contlog_t n, contlog_t d);
contlog_t contlog_sqrt(contlog_t operand);
contlog_t contlog_recip_hypot1(contlog_t operand);
contlog_t contlog_log1p(contlog_t operand);
contlog_t contlog_exp(contlog_t operand);
contlog_t contlog_cosqrt(contlog_t operand);
contlog_t contlog_sisqrt(contlog_t operand);

static contlog_t
contlog_incr(contlog_t operand)
{
  contlog_t frac[2];
  contlog_load_arg(operand, frac);
  if (sum_overflows(frac[0], frac[1])) {
    frac[1] = (frac[1] + 1) >> 1;
    frac[0] >>= 1;
  }
  frac[1] += frac[0];
  return frac_to_contlog(frac[1], frac[0]);
}

static contlog_t
contlog_add(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  contlog_load_arg(op1, frac);
  contlog_t quad[] = {frac[0], frac[1], 0, frac[0]};
  return (contlog_fold(op0, quad));
}

static contlog_t
contlog_sub(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  contlog_load_arg(op1, frac);
  contlog_t quad[] = {frac[0], -frac[1], 0, frac[0]};
  return (contlog_fold(op0, quad));
}

static contlog_t
contlog_mult(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  contlog_load_arg(op1, frac);
  contlog_t quad[] = {frac[0], 0, 0, frac[1]};
  return (contlog_fold(op0, quad));
}

static contlog_t
contlog_div(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  contlog_load_arg(op1, frac);
  contlog_t quad[] = {frac[1], 0, 0, frac[0]};
  return (contlog_fold(op0, quad));
}

static contlog_t
contlog_backdiv(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  contlog_load_arg(op1, frac);
  contlog_t quad[] = {0, frac[1], frac[0], 0};
  return (contlog_fold(op0, quad));
}

static contlog_t
contlog_atnsum(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  contlog_load_arg(op1, frac);
  contlog_t quad[] = {frac[0], frac[1], -frac[1], frac[0]};
  return (contlog_fold(op0, quad));
}

static contlog_t
contlog_harmsum(contlog_t op0, contlog_t op1)
{
  contlog_t frac[2];
  contlog_load_arg(op1, frac);
  contlog_t quad[] = {frac[1], 0, frac[0], frac[1]};
  return (contlog_fold(op0, quad));
}

static contlog_t
contlog_hypot(contlog_t op0, contlog_t op1)
{
  if (op0 < 0)
    op0 = -op0;
  if (op1 < 0)
    op1 = -op1;
  if (op0 < op1) {
    contlog_t tmp = op0;
    op0 = op1;
    op1 = tmp;
  }
  return (contlog_div(op1, contlog_recip_hypot1(contlog_div(op0, op1))));
}
