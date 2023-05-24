
/* How big should the representation be? */
#if !defined(CONTLOG_BASE)
#define CONTLOG_BASE int
#endif

/* Should the numbers represented include all rationals (0), nonnegative
 * rationals (1) or unit-interval rationals (2) ?
 */
#if !defined(CONTLOG_RANGE)
#define CONTLOG_RANGE 0
#endif

#if (CONTLOG_RANGE == 0)
#define CONTLOG_SIGNED
#else
#define CONTLOG_SIGNED unsigned
#endif

typedef CONTLOG_SIGNED CONTLOG_BASE contlog_t;
typedef CONTLOG_BASE fracpart_t;

/* use continued logrithms, as described by Gosper,
 * but with logs a, b, c, d, ... represented as bits
 * 11^a 00^b 11^c 00^d ...
 * so that x < y for rationals x and y iff rep(x) < rep(y),
 * where rep(x) is the integer that stores the representation of x.
 */

void contlog_decode_frac(contlog_t operand, fracpart_t pair[]);
contlog_t contlog_encode_frac(fracpart_t pair[]);
contlog_t contlog_add(contlog_t op0, contlog_t op1);
contlog_t contlog_sub(contlog_t op0, contlog_t op1);
contlog_t contlog_mult(contlog_t op0, contlog_t op1);
contlog_t contlog_div(contlog_t op0, contlog_t op1);
contlog_t contlog_atnsum(contlog_t op0, contlog_t op1);
contlog_t contlog_harmsum(contlog_t op0, contlog_t op1);
contlog_t contlog_sqrt(contlog_t operand);
contlog_t contlog_sinarctan(contlog_t operand);
contlog_t contlog_hypot(contlog_t op0, contlog_t op1);
contlog_t contlog_log1p(contlog_t operand);
contlog_t contlog_expm(contlog_t operand);
contlog_t contlog_cosqrt(contlog_t operand);
contlog_t contlog_sisqrt(contlog_t operand);
