
typedef int contlog_t;
typedef int fracpart_t;


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
contlog_t contlog_recip_hypot1(contlog_t operand);
contlog_t contlog_hypot(contlog_t op0, contlog_t op1);
contlog_t contlog_log1p(contlog_t operand);
contlog_t contlog_exp(contlog_t operand);
contlog_t contlog_cosqrt(contlog_t operand);
contlog_t contlog_sisqrt(contlog_t operand);
