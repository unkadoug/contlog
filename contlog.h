/* Use continued logrithms, as described by Gosper
 * (https://perl.plover.com/yak/cftalk/INFO/gosper.txt) to implement basic
 * arithmetic and a few useful functions.
 *
 * Author: Doug Moore (unkadoug@gmail.com)
 * Use it, share it, make it replace floating point everywhere. Just keep my
 * name in it somewhere.
 */

/* How big should the representation be? */
#if !defined(CONTLOG_BASE)
#define CONTLOG_BASE int
#endif

/* Should the numbers represented include negative numbers? */
#if !defined(CONTLOG_SIGNED)
#define CONTLOG_SIGNED 1	/* yes, by default */
#endif
#if CONTLOG_SIGNED == 1
#define CONTLOG_SIGNED_MOD signed
#else
#define CONTLOG_SIGNED_MOD unsigned
#endif

/* Should the numbers represented include values outside the unit interval? */
#if !defined(CONTLOG_UNBOUNDED)
#define CONTLOG_UNBOUNDED 1	/* yes, by default */
#endif

typedef CONTLOG_SIGNED_MOD CONTLOG_BASE contlog_t;
typedef unsigned CONTLOG_BASE ufracpart_t;
typedef CONTLOG_BASE fracpart_t;

int contlog_decode_frac(contlog_t operand, ufracpart_t pair[]);
contlog_t contlog_encode_frac(fracpart_t pair[]);
contlog_t contlog_add(contlog_t op0, contlog_t op1);
contlog_t contlog_sub(contlog_t op0, contlog_t op1);
contlog_t contlog_mult(contlog_t op0, contlog_t op1);
contlog_t contlog_compmult(contlog_t op0, contlog_t op1);
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
