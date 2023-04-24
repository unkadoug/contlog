#include "contlog.h"
#include <stdio.h>

contlog_t
sscan_frac(const char *s)
{
  contlog_t operand;
  if (strchr(s, '/')) {
    __intmax_t nx, dx;
    sscanf(s, "%jd/%jd", &nx, &dx);
    contlog_t n = nx, d = dx;
    operand = frac_to_contlog(n, d);
  }
  else {
    __intmax_t opx;
    sscanf(s, "%jx", &opx);
    operand = opx;
    operand <<= 8*sizeof(contlog_t) - 4*strlen(s);
  }
  return operand;
}

void
print_frac(contlog_t operand)
{
  contlog_t n, d;
  contlog_to_frac(operand, &n, &d);
  __intmax_t nx = n;
  __intmax_t dx = d;
  __intmax_t opx = operand;
  int sh = 4*sizeof(contlog_t);
  opx -= opx >> sh >> sh << sh << sh;
  printf("%jd/%jd (%0*jx) = %18.12f\n", nx, dx, (int)(2*sizeof(contlog_t)), opx, (double)n/d);
}

static int usage(const char *cmd)
{
  printf("Usage: '%s x' or '%s x op` or '%s x op y', where:\n"
       "x and y are rationals, expressed as fractions like 3/7 or as hex representations,\n"
       "and\n"
       "op is from '+', '-', '*', '/', 'H' (harmonic sum), 'T' (arctangent sum), '@' (pythagorean sum).\n"
       "where unary ops are '/' (square root), 'e' (exponential), 'l' (log1p), '@' (1/sqrt(1+x*x)).\n",

       cmd, cmd, cmd);
    return 0;
}

int main(int argc, char *argv[])
{
#if 0
  unsigned i = 0x40000000;
  do {
    (void)contlog_log1p(i);
    i += 4;
  } while (i != 0x80000000);
#endif
  if (argc == 1)
    return usage(argv[0]);

  contlog_t a = sscan_frac(argv[1]);
  printf("x: ");
  print_frac(a);
  if (argc == 2)
    return 0;

  char op = argv[2][0];

  if (argc == 3) {
    contlog_t func_val;
    switch (op) {
    case '/':
      printf("sqrt(x): ");
      func_val = contlog_sqrt(a);
      print_frac(func_val);

      printf("x/sqrt(x): ");
      print_frac(contlog_div(a, func_val));
      break;

    case '@':
      printf("1/hypot(1,x): ");
      func_val = contlog_recip_hypot1(a);
      print_frac(func_val);
      break;

    case 'l':
      printf("log1p(x): ");
      func_val = contlog_log1p(a);
      print_frac(func_val);
      break;

    case 'e':
      printf("exp(x): ");
      func_val = contlog_exp(a);
      print_frac(func_val);
      break;

    case 'c':
      printf("cosqrt(x): ");
      func_val = contlog_cosqrt(a);
      print_frac(func_val);
      break;

    case 's':
      printf("sisqrt(x): ");
      func_val = contlog_sisqrt(a);
      print_frac(func_val);
      break;

    default:
      break;
    }
    return 0;
  }
  contlog_t b = sscan_frac(argv[3]);
  printf("y: ");
  print_frac(b);
  
  switch (op) {
  case '+':
    printf("x+y: ");
    print_frac(contlog_add(a, b));
    break;
    
  case '-':
    printf("x-y: ");
    print_frac(contlog_sub(a, b));
    break;
    
  case '*':
    printf("x*y: ");
    print_frac(contlog_mult(a, b));
    break;
    
  case '/':
    printf("x/y: ");
    print_frac(contlog_div(a, b));
    break;
    
  case 'T':
    printf("(x+y)/(1-xy): ");
    print_frac(contlog_atnsum(a, b));
    break;
    
  case 'H':
    printf("xy/(x+y): ");
    print_frac(contlog_harmsum(a, b));
    break;
    
  case '@':
    printf("sqrt(x*x+y*y): ");
    print_frac(contlog_hypot(a, b));
    break;
    
  default:
    break;
  }

  return 0;
}
