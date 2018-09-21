#include "contlog.h"

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
  printf("%jd/%jd (%0*jx) = %18.12f\n", nx, dx, 2*sizeof(contlog_t), opx, (double)n/d);
}

static int usage(const char *cmd)
{
  printf("Usage: '%s x' or '%s x op` or '%s x op y', where:\n"
       "x and y are rationals, expressed as fractions like 3/7 or as hex representations,\n"
       "and\n"
       "op is from '+', '-', '*', '/', '\\' (backward division), 'H' (harmonic sum), 'T' (arctangent sum).\n"
       "where unary ops are '+' (increment) and '/' (square root).\n",

       cmd, cmd, cmd);
    return 0;
}

int main(int argc, char *argv[])
{
  if (argc == 1)
    return usage(argv[0]);

  int a = sscan_frac(argv[1]);
  printf("x: ");
  print_frac(a);
  if (argc == 2)
    return 0;

  char op = argv[2][0];

  if (argc == 3) {
    switch (op) {
    case '+':
      printf("x+1: ");
      print_frac(contlog_incr(a));
      break;

    case '/':
      printf("sqrt(x): ");
      contlog_t r = contlog_sqrt(a);
      print_frac(r);

      printf("x/sqrt(x): ");
      print_frac(contlog_div(a, r));
      break;

    default:
      break;
    }
    return 0;
  }
  int b = sscan_frac(argv[3]);
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
    
  case '\\':
    printf("x\\y: ");
    print_frac(contlog_backdiv(a, b));
    break;
    
  case 'T':
    printf("(x+y)/(1-xy): ");
    print_frac(contlog_atnsum(a, b));
    break;
    
  case 'H':
    printf("xy/(x+y): ");
    print_frac(contlog_harmsum(a, b));
    break;
    
  default:
    break;
  }

  return 0;
}
