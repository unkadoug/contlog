#include "contlog.h"
#include <stdio.h>
#include <string.h>

contlog_t
sscan_frac(const char *s)
{
	contlog_t operand;
	if (strchr(s, '/')) {
		__intmax_t nx, dx;
		sscanf(s, "%jd/%jd", &nx, &dx);
		fracpart_t pair[] = {nx, dx};
		operand = contlog_encode_frac(pair);
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
	fracpart_t pair[2];
	contlog_decode_frac(operand, pair);
	__intmax_t nx = pair[0];
	__intmax_t dx = pair[1];
	__intmax_t opx = operand;
	int sh = 4*sizeof(contlog_t);
	opx -= opx >> sh >> sh << sh << sh;
	printf("%jd/%jd (%0*jx) = %18.12f\n", nx, dx,
	       (int)(2*sizeof(contlog_t)), opx,
	       (double)pair[0]/pair[1]);
}

static int usage(const char *cmd)
{
	printf("Usage: '%s x' or '%s x op` or '%s x op y', where:\n"
	       "x and y are rationals, expressed as fractions like 3/7 "
	       "or as hex representations,\n"
	       "and\n"
	       "op is from\n"
	       "'+' (arithmetic sum)\n"
	       "'-' (difference)\n"
	       "'*' (product x*y)\n"
	       "'~' (product with complement x*(1-y)\n"
	       "'/' (ratio)\n"
	       "'H' (harmonic sum 1/(1/x + 1/y))\n"
	       "'T' (arctangent sum (x+y)/(1-x*y))\n"
	       "'@' (pythagorean sum sqrt(x*x + y*y))\n"
	       "and where unary ops are\n"
	       "'/' (square root)\n"
	       "'e' (inverse exponential 1/e**x)\n"
	       "'l' (log1p ln(1+x))\n"
	       "'@' (sin(arctan(x)) x/sqrt(1+x*x)\n",
	       cmd, cmd, cmd);
	return 0;
}

int main(int argc, char *argv[])
{
#if 0
	fracpart_t pair[2];
	double last = 0.;
	unsigned i = 0x10000000;
	do {
		contlog_decode_frac(i, pair);
		double next = (double)pair[0]/pair[1];
		if (next <= last)
			printf("i: %u, %g, %g\n", i, last, next);
		unsigned j = contlog_encode_frac(pair);
		if (i != j)
			printf("i: %u, j: %u\n", i, j);
		last = next;
		i += 1;
	} while (i != 0xf0000000);
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
			printf("x/hypot(1,x): ");
			func_val = contlog_sinarctan(a);
			print_frac(func_val);
			break;

		case 'l':
			printf("log1p(x): ");
			func_val = contlog_log1p(a);
			print_frac(func_val);
			break;

		case 'e':
			printf("1/exp(x): ");
			func_val = contlog_expm(a);
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

	case '~':
		printf("x*(1-y): ");
		print_frac(contlog_compmult(a, b));
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
