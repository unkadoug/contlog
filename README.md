# Contlog - representing rational numbers with continued logarithms

## Motivation

The standard way that computers represent numerical values that
include not just integers, but also other rational numbers, is with
floating point, a binary imitation of scientific notation that uses a
fixed number of bits to represent the sign and the significant digits
of a value, and the rest of the available bits to represent an
exponent, or scaling factor, to express where the radix point appears
relative to those significant bits.  How many bits to assign to each
purpose is standardized in the [IEEE floating
point](https://en.wikipedia.org/wiki/IEEE_754) formats for half-,
single- and double-precision calculations.  A bit is reserved for
expressing the sign of a value, which leads to distinct negative and
nonnegative zeroes.  A natural 'hole' in the range of representable
numbers is filled in by 'subnormal' numbers that cover a range where
exponents are too small to be represented.  There are also positive
and negative infinities, and the two 'not a number' numbers.  For all
its complexity, though, the floating point representation cannot
express the result of dividing 1 by 3, or 5, or anything that isn't a
power of two.

An alternative represention of rational numbers that avoids these
complexities and inaccuracies is one based on continued fractions,
where the binary encoding of a number is a description of the path
from the root of the [Stern-Brocot
tree](https://en.wikipedia.org/wiki/Stern%E2%80%93Brocot_tree) to that
number.  That representation, though, amounts to a unary, simple
bit-counting representation for integers, and limiting the largest
representable number to a value like 32 or 64 in a proposed
alternative to floating point makes that alternative impractical.

Instead, the software presented here implements a representation and
some basic arithmetic based on the ideas of continued binary
logarithms first expressed by Bill Gosper in his document describing
[arithmetic with continued
fractions](https://perl.plover.com/yak/cftalk/INFO/gosper.txt).
Continued logarithms allow any rational number to be represented
exactly with a finite number of bits, allows some large integers to be
expressed with no more bits than are necessary in their binary integer
expression, and avoids the strangeness around zero inherent in any
floating point representation.

## Introduction

Consider this version of the binary GCD algorithm, for finding the
greatest common divisor of two positive numbers using only shifting
and subtraction.

```
int gcd(int pair[2])
{
	int numer = pair[1] <= pair[0];
	int shift = ffs(pair[0] | pair[1]) - 1;
	pair[0] >>= shift;
	pair[1] >>= shift;
	while (pair[0] != pair[1]) {
		while (pair[numer] << 1 < pair[!numer]) {
			pair[numer] <<= 1;
			printf("%d", numer);
		}
		pair[!numer] -= pair[numer];
		printf("%d", numer);
		numer = !numer;
	}
	printf("1\n");
	return (pair[0] >> (ffs(pair[0]) - 1) << shift);
}
```

The print statements in the function log the transformations made to
the values that ultimately caused them to match, to end the iteration.
They also write a binary representation of the rational number that is
the ratio of the two inputs.  So, for example computing the gcd of 14
and 9 prints the sequence "101001", which is a binary representation
of the rational number 14/9 in the 'contlog' format.

Translating from a nonzero binary representation back to a ratio of
integers can be understood as walking a path from the root of an
infinite tree, the Gosper tree, that includes all the positive
rational numbers.  The root is labelled with the ratio 1/1, and the
final result lies somewhere in a range bounded by limits 0/1 and 1/0.
The walk ends when the bits unexamined consist of a leading 1-bit and
only 0-bits thereafter; thus "1" or "10" or "1000000" represents the
number 1/1.  Otherwise, processing a '0' bit means computing the next
node label by computing the mediant (sum of numerators over sum of
divisors) of the lower bound and the current node label, and then
doubling the values in the lower bound and replacing the upper bound
with the current node label.  Processing a '1' bit is similar,
computing a new label from the current label and the upper bound,
doubling the values of the upper bound, and replacing the lower bound
with the current value.  For the sequence "101001", the steps are:

```
1. low 0/1, current 1/1, high 1/0.  Process 1 bit.
2. low 1/1, current 2/1, high 2/0.  Process 0 bit.
3. low 2/2, current 3/2, high 2/1.  Process 1 bit.
4. low 3/2, current 5/3, high 4/2.  Process 0 bit.
5. low 6/4, current 8/5, high 5/3.  Process 0 bit.
6. low 12/8, current 14/9, high 8/5.  Done.
```

<details>
<summary>The first few levels of the Gosper tree, bounded
by 0/1 and 1/0</summary>

```
0/1
						1/32
					1/16----+
					|	1/12
				1/8-----+
				|	|	3/20
				|	1/6-----+
				|		1/5
			1/4-----+
			|	|		5/18
			|	|	3/10----+
			|	|	|	5/16
			|	1/3-----+
			|		|	3/8
			|		2/5-----+
			|			4/9
		1/2-----+
		|	|			9/17
		|	|		5/9-----+
		|	|		|	4/7
		|	|	3/5-----+
		|	|	|	|	8/13
		|	|	|	5/8-----+
		|	|	|		9/14
		|	2/3-----+
		|		|		5/7
		|		|	3/4-----+
		|		|	|	10/13
		|		4/5-----+
		|			|	6/7
		|			8/9-----+
		|				16/17
	1/1-----+
		|				17/16
		|			9/8-----+
		|			|	7/6
		|		5/4-----+
		|		|	|	13/10
		|		|	4/3-----+
		|		|		7/5
		|	3/2-----+
		|	|	|		14/9
		|	|	|	8/5-----+
		|	|	|	|	13/8
		|	|	5/3-----+
		|	|		|	7/4
		|	|		9/5-----+
		|	|			17/9
		2/1-----+
			|			9/4
			|		5/2-----+
			|		|	8/3
			|	3/1-----+
			|	|	|	16/5
			|	|	10/3-----+
			|	|		18/5
			4/1-----+
				|		5/1
				|	6/1-----+
				|	|	20/3
				8/1-----+
					|	12/1
					16/1----+
						32/1
1/0
```

</details>

In his notes on continued fraction arithmetic, Gosper introduced
continued logarithms as a binary alternative to continued fractions
for representing real numbers.  The representation here is only
slightly different from the one he proposed.  He observed that a
string of bits like

```
101001
```

can represent a sequence of exponents 0, 0, 0, 1, 0, that appear in
the continued fraction

$$2^0 + 2^0 / (
       2^0 + 2^0 / (
              2^0 + 2^0 / (
	             2^1 + 2^1 / (
		            2^0))))$$

which is 14/9.  An initial substring of exactly $n+1$ 1-bits shows
that the value is between $2^n$ and $2^{n+1}$.  After subtracting
$2^n$, and dividing by $2^n$, the value left is expressed in the bits
that follow, beginning with a string of one or more 0-bits.  Each
maximal set of matching bits represents the integer part of a binary
logarithm, thus the name 'continued logarithms'.

The representation here has the virtue that the mapping from a set of
bits, interpreted as a binary integer, to a rational number, is
monotonic.  That is, if integer i maps to rational number p, and
integer j maps to rational number q, then i < j if and only if p < q.
This is not a characteristic of Gosper's initially proposed
representation.

Strictly speaking, Gosper's encoding only represented numbers greater
than or equal to one.  The representation presented here can represent
proper fractions by beginning the bitstring with 0-bits instead of
1-bits.  If the bits of integer x represent the value y in contlog
format, then -x represents 1/y, when x is nonzero.  To extend this
representation further, to include negative numbers, add a leading
0-bit, and the twos-complement of the resulting bitstring represents
the negative of the initial value. This means that if the integer x
has bits set to represent the value y in contlog format, then -x
represents -y.  It also means that a reciprocal for zero can exist;
the integer with only its highest order bit set represents a negative
infinity, the only not-a-number that this format can represent.

Just as some langauges offer an integer type without a sign bit, to
represent a wider range of nonnegative integers, the sign bit can be
discarded from this format, for more accuracy in representing
nonnegtive numbers only.  The next bit, a 'reciprocal' bit, can also
be discarded for still more accuracy in representing numbers in the
unit interval, including 0 but excluding 1.

The farther from 1 that a number gets, the more bits it uses for the
initial scaling, and the fewer it has left for anything else.  Thus,
half the numbers that can be represented have absolute values at least
1, a quarter are at least 2, an eighth are at least 4, and so on.
With half the positive representable values between 1/2 and 2/1,
calculating with values in that range would likely preserve accuracy
better than would working with very large or small values.

## Extracting ratios from this representation

While this representation is great for representing exactly rational
numbers with small numerator and denominator in just a few bits, it
sacrifices the representation of large numbers, especially ones near
powers of two.  For example 127/1 is represented by

1111111011111101111101111011101101

so the signed contlog representation of 127 requires 36 bits.
Computing a 32-bit approximation to 127/1 requires rounding; this
implementation rounds to nearest, or to even integer representations
in case of a tie.

Translating that truncated bit string back to a rational number, leads
to the result 37722176/297025, which is not likely to be an expected
result.  To provide a more expected result, this implementation
considers a fixed size bit string to represent the interval that
includes all the ratios that would be converted to that bit string,
and picks the one to represent the interval that has the smallest
numerator and denominator.  For a fixed-size bit string that ends in a
zero, that includes the values at the endpoints of the interval, and
for the others, it does not.  Otherwise, two consecutive bit strings
could map to the same ratio.  For the case of 127/1, this produces the
expected result 127/1.  For 1000/999, the result looks like 1000/999,
and not 470999/470528.

## Arithmetic

Where Gosper seeks to compute with unbounded continued fractions, this
implementation aspires only to fixed-width operands and results.
Thus, it is feasible to completely decode one operand into a numerator
and denominator before considering the other operand at all, and this
implementation does that.  Those two values are scaled up, doubling
both until one or both of them would overflow if doubled again.  The
values are used to form a homographic function that is evaluated at
the value of the other operand to compute the composition of the two
operands.  So, for example, to compute x + y, first find that y=n/d, and then form

(n + dx) / (d + 0x)

and update the four coefficients as the program consumes the leading
bits of x.  When such an update leads to overflow, scale the
coefficents down, to fit within the fixed number of bits available.

After every incremental update, examine the coefficients to see if
some of the output bits can be extracted.


## Test program

The 'contlog' program produced from this source code is a simple
calculator that uses the contlog representation, and allows both
integer ratios and hexadecimal bitstrings to be used to represent
operands, and displays results in both forms, and also in a familiar
floating point form.  Sample usage:
```
$ contlog 4/7
x: 4/7 (26000000) =     0.571428571429
$ contlog 4/7 - 5/9
x: 4/7 (26000000) =     0.571428571429
y: 5/9 (24000000) =     0.555555555556
x-y: 1/63 (01042260) =     0.015873015873
$ contlog 55555555
x: 2178309/1346269 (55555555) =     1.618033988750
$ contlog 2/1 /
x: 2/1 (60000000) =     2.000000000000
sqrt(x): 8119/5741 (4e38e38e) =     1.414213551646
x/sqrt(x): 8119/5741 (4e38e38e) =     1.414213551646
```
