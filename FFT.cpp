//
//  Fast Fourier Transform in N * (sum of factors of N) Time
//  D.P. Mitchell  1985/08/11.
//
#include "stdafx.h"
#include "Venera15.h"

static Complex zero;
static int AB;
static int prime[] = {
   4,   2,   3,   5,   7,  11,  13,  17,  19,  23,  29,  31,  37,  41,
  43,  47,  53,  59,  61,  67,  71,  73,  79,  83,  89,  97, 101, 103,
 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,   1
};
static Complex z[10000];

static void
slow_transform(Complex x0[], Complex x1[], int A, int B, Complex root)
{
	int aB, alphaB;
	int B2, B3, B4;
	Complex w0, w1, t, t1, t2, t3, t4, t5, m0, m1, m2, m3, m4, m5;

	switch (A) {

	case 2:
		x1[0].re = x0[0].re + x0[B].re;
		x1[0].im = x0[0].im + x0[B].im;
		x1[B].re = x0[0].re - x0[B].re;
		x1[B].im = x0[0].im - x0[B].im;
		break;
	case 3:
		B2 = B + B;
		t1.re = x0[B].re + x0[B2].re;
		t1.im = x0[B].im + x0[B2].im;
		m0.re = x0[0].re + t1.re;
		m0.im = x0[0].im + t1.im;
		m1.re = (root.re - 1.0)*t1.re + m0.re;
		m1.im = (root.re - 1.0)*t1.im + m0.im;
		m2.re = root.im * (x0[B2].im - x0[B].im);
		m2.im = root.im * (x0[B].re - x0[B2].re);
		x1[0] = m0;
		x1[B].re = m1.re + m2.re;
		x1[B].im = m1.im + m2.im;
		x1[B2].re = m1.re - m2.re;
		x1[B2].im = m1.im - m2.im;
		break;
	case 4:
		B2 = B + B;
		B3 = B2 + B;
		t1.re = x0[0].re + x0[B2].re;
		t1.im = x0[0].im + x0[B2].im;
		t2.re = x0[B].re + x0[B3].re;
		t2.im = x0[B].im + x0[B3].im;
		x1[0].re = t1.re + t2.re;
		x1[0].im = t1.im + t2.im;
		x1[B2].re = t1.re - t2.re;
		x1[B2].im = t1.im - t2.im;
		m2.re = x0[0].re - x0[B2].re;
		m2.im = x0[0].im - x0[B2].im;
		m3.re = x0[B3].im - x0[B].im;
		m3.im = x0[B].re - x0[B3].re;
		x1[B].re = m2.re + m3.re;
		x1[B].im = m2.im + m3.im;
		x1[B3].re = m2.re - m3.re;
		x1[B3].im = m2.im - m3.im;
		break;
	case 5:
		B2 = B + B;
		B3 = B2 + B;
		B4 = B3 + B;
		w0.re = root.re * root.re - root.im * root.im;
		w0.im = root.re * root.im;
		w0.im += w0.im;
		t1.re = x0[B].re + x0[B4].re;
		t1.im = x0[B].im + x0[B4].im;
		t2.re = x0[B2].re + x0[B3].re;
		t2.im = x0[B2].im + x0[B3].im;
		t3.re = x0[B].re - x0[B4].re;
		t3.im = x0[B].im - x0[B4].im;
		t4.re = x0[B3].re - x0[B2].re;
		t4.im = x0[B3].im - x0[B2].im;
		t5.re = t1.re + t2.re;
		t5.im = t1.im + t2.im;
		m0.re = x0[0].re + t5.re;
		m0.im = x0[0].im + t5.im;
		m1.im = m1.re = 0.5 * (root.re + w0.re) - 1.0;
		m1.re *= t5.re;
		m1.im *= t5.im;
		m2.re = m2.im = 0.5 * (root.re - w0.re);
		m2.re *= (t1.re - t2.re);
		m2.im *= (t1.im - t2.im);
		m3.re = root.im * (t3.im + t4.im);
		m3.im = root.im * (-t3.re - t4.re);
		m4.re = (root.im + w0.im) * t4.im;
		m4.im = - (root.im + w0.im) * t4.re;
		m5.re = - (root.im - w0.im) * t3.im;
		m5.im = (root.im - w0.im) * t3.re;
		t3.re = m3.re - m4.re;
		t3.im = m3.im - m4.im;
		t5.re = m3.re + m5.re;
		t5.im = m3.im + m5.im;
		t1.re = m0.re + m1.re;
		t1.im = m0.im + m1.im;
		t2.re = t1.re + m2.re;
		t2.im = t1.im + m2.im;
		t4.re = t1.re - m2.re;
		t4.im = t1.im - m2.im;
		x1[0]   = m0;
		x1[B4].re = t2.re + t3.re;
		x1[B4].im = t2.im + t3.im;
		x1[B3].re = t4.re + t5.re;
		x1[B3].im = t4.im + t5.im;
		x1[B2].re = t4.re - t5.re;
		x1[B2].im = t4.im - t5.im;
		x1[B].re = t2.re - t3.re;
		x1[B].im = t2.im - t3.im;
		break;
	default:
		t = zero;
		for (aB = 0; aB < AB; aB += B) {
			t.re += x0[aB].re;
			t.im += x0[aB].im;
		}
		x1[0] = t;
		for (w0 = root, alphaB = B; alphaB < AB; alphaB += B) {
			x1[alphaB] = x0[0];
			for (w1 = w0, aB = B; aB < AB; aB += B) {
				x1[alphaB].re += w1.re*x0[aB].re - w1.im*x0[aB].im;
				x1[alphaB].im += w1.re*x0[aB].im + w1.im*x0[aB].re;
				t.re = w1.re*w0.re - w1.im*w0.im;
				t.im = w1.re*w0.im + w1.im*w0.re;
				w1 = t;
			}
			t.re = root.re*w0.re - root.im*w0.im;
			t.im = root.re*w0.im + root.im*w0.re;
			w0 = t;
		}
	}
}

static void
fast_transform(Complex x[], int factors[], int A, int C, Complex root[], Complex slowroot[])
{
	int aC, bC, aBC, bAC, BC, AC, B;
	Complex w0, w1, t;

	B = A / *factors;
	BC = B * C;
	A = *factors++;
	AC = A * C;
	for (bC = 0; bC < BC; bC += C)
		slow_transform(x + bC, z + bC, A, BC, *slowroot);
	for (bC = bAC = 0; bC < BC; bC += C, bAC += AC)
		x[bAC] = z[bC];
	for (w0 = *root, aC = C, aBC = BC; aC < AC; aC += C, aBC += BC) {
		x[aC] = z[aBC];
		for (w1 = w0, bC = C, bAC = AC; bC < BC; bC += C, bAC += AC) {
			t = z[bC + aBC];
			x[aC + bAC].re = w1.re*t.re - w1.im*t.im;
			x[aC + bAC].im = w1.re*t.im + w1.im*t.re;
			t.re = w1.re*w0.re - w1.im*w0.im;
			t.im = w1.re*w0.im + w1.im*w0.re;
			w1 = t;
		}
		t.re = root->re*w0.re - root->im*w0.im;
		t.im = root->re*w0.im + root->im*w0.re;
		w0 = t;
	}
	if (B > 1) {
		root++;
		slowroot++;
		for (aC = 0; aC < AC; aC += C)
			fast_transform(x + aC, factors, B, AC, root, slowroot);
	}
}

void
fft(Complex x[], int n)
{
	int i, j, k, lastfactor;
	int factors[32];
	Complex root[32], slowroot[32];

	k = n;
	AB = n;
	lastfactor = 0;
	for (j = 0; k > 1; j++) {
		for (i = 0; k % prime[i]; i++)
			;
		factors[j] = prime[i];
		root[j].re = cos(D_2PI / double(k));
		root[j].im = sin(D_2PI / double(k));
		if (factors[j] == lastfactor)
			slowroot[j] = slowroot[j - 1];
		else {
			slowroot[j].re = cos(D_2PI / double(factors[j]));
			slowroot[j].im = sin(D_2PI / double(factors[j]));
		}
		lastfactor = factors[j];
		k = k / prime[i];
	}
	factors[j] = 1;
	fast_transform(x, factors, n, 1, root, slowroot);
}
