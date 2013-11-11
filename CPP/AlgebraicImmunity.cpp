/*
	Calculates algebraic immunity of the vectorial boolean function given by
	the component functions

AUTHORS:

- Oleksandr Kazymyrov (date in ISO year-month-day format): initial version

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include "AlgebraicImmunity.h"

unsigned long int  binomial(int n, int k);
int                twiddle(int *x, int *y, int *z, int *p);
void               inittwiddle(int m, int n, int *p);

// ***
void algebraic_immunity(unsigned long long *sbox, algebraic_properties* AP, unsigned long long n, unsigned long long m, unsigned long long sparseness)
{
    unsigned long int cols = 0;
    unsigned long long rows = 1<<n, AI = 0, d = 0, c = 0, N = m+n, v = 0, *SBOX = NULL, column;
    mzd_t *A = NULL, *X = NULL, *B = NULL;
    mzp_t *P = NULL, *Q = NULL;
    rci_t i = 0, j = 0, rank = 0;
    int x, y, z, *p = NULL, *b = NULL; // additional variables for "twiddle" (generation all combinations)

    p = (int*)calloc(N+2,sizeof(int)); // buffer for generating all combinationss
    b = (int*)calloc(N,sizeof(int));   // current vector for combinations

    SBOX = (unsigned long long *)calloc(rows,sizeof(unsigned long long));   // current vector for combinations

    for(i=0;i<rows;i++)
    {
        SBOX[i] = (i<<m)^sbox[i]; // "<<m" because the bit length of sbox[i] is "m"
    }

    for(AI=1;AI<n;AI++)
    {
        cols = 0;

        for(c=0;c<=AI;c++)
        {
            cols += binomial(N,c);
        }

        A = mzd_init(rows, cols);
        P = mzp_init(rows);
        Q = mzp_init(cols);

        for(j=0;j<rows;j++)
        {
            mzd_write_bit(A,j,0,1);
        }

        for(d=AI,column=1;d != 0;d--)
        {
            inittwiddle(d, N, p);

            for(i = 0; i != N-d; i++)
            {
                b[i] = 0;
            }

            while(i != N)
            {
                b[i++] = 1;
            }

            for(j=0;j<rows;j++)
            {
                v = 1;
                for(i = 0; i != N; i++)
                    if (b[i])
                    {
                        v &= GET_BIT(SBOX[j],i);
                    }
               mzd_write_bit(A,j,column,v);
            }

            column++;

            while(!twiddle(&x, &y, &z, p))
            {
                b[x] = 1;
                b[y] = 0;

                for(j=0;j<rows;j++)
                {
                    v = 1;
                    for(i = 0; i != N; i++)
                        if (b[i])
                        {
                            v &= GET_BIT(SBOX[j],i);
                        }
                   mzd_write_bit(A,j,column,v);
                }
                column++;
            } // while
        } // d

        B = mzd_copy(NULL, A);
        rank = mzd_ple(A, P, Q, 0);

        mzd_free(A);
        mzp_free(P);
        mzp_free(Q);

        if ((cols-rank) > 0)
        {
            AP->AI = AI;
            AP->NE = cols-rank;

            if (sparseness == 1)
            {
                X = mzd_kernel_left_pluq(B,0);
                AP->SP = 1 - mzd_density(X,0);

                mzd_free(X);
            }

            mzd_free(B);

            break;
        }
        else
        {
            mzd_free(B);
        }

        if (AI > (n>>1))
        {
            perror(">>> You have found a bug in 'algebraic_immunity' <<<");

            if (SBOX)
                free(SBOX);
            if (p)
                free(p);
            if (b)
                free(b);

            return;
        }
    }

    if (SBOX)
        free(SBOX);
    if (p)
        free(p);
    if (b)
        free(b);
}
// ***
void factorial(mpz_t result, unsigned long input) {
    mpz_set_ui(result, 1);
    while (input > 1) {
        mpz_mul_ui(result, result, input--);
    }
}
// ***
unsigned long int binomial(int n, int k)
{
	unsigned long int ret = 0;
	mpz_t tmpN,tmpK,tmpNK;

    if (k>n)
        return 0;

	mpz_init(tmpN);
	mpz_init(tmpK);
	mpz_init(tmpNK);

	factorial(tmpN,n);
	factorial(tmpK,k);
	factorial(tmpNK,n-k);

	mpz_mul(tmpNK,tmpK,tmpNK);
	mpz_cdiv_q(tmpN,tmpN,tmpNK);

	ret = mpz_get_ui(tmpN);

	mpz_clear(tmpN);
	mpz_clear(tmpK);
	mpz_clear(tmpNK);

	return ret;
}
// ***
int twiddle(int *x, int *y, int *z, int *p)
{
	register int i, j, k;
	j = 1;

	while(p[j] <= 0)
		j++;

	if(p[j-1] == 0)
	{
		for(i = j-1; i != 1; i--)
			p[i] = -1;

		p[j] = 0;
		*x = *z = 0;
		p[1] = 1;
		*y = j-1;
	}
	else
	{
		if(j > 1)
			p[j-1] = 0;
		do
		{
			j++;
		}
		while(p[j] > 0);

		k = j-1;
		i = j;

		while(p[i] == 0)
			p[i++] = -1;

		if(p[i] == -1)
		{
			p[i] = p[k];
			*z = p[k]-1;
			*x = i-1;
			*y = k-1;
			p[k] = -1;
		}
		else
		{
			if(i == p[0])
				return(1);
			else
			{
				p[j] = p[i];
				*z = p[i]-1;
				p[i] = 0;
				*x = j-1;
				*y = i-1;
			}
		}
	}

	return(0);
}
// ***
void inittwiddle(int m, int n, int *p)
{
	int i;
	p[0] = n+1;

	for(i = 1; i != n-m+1; i++)
		p[i] = 0;

	while(i != n+1)
	{
		p[i] = i+m-n;
		i++;
	}

	p[n+1] = -2;

	if(m == 0)
		p[1] = 1;
}
// ***