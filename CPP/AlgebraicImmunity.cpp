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
void init_comb(unsigned long long* state, unsigned long long state_length, unsigned long long k);
unsigned long long next_comb(unsigned long long* state, unsigned long long n, unsigned long long k);

// ***
int algebraic_immunity(unsigned long long *sbox, algebraic_properties* AP, unsigned long long n, unsigned long long m, unsigned long long sparseness)
{
    unsigned long int cols = 0;
    unsigned long long rows = 1<<n, AI = 0, d = 0, c = 0, N = m+n, v = 0, *SBOX = NULL, column = 0, *state = NULL, mask = 0;
    mzd_t *A = NULL, *X = NULL, *B = NULL;
    mzp_t *P = NULL, *Q = NULL;
    rci_t i = 0, j = 0, rank = 0;

    if(n+m > 64)
        return 1;

    state = (unsigned long long*)calloc(n,sizeof(unsigned long long)); // the state for combinations
    SBOX = (unsigned long long *)calloc(rows,sizeof(unsigned long long)); // current vector for combinations

    for(i=0;i<rows;i++)
    {
        SBOX[i] = (i<<m)^sbox[i]; // "<<m" because the bit length of sbox[i] is "m"
    }

    for(AI=0;AI<n;AI++)
    {
        cols = 0;

        for(c=0;c<=AI;c++)
        {
            cols += binomial(N,c);
        }

        A = mzd_init(rows, cols);
        P = mzp_init(rows);
        Q = mzp_init(cols);

        for(d = AI, column = 0; d != 0; d--)
        {
            init_comb(state,n,d);

            mask = next_comb(state,N,d);

            while(mask != 0)
            {
                for(j=0;j<rows;j++)
                {
                   mzd_write_bit(A,j,column,__builtin_popcountll(SBOX[j]&mask));
                }
                mask = next_comb(state,N,d);
                column++;
            } // while
        } // d
        
        for(j=0;j<rows;j++)
        {
            mzd_write_bit(A,j,column,1);
        }

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

            if (state)
                free(state);

            if (SBOX)
                free(SBOX);

            return 0;
        }
    }

    if (state)
    	free(state);

    if (SBOX)
        free(SBOX);

    return 0;
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
void init_comb(unsigned long long* state, unsigned long long state_length, unsigned long long k)
{
	unsigned long long i = 0;

	memset(state,0,state_length*sizeof(unsigned long long));

	for(i=0;i<k;i++)
	{
		state[i] = 1<<i;
	}
}
// ***
// n is limiten by 63
unsigned long long next_comb(unsigned long long* state, unsigned long long n, unsigned long long k)
{
	unsigned long long i = 0, ret = 0, j = 0, skip = 0;

	if( n > 63 )
	{
		return 0;
	}

	ret = 0;
	for(i=0;i<k;i++)
	{
		ret ^= state[i];
	}
	if(ret == 0)
		return 0;

	for(i=k-1;i<k;i--)
	{
		if (!skip)
			state[i] <<= 1;

		if ( (state[i] ^ ((unsigned long long)1<<(n-(k-1-i))) ) == 0 )
		{
			if(i == 0)
			{
				if ( (state[i] ^ ((unsigned long long)1<<(n-(k-1-i)))) == 0 )
				{
					memset(state,0,k*sizeof(unsigned long long));
					return ret;
				}
				else
				{
					state[i] <<= 1;

					for(j=1;j<k;j++)
					{
						state[j] = state[j-1] << 1;
					}
				}
			}
			else
			{
				state[i-1] <<= 1;
				for(j=i;j<k;j++)
				{
					state[j] = state[j-1] << 1;
				}
				skip = 1;
			}
		}
		else
		{
			break;
		}
	}

	return ret;
}
// ***