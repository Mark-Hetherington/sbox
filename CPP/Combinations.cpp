#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define GET_BIT(x,b)	(((x)>>(b))&(1))

void init_comb(unsigned long long* state, unsigned long long k)
{
	unsigned long long i = 0;

	for(i=0;i<k;i++)
	{
		state[i] = 1<<i;
	}
}

// void print_state(unsigned long long* state,unsigned long long k)
// {
// 	unsigned long long j = 0, debug = 0;

// 	debug = 0;

// 	for(j=0;j<k;j++)
// 	{
// 		debug ^= state[j];
// 	}

// 	printf(">>> ");
// 	for(j=6;j<7;j--)
// 	{
// 		printf("%lld", GET_BIT(debug,j));
// 	}
// 	printf("\n");
// }

unsigned long long next_comb(unsigned long long* state, unsigned long long n, unsigned long long k)
{
	unsigned long long i = 0, ret = 0, j = 0, skip = 0;

	ret = 0;
	for(i=0;i<k;i++)
	{
		//printf("state[%lld] = %lld\n",i,state[i]);
		ret ^= state[i];
	}

	for(i=k-1;i<k;i--)
	{
		// printf("(i,n,k) = (%lld,%lld,%lld)\n",i,n,k);

		// print_state(state,k);

		if (!skip)
			state[i] <<= 1;

		// print_state(state,k);

		if ( (state[i] ^ ((unsigned long long)1<<(n-(k-1-i))) ) == 0 )
		{
			// printf(">>> first if (%lld) <<<\n",i);
			if(i == 0)
			{
				// printf(">>> second if <<<\n");
				state[i] <<= 1;

				if ( (state[i] ^ ((unsigned long long)1<<(n-(k-1-i)))) == 0 )
				{
					return 0;
				}
				else
				{
					for(j=1;j<k;j++)
					{
						state[j] = state[j-1] << 1;
					}
				}
			}
			else
			{
				// print_state(state,k);
				state[i-1] <<= 1;
				for(j=i;j<k;j++)
				{
					state[j] = state[j-1] << 1;
				}
				// print_state(state,k);
				skip = 1;
			}
		}
		else
		{
			break;
		}
	}

	// printf("=====\n");

	return ret;
}

int main()
{
	unsigned long long d = 0, i = 0, n = 6, k = 0, m = 8, j = 0;
	unsigned long long *state = NULL;

	state = (unsigned long long*)calloc(k,sizeof(unsigned long long));

	init_comb(state,k);
    
    for(i=0;i<20;i++)
    {
    	d = next_comb(state,n,k);

    	for(j=5;j<6;j--)
    	{
    		printf("%lld", GET_BIT(d,j));
    	}
    	printf("\n");
    }

    if (state)
    	free(state);

    return 0;
}