#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GET_BIT(x,b)	(((x)>>(b))&(1))

void init_comb(unsigned long long* state, unsigned long long state_length, unsigned long long k)
{
	unsigned long long i = 0;

	memset(state,0,state_length*sizeof(unsigned long long));

	for(i=0;i<k;i++)
	{
		state[i] = 1<<i;
	}
}

// void print_state(unsigned long long* state, unsigned long long k)
// {
// 	unsigned long long j = 0, debug = 0;

// 	debug = 0;

// 	for(j=0;j<k;j++)
// 	{
// 		debug ^= state[j];
// 	}

// 	for(j=63;j<64;j--)
// 	{
// 		printf("%lld", GET_BIT(debug,j));
// 	}
// 	printf(" <<<\n");
// }

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
		//printf("state[%lld] = %lld\n",i,state[i]);
		ret ^= state[i];
	}
	if(ret == 0)
		return 0;

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

				// print_state(state,k);
				// printf("%lld\n", ((unsigned long long)1<<(n-(k-1-i)+1)));

				if ( (state[i] ^ ((unsigned long long)1<<(n-(k-1-i)))) == 0 )
				{
					memset(state,0,k*sizeof(unsigned long long));
					//print_state(state,k);
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
	unsigned long long d = 0, i = 0, n = 16, k = 2, j = 0, Dmax = 2, cnt = 0;
	unsigned long long *state = NULL;

	state = (unsigned long long*)calloc(Dmax,sizeof(unsigned long long));

	for(k=0;k<Dmax+1;k++)
	{
		printf("k = %lld\n",k);
		init_comb(state,Dmax,k);
	    
	    for(i=0;i<100000;i++)
	    {
	    	d = next_comb(state,n,k);

	    	if(d == 0)
	    		break;

	    	cnt++;

	    	for(j=63;j<64;j--)
	    	{
	    		printf("%lld", GET_BIT(d,j));
	    	}
	    	printf("\n");
	    }
	    printf("\n");
	}

	printf("cnt = %lld\n", cnt);

    if (state)
    	free(state);

    return 0;
}