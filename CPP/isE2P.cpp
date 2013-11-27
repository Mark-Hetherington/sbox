//============================================================================
// Name        : ThesisAlgorithm.cpp
// Author      : Maxim Storetvedt, Oleksandr Kazymyrov
// Description : An implementation of Browning K. A. et al. algorithm.
//============================================================================
#include "isE2P.h"

void										clearColumn(mzd_t *, unsigned long long column);
vector< mzd_t* >							findMatrix(map<unsigned long long, vector< mzd_t* > >, E2P_parameters);
map<unsigned long long, vector< mzd_t* > >	findSigmas(const E2P_parameters);
bool 										matrixCheck(const mzd_t*, unsigned long long, map<unsigned long long, vector< mzd_t* > >);
void										print_matrix(const mzd_t*, string);
void										print_Sigmas(map<unsigned long long, vector< mzd_t* > >);
mzd_t* 										tryCombine(const mzd_t*, vector< mzd_t* >);
template<typename VectorOrMap> int 			unique(VectorOrMap list, mzd_t* el);
void										updatedMat(unsigned long long *, unsigned long long, mzd_t *);

/*
 * The main function. Interface between Cython and C++.
 */
vector< mzd_t* > is_E2P(E2P_parameters io)
{
	vector< mzd_t* > M;
	map<unsigned long long, vector< mzd_t* > > Sigmas;
	unsigned long long i = 0, j = 0;

	if(io.sbox[0] != 0)
	{
		for(i=1;i<io.length;i++)
			io.sbox[i] = io.sbox[0] ^ io.sbox[i];
		io.sbox[0] = 0;
	}

	// printf("length = %lld\n",io.length);
	printf("findSigmas: \n");

	Sigmas = findSigmas(io);

	// for (i = 0; i < Sigmas.size(); i++)
	// {
	// 	printf("Sigmas[%lld] = %lld\n",i,Sigmas[i].size());
	// }

	// Print sbox
	// printf("sbox:\n");
	// for(i = 0; i < length; i++)
	// {
	// 	printf("%02X ", sbox[i]);
	// 	if((i+1)%8 == 0)
	// 		printf("\n");
	// }
	// printf("\n");

	// Print sigmas
	//print_Sigmas(Sigmas);
	// printf("Sigmas:\n");
	// for(i = 0; i < Sigmas.size(); i++)
	// {
	// 	printf("Sigmas[%lld]: %lld\n",i,Sigmas[i].size());
	// }
	// printf("\n");	

	printf("findMatrix: \n");

	M = findMatrix(Sigmas,io);

	// Print the matrix
	// if(!mzd_is_zero(M))
	// {
	// 	print_matrix(M,"M");
	// }

	for(i=0;i<Sigmas.size();i++)
	{
		for(j=0;j<Sigmas[i].size();j++)
		{
			mzd_free(Sigmas[i][j]);
		}
	}

	Sigmas.clear();

	return M;
}

/*
 * Clears the column of the matrix
 */
void clearColumn(mzd_t *L, unsigned long long column)
{
	unsigned long long i = 0;

	for(i=0; i < L->nrows; i++)
		mzd_write_bit(L,i,column,0);
}

vector< mzd_t* > findMatrix(map<unsigned long long, vector< mzd_t* > > Sigmas, E2P_parameters io)
{
	// Variables for testing performance
	//clock_t time_updatedMat[3] = {0,0,0}, time_gauss[3] = {0,0,0}, time_matrixCheck[3] = {0,0,0}, time_find[3] = {0,0,0}, time_tryCombine[3] = {0,0,0};

	//  The working matrices
	mzd_t *L, *M, *T;

	// The found matrices will be stored here
	vector< mzd_t* > foundM;

	// Additional variables
	unsigned long long i = 0, rank = 0, column = 0;

	// The progress tracker and the tracker of max values 
	unsigned long long *maxValueTracker = NULL;

	// Initialises the empty matrices
	T = mzd_init(io.n, io.n<<1);
	L = mzd_init(io.n, io.n<<1);
	M = mzd_init(io.n<<1, io.n<<1);

	// Defines the highest achieveable value for a column. This is 2^(n-1) for most
	// except the identity matrix part.
	maxValueTracker = (unsigned long long*)calloc(io.n<<1,sizeof(unsigned long long));

	for(i = 0; i < io.n<<1; i++)
	{
		if(i < io.n)
		{
			maxValueTracker[i] = (1 << (i+1)) - 1;
		}
		else
		{
			maxValueTracker[i] = (1 << io.n) - 1;
		}
	}

	// Several prints for debugging
	// for(i=0;i<foundL.size();i++)
	// 	print_matrix(foundL[i],"L");

	// printf("progressTracker.start\t: [");
	// for(i = 0; i < io.n<<1; i++)
	// {
	// 	if (i != ((io.n<<1)-1) )
	// 		printf("%lld,",io.progressTracker.start[i]);
	// 	else
	// 		printf("%lld]\n",io.progressTracker.start[i]);
	// }

	// printf("progressTracker.end\t: [");
	// for(i = 0; i < io.n<<1; i++)
	// {
	// 	if (i != ((io.n<<1)-1) )
	// 		printf("%lld,",io.progressTracker.end[i]);
	// 	else
	// 		printf("%lld]\n",io.progressTracker.end[i]);
	// }

	// for(i = 0; i < n<<1; i++)
	// 	printf("maxValueTracker[%lld] = %lld\n",i,maxValueTracker[i]);

	while(true)
	{
		//if(io.progressTracker.start[(io.n<<1)-(io.n>>1)] == 0)
		if(io.progressTracker.start[15] == 0)
		{
			printf("progressTracker\t\t: [");
			for(i = 0; i < io.n<<1; i++)
			{
				if (i != ((io.n<<1)-1) )
					printf("%lld,",io.progressTracker.start[i]);
				else
					printf("%lld] (%ld) // (%lld,%d)\n", io.progressTracker.start[i],io.foundL.size(),column,maxValueTracker[column]);
			}

			// printf("time_updatedMat\t\t: %f\n", (double)(time_updatedMat[2]) / CLOCKS_PER_SEC);
			// printf("time_gauss\t\t: %f\n", (double)(time_gauss[2]) / CLOCKS_PER_SEC);
			// printf("time_matrixCheck\t: %f\n", (double)(time_matrixCheck[2]) / CLOCKS_PER_SEC);
			// printf("time_find\t\t: %f\n", (double)(time_find[2]) / CLOCKS_PER_SEC);
			// printf("time_tryCombine\t\t: %f\n", (double)(time_tryCombine[2]) / CLOCKS_PER_SEC);
			// printf("~~~~~~~~~~~~~~~~~~~~~\n");
		}

		if (io.progressTracker.start[column] > maxValueTracker[column])
		{
			if( memcmp(io.progressTracker.start,io.progressTracker.end,(io.n<<1)*sizeof(unsigned long long)) >= 0 )
			{
				// printf("progressTracker\t\t: [");
				// for(i = 0; i < io.n<<1; i++)
				// {
				// 	if (i != ((io.n<<1)-1) )
				// 		printf("%lld,",io.progressTracker.start[i]);
				// 	else
				// 		printf("%lld]\n",io.progressTracker.start[i]);
				// }
				break;
			}
			io.progressTracker.start[column] = 0;
			clearColumn(L, column);
			column--;
			if( column == (unsigned long long)(-1) )
			{
				break;
			}
			io.progressTracker.start[column]++;
			continue;
		}
	
		//time_updatedMat[0] = clock();
		updatedMat(io.progressTracker.start, column, L);
		//time_updatedMat[1] = clock();
		//time_updatedMat[2] += (time_updatedMat[1] - time_updatedMat[0]);

		// print_matrix(L,"L");

		//time_gauss[0] = clock();
		mzd_copy(T,L);
		rank = mzd_echelonize_pluq(L,1);
		//time_gauss[1] = clock();
		//time_gauss[2] += (time_gauss[1] - time_gauss[0]);

		// print_matrix(L,"L");
		// print_matrix(T,"T");
		// printf("~~~~~~~~~~~~~~~~~~~~~\n");

		// printf("(rank,column) = (%lld,%lld)\n", rank, column);

		if( (!mzd_equal(T,L)) or ( ( column == (2*io.n-1) )  and (rank != io.n)  ) ) 
		{
			io.progressTracker.start[column]++;
			continue;
		}

		//time_matrixCheck[0] = clock();
		if(matrixCheck(L, column, Sigmas))
		{
			//time_matrixCheck[1] = clock();
			//time_matrixCheck[2] += (time_matrixCheck[1] - time_matrixCheck[0]);
			if( column == (2*io.n-1) )
			{
				//time_find[0] = clock();
				// Is the following if redundant? It make sence when foundL is predefined.
				// It is nessessary to check this.
				//if(find(foundL.begin(), foundL.end(), L) == foundL.end())
				if(unique(io.foundL,L))
				{
					//time_find[1] = clock();
					//time_find[2] += (time_find[1] - time_find[0]);
					// printf("progressTracker = [");
					// for(i = 0; i < 2*io.n; i++)
					// {
					// 	if(i != (2*io.n - 1) )
					// 		printf("%d,", io.progressTracker.start[i]);
					// 	else
					// 		printf("%d] (%ld)\n", io.progressTracker.start[i],io.foundL.size() + 1);
					// }
					// printf(">> %d\n",__LINE__);
					//time_tryCombine[0] = clock();
					M = tryCombine(L,io.foundL);
					//time_tryCombine[1] = clock();
					//time_tryCombine[2] += (time_tryCombine[1] - time_tryCombine[0]);

					if (mzd_is_zero(M))
					{
						// printf(">> %d\n",__LINE__);
						io.foundL.push_back(mzd_copy(NULL,L));
					}
					else
					{
						foundM.push_back(mzd_copy(NULL,M));
						// printf(">> %d\n",__LINE__);
						if(!io.full)
						{
							if(maxValueTracker)
								free(maxValueTracker);

							mzd_free(T);
							mzd_free(L);
							mzd_free(M);

							return foundM;
						}
					}
				}
				io.progressTracker.start[column]++;
			}
			else
			{
				// printf(">> %d\n",__LINE__);
				column++;
			}
		}
		else
		{
			//time_matrixCheck[1] = clock();
			//time_matrixCheck[2] += (time_matrixCheck[1] - time_matrixCheck[0]);
			io.progressTracker.start[column]++;
		}
	}

	// printf("Done. Number of linear functions is %ld\n", foundL.size());
	// printf("Number of matriceis is %ld\n", foundM.size());

	if(maxValueTracker)
		free(maxValueTracker);

	mzd_free(T);
	mzd_free(L);
	mzd_free(M);

 	return foundM;
}

/*
 * Returns the map of Sigmas
 */
map<unsigned long long, vector< mzd_t* > > findSigmas(const E2P_parameters io)
{
	map<unsigned long long, vector< mzd_t* > > Sigmas;

	unsigned long long xy = 0, Fxy = 0, i = 0, j = 0, k = 0, nbits = 0;
	mzd_t* sigma;

	sigma = mzd_init(io.n<<1, 1);

	for (i = 0; i < (unsigned long long)(1<<io.n); i++)
	{
		for (j = i+1; j < (unsigned long long)(1<<io.n); j++)
		{
			//printf("j = %lld\n",j);
			nbits = 0;
			xy = i^j;
			Fxy = (io.sbox[i]^io.sbox[j]);

			for(k = 0; k < io.n; k++)
			{
				mzd_write_bit(sigma,k,0,(xy >> k) & 1);
				mzd_write_bit(sigma,io.n+k,0,(Fxy >> k) & 1);
			}

			for(k = 2*io.n - 1; k < 2*io.n; k--)
			{
				if (mzd_read_bit(sigma,k,0) == 1)
				{
					nbits = k;
					break;
				}
			}

			// if (nbits == 2)
			// {
			// 	printf("================\n");
			// 	printf(" xy = %02llX\n", xy);
			// 	printf("Fxy = %02llX\n", Fxy);
			// 	mzd_print(sigma);
			// 	printf("nbits = %lld\n", nbits);
			// 	printf("================\n");
			// }

			// Add only unique vectors
			Sigmas[nbits].push_back(mzd_copy(NULL,sigma));
		}
	}

	// for (i = 0; i < Sigmas.size(); i++)
	// {
	// 	printf("Sigmas[%lld] = %lld\n",i,Sigmas[i].size());
	// }

	// for (i = 0; i < Sigmas.size(); i++)
	// {
	// 	printf("i = %lld\n",i);
	// 	for (j = 0; j < Sigmas[i].size(); j++)
	// 	{
	// 		k = j+1;
	// 		while(k < Sigmas[i].size())
	// 		{
	// 			if(mzd_equal(Sigmas[i][j],Sigmas[i][k]))
	// 			{
	// 				mzd_free(Sigmas[i][k]);
	// 				Sigmas[i].erase(Sigmas[i].begin() + k);
	// 				k--;
	// 			}
	// 			k++;
	// 		}
	// 	}
	// }

	// for (i = 0; i < Sigmas.size(); i++)
	// {
	// 	printf("Sigmas[%lld] = %lld\n",i,Sigmas[i].size());
	// }

	mzd_free(sigma);

	return Sigmas;
}

/**
 * Iterates through the columns of a matrix, and checks if they all are accepted with their corresponding sigmas
 */
bool matrixCheck(const mzd_t *L, unsigned long long column, map<unsigned long long, vector< mzd_t* > > Sigmas)
{
	unsigned long long s = 0;
	mzd_t *T = NULL;

	T = mzd_init(L->nrows,1);

	for(s=0;s<Sigmas[column].size();s++)
	{
		// printf("T:\t\t");
		// mzd_info(T,0);
		// printf("L:\t\t");
		// mzd_info(L,0);
		// printf("Sigma[%lld][%lld]:\t",column,s);
		// mzd_info(Sigmas[column][s],0);
		mzd_mul(T, L, Sigmas[column][s], 0);

		if (mzd_is_zero(T) != 0)
		{
			mzd_free(T);
			return false;
		}
	}

	mzd_free(T);

	return true;
}

// Print the given matrix.
void print_matrix(const mzd_t *L, string str)
{
	unsigned long long i = 0, j = 0;

	cout << str << ":" << endl;

	for (i = 0; i < L->nrows; i++)
	{
		for(j = 0; j < L->ncols; j++)
		{
			cout << mzd_read_bit(L,i,j) << " ";
		}
		cout << endl;
	}
}

// Print Sigmas
void print_Sigmas(map<unsigned long long, vector< mzd_t* > > Sigmas)
{
	unsigned long long i = 0, j = 0, s = 0, m = 0;

	cout << "Sigmas:" << endl;

	for (s = 0; s < Sigmas.size(); s++)
	{
		cout << "Sigmas[" << s << "] (" << Sigmas[s].size() << ")" << " = [";
		for (m = 0; m < Sigmas[s].size(); m++)
		{
			cout << "[";
			for (i = 0; i < Sigmas[s][m]->nrows; i++)
			{
				for(j = 0; j < Sigmas[s][m]->ncols; j++)
				{
					cout << mzd_read_bit(Sigmas[s][m],i,j) << " ";
				}
			}
			if (m != Sigmas[s].size()-1)
				cout << "],";
			else
				cout << "]";
		}
		cout << "]" << endl;
	}
}

/*
 * Try to combines two n x 2n matrices into one invertible 2n x 2n matrix.
 */
mzd_t* tryCombine(const mzd_t *L, vector< mzd_t* > foundL)
{
	unsigned long long rank = 0, k = 0, n = L->nrows;
	mzd_t *M = NULL, *T = NULL;

	M = mzd_init(n<<1,n<<1);
	T = mzd_init(n<<1,n<<1);

	for(k = 0; k < foundL.size(); k++)
	{
		mzd_stack(M,foundL[k],L);
		mzd_copy(T,M);

		rank = mzd_echelonize_pluq(T,1);

		if( rank == (n<<1) )
		{
			mzd_free(T);
			return M;
		}
	}

	mzd_free(T);
	mzd_free(M);

	return mzd_init(n<<1,n<<1);
}

/*
 * Find the element in the list. Return 1 if found and 0 otherwise.
 */
template <typename VectorOrMap>
int unique(VectorOrMap list, mzd_t* el)
{
    for(typename VectorOrMap::iterator it = list.begin(); it != list.end(); ++it )
    	if(mzd_equal(*it,el))
			return 0;

	return 1;
}

/*
 * Updates the column of the matrix
 */
void updatedMat(unsigned long long *progressTracker, unsigned long long column, mzd_t *L)
{
	unsigned long long i = 0;

	for(i=0; i < L->nrows; i++)
		mzd_write_bit(L,i,column,(progressTracker[column] >> i) & 1);
}
