/*


AUTHORS:

- Oleksandr Kazymyrov (date in ISO year-month-day format): initial version

- person (date in ISO year-month-day format): short desc

*/

/*****************************************************************************
 *       Copyright (C) 2013 Oleksandr Kazymyrov <oleksandr.kazymyrov@ii.uib.no>
 *
 *  Distributed under the terms of the GNU General Public License (GPL)
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *                  http://www.gnu.org/licenses/
 ****************************************************************************/

#include <iostream>
#include <vector>
#include <map>
#include <stdio.h>
#include <stdio.h>
#include <algorithm>
#include <sstream>
#include <time.h>

#include <m4ri/m4ri.h>

using namespace std;

struct pt
{
	unsigned long long *start;
	unsigned long long *end;
};

struct E2P_parameters
{
	unsigned long long *sbox;
	unsigned long long length;
	unsigned long long n;
	vector< mzd_t* > foundL;
	pt progressTracker;
	unsigned long long full;
	unsigned long long ncpu;
	unsigned long long cpu;
	unsigned long long debug;
	FILE *output;
};

vector< mzd_t* > is_E2P(E2P_parameters io);
