/*
	Auxiliary functions

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

#include "tools_functions.h"

long long sumArray(char *arr, unsigned long long len)
{
	unsigned long long i=0;
	unsigned long long sum=0;

	for(i=0;i<len;i++)
		sum+=arr[i];

	return sum;
}
