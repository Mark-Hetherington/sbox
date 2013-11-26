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
#include <algorithm>
#include <sstream>
#include <time.h>

#include <m4ri/m4ri.h>
//#include <m4ri/io.h>

using namespace std;

int test();
int is_E2P(unsigned long long *sbox, unsigned long long length, unsigned long long n);
