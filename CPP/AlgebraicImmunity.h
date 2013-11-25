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

#include <m4ri/m4ri.h>
#include <m4ri/io.h>
#include <stdio.h>
#include <gmp.h>
#include "../C/tools_functions.h"

struct algebraic_properties {
    unsigned long long AI;
    double SP;
    unsigned long long NE;
};

int algebraic_immunity(unsigned long long *sbox, algebraic_properties* AP, unsigned long long n, unsigned long long m, unsigned long long sparseness);