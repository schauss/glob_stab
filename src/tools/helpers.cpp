/*
    Copyright (C) 2017 Thomas Schauss

    This file is part of glob_stab.

    glob_stab is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    glob_stab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with glob_stab. If not, see <http://www.gnu.org/licenses/>.
*/

#include "helpers.h"
#include <cstdio>

#include <sys/time.h>

double sec()
{
    struct timeval tv;
    gettimeofday(&tv,0 );

        return(tv.tv_sec + 1e-6*tv.tv_usec);
}

void wait_key()
{
    printf("\n\nPush Key to continue...\n\n");
    getchar();
}
