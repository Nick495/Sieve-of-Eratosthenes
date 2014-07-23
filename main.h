//
//  main.h
//  Sieve of Eratosthenes
//
//  Created by Nick on 7/22/14.
//  Copyright (c) 2014 nick. All rights reserved.
//

#ifndef Sieve_of_Eratosthenes_main_h
#define Sieve_of_Eratosthenes_main_h

#include <stdio.h>   /* printf().                                             */
#include <stdlib.h>  /* malloc(), calloc(), free().                           */
#include <limits.h>  /* UINT_MAX for overflow checking.                       */
#include <math.h>    /* Sqrt(). For tighter array bounds.                     */
#include <string.h>  /* Stoll() for user input.                               */

/* pthreads for parallel solving! */
#include <unistd.h>  /* sysconf(), gives us the thread count of the machine.  */
#include <pthread.h> /* pthreads!                                             */

#endif
