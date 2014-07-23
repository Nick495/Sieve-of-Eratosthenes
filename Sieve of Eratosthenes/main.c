/*
*  main.c
*  Sieve of Eratosthenes
*
*  Created by Nick on 7/22/14.
*  Copyright (c) 2014 nick. All rights reserved.
*/

typedef unsigned long ulong;

#include "main.h"

#define INT_BIT 8*sizeof(int)
/* These defines are taken from the C-FAQ. http://c-faq.com/misc/bitsets.html */
/* 
 * Adjusted for int to operate on word size, which was suggested to 
 * improve performance. 
*/
#define BITMASK(b) (1 << ((b) % CHAR_BIT))
#define BITSLOT(b) ((b) / CHAR_BIT)
#define BITSET(a, b) ((a)[BITSLOT(b)] |= BITMASK(b))
#define BITCLEAR(a, b) ((a)[BITSLOT(b)] &= ~BITMASK(b))
#define BITTEST(a, b) ((a)[BITSLOT(b)] & BITMASK(b))
#define BITNSLOTS(nb) ((nb + CHAR_BIT - 1) / CHAR_BIT)

/* Icky global variables for our thread syncronization */
pthread_mutex_t top_sq_lock;
pthread_cond_t top_sq;

#define START 5

/* This is partly for the threads, and partly for our work allocation scheme.
 * start is the starting offset for all threads, which needs to be changed based
 * on how many threads one is running. In our case, we're using 8 threads, so
 * we need a start value of 5. (The minimum start value is 3.) Note that the
 * program will NOT find primes under this start value, and that this start 
 * value's square (i.e. 5^2 = 25) must be greater than all of the threads
 * first rows, in order to ensure correct execution. Assuming that one is 
 * storing all odd numbers, the starting value is given by the inequality
 * tcnt < start(start-1)/2 + 1 where start is the start value (an odd number)
 * and tcnt is the total number of threads on which the program is run. 
 * 
 * tid is the thread id of the thread, which is also it's starting offset in
 * the array. tcount is the total number of threads being run, which could be
 * a global variable (it tells the thread how many positions to jump over in the
 * loop). imax is the maximum number of iterations of the i loop, and jmax the 
 * max the number for the jloop. int* sieve points to the bit array, which is
 * registered as an int pointer because operating on words is most efficient.
 * max_sieved_sq could also be a global variable, but in this case is a pointer
 * to the square of the largest number that we've sieved successfully ( gone
 * through the i-loop successfully for ). This
 */
struct psoe_args {
    ulong start;
    ulong tid;
    ulong tcount;
    ulong imax;
    ulong jmax;
    char* sieve;
    ulong *max_sieved_sq;
};

/* Returns the real value for a checked index of the sieve. */
ulong itoval(ulong start, ulong i)
{
    return start + 2 * i;
}

ulong valtoi(ulong start, ulong val)
{
    return (val - start) / 2;
}

/* One instance of the parallel sieve, which will be run on numerous threads */
void *psoe(void *arg_ptr)
{
    struct psoe_args args = *(struct psoe_args*)arg_ptr;
    
    for (ulong i = args.tid; i < args.imax; i += args.tcount) {
        // Find a way to spin until max_sieved_sq is > i.
        ulong ival = itoval(args.start, i);
        ulong isq = ival * ival;
        /* We can't be certain of isq's primality until max_sieved_sq is > it.*/
        pthread_mutex_lock(&top_sq_lock);
        while (i > *(args.max_sieved_sq))
            pthread_cond_wait(&top_sq, &top_sq_lock);
        pthread_mutex_unlock(&top_sq_lock);
        
        if (BITTEST(args.sieve, i)) /* Non prime, do nothing. */
            goto up_max;
        
        for (ulong j = isq; j < args.jmax; j+=ival) {
            /*
            * It's a very important part of our algorithm that we DON'T
            * need to lock here. Due to the way we've written our algorithm
            * (with max_sieved_sq, which is global to all threads), as long
            * as the spin condition above is true, we're fine to check
            * sieve[i] (as it's in the range of variables that have already
            * been crossed off, and since and since that condition is true
            * for all threads, we guarentee that no thread will be reading
            * sieve[j] while we're writing to it. This is why this algorithm
            * is more performant than other variants.
            */
            BITSET(args.sieve, valtoi(args.start, j));
        }
        
        /* Now we can update max_sieved_sq */
        /* Make sure i*i isn't an overflow */
        
    up_max: if (isq < *(args.max_sieved_sq))/* If less than the max, ignore. */
            continue;
        
        /* Otherwise, it might be the new max, so lock and check. */
        pthread_mutex_lock(&top_sq_lock);
        if (isq > *(args.max_sieved_sq))
            *(args.max_sieved_sq) = isq;
        pthread_cond_signal(&top_sq);
        pthread_mutex_unlock(&top_sq_lock);
    }
    
    return NULL;
}

/* Wrapper around psoe, initializes data needed for parallel solving. */
char *psoe_wrapper(ulong start, ulong max, ulong tcount)
{
    /* Initialize our thread synchronization */
    pthread_mutex_init(&top_sq_lock, NULL);
    pthread_cond_init(&top_sq, NULL);
    
    /*
     * Chunk is the amount of work that each thread will do, leftovers
     * deals with remaining work that doesn't evenly divide. We use a simple
     * scheme ( which should work well for small thread counts relative to
     * the amount of work ) where we allocate leftover work to the earliest
     * threads. (So threads 1-7 can get at most 1 extra unit of work). There
     * is something to be said for using a queue like structure, but considering
     * that most numbers don't require doing any work at all (we just skip them,
     * since most numbers are composite ), it seems that that methoud would be
     * both overcomplicated, and inefficent, since syncing threads requires a
     * nontrivial amount of work, which this method avoids in large part.
    */
    ulong itop = valtoi(start, (ulong) sqrt(max) + 1);
    ulong leftovers = itop % tcount;
    ulong max_sieved_sq = start * start;
    
    char *sieve = calloc(BITNSLOTS(valtoi(start, max)), sizeof(char));
    if (!sieve) {
        perror("malloc");
        return NULL;
    }
    
    pthread_t *tpool = malloc(sizeof(pthread_t) * tcount);
    if (!tpool) {
        free(sieve);
        perror("malloc");
        return NULL;
    }
    
    struct psoe_args *argpool = malloc(sizeof(struct psoe_args) * tcount);
    if(!argpool) {
        perror("malloc");
        free(tpool);
        free(sieve);
        return NULL;
    }
    
    for (ulong i = 0; i < tcount; ++i ) {
        /* Nonchanging parameters */
        argpool[i].start = start;
        argpool[i].tid = i;
        argpool[i].tcount = tcount;
        argpool[i].jmax = max;
        argpool[i].max_sieved_sq = &max_sieved_sq;
        argpool[i].sieve = sieve;
        
        /* Parameters subject to change */
        argpool[i].imax = itop;
        
        /*
         * This is part of our work distribution scheme. It allocates the work
         * over all of our threads, and gives the leftover work ( that which
         * can't be evenly distributed ) to the earliest threads, adjusting the
         * range of data that they operate on in step.
         */
        if (i < leftovers )
            argpool[i].imax += tcount;
        
        if (pthread_create(&tpool[i], NULL, psoe, &argpool[i])){
            perror("pthread_create");
            free(argpool);
            free(tpool);
            free(sieve);
            return NULL;;
        }
    }
    
    for (ulong i = 0; i < tcount; ++i )
        pthread_join(tpool[i], NULL);
   
    pthread_mutex_destroy(&top_sq_lock);
    pthread_cond_destroy(&top_sq);
    free(argpool);
    free(tpool);
    
    return sieve;
}

void usage(const char *progname)
{
    fprintf(stderr, "Usage: %s top.\nWhere top is an integer to which you'd "
            "like to enumerate the primes.\n", progname);
    exit(EXIT_SUCCESS);
    return;
}

int main(int argc, const char * argv[])
{
#define TOP 1000000000
#define START 5
    
    if (argc != 2) {
        usage(argv[0]);
    }
    
    ulong max = atoll(argv[1]);
    
    ulong nCPU = sysconf( _SC_NPROCESSORS_ONLN );
    
    if ( max > ULONG_MAX ) { /* We'll overflow */
        printf("Sorry, too big.\n");
        return 0;
    }
    
    printf("2, 3, ");
    
    char* sieve = psoe_wrapper(START, max, nCPU);
    if (!sieve) {
        return 1;
    }
    
    for (ulong i = 0; i < valtoi(START, max); ++i) {
        if(!BITTEST(sieve, i))
            printf("%lu\n", itoval(START, i));
    }
    
    return 0;
}

