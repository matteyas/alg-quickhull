/*------------------------------------------------------------------------------
 * File: benchmark.c
 * Created: May 16, 2015
 * Last changed: May 16, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Benchmark-funktioner f�r quickhull.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "benchmark.h"
#include "common.h"
#include "math.h"

#include <stdio.h>
#include <time.h>

#include <Windows.h>

/*------------------------------------------------
 * CONSTANTS
 *----------------------------------------------*/

/*--------------------------------------
 * Constant: NumIterations
 *
 * Description:
 *   Antal iterationer som ska k�ras i benchmarkl�ge.
 *------------------------------------*/
#define NumIterations 5000

/*------------------------------------------------
 * GLOBALS
 *----------------------------------------------*/

LARGE_INTEGER stopwatch;

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

static benchmarkdataT InitBenchmarkData() {
    benchmarkdataT bmd = { 0 };

    bmd.minCritOps = INT_MAX;
    bmd.maxCritOps = INT_MIN;
    bmd.minTime    = INT_MAX;
    bmd.maxTime    = INT_MIN;

    return bmd;
}

static void PrintStatistics(string s, benchmarkdataT *bmd) {
    printf("STATISTICS (%s)\n"
           "------------------------------------------------------------\n"
           "                           Min.       Max.       Avg.\n"
           " Critical Operations       %-7d    %-7d    %-7d\n"
           " Memory Allocated (bytes)  0          0          0\n"
           " Execution Time (usecs)    %-5d      %-5d      %-5d\n"
           "------------------------------------------------------------\n"
           "\n",
           s,
           bmd->minCritOps, bmd->maxCritOps, (int)bmd->avgCritOps,
           bmd->minTime   , bmd->maxTime   , (int)bmd->avgTime);
}

static void StopwatchStart() {
    QueryPerformanceCounter(&stopwatch);
}

static int StopwatchStop() {
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);

    LARGE_INTEGER pc;
    QueryPerformanceCounter(&pc);

    pc.QuadPart -= stopwatch.QuadPart;
    pc.QuadPart *= 1000000;
    pc.QuadPart /= freq.QuadPart;

    return (int)pc.QuadPart;
}

/*--------------------------------------
 * Function: RunBenchmark()
 * Parameters:
 *   numPoints  Antal punkter att anv�nda vid utr�kning av det konvexa
 *              h�ljet.
 *
 * Description:
 *   K�r programmet i benchmark-l�ge.
 *------------------------------------*/
void RunBenchmark(int numPoints) {
    printf("Creating random points...");
    pointsetT ps = CreatePoints(numPoints);
    printf(" done.\n");

    printf("Initializing hull...");
    hullT hull = InitHull(ps);
    printf(" done.\n");

    printf("Benchmark will now run.\n");
    printf("Benchmarking...\n");

    benchmarkdataT bmdbf = InitBenchmarkData(),
                   bmdqh = InitBenchmarkData();

    clock_t start = clock();
    for (int i = 1; i <= NumIterations; i++) {
        RandomizePoints(ps);

        /*--------------------------------------------------------------------*/
        /* Bruteforce                                                         */
        /*--------------------------------------------------------------------*/

        StopwatchStart();
        int numCritOps = BruteforceHull(ps, &hull);
        int microSecs  = StopwatchStop();

        if (numCritOps < bmdbf.minCritOps) bmdbf.minCritOps = numCritOps;
        if (numCritOps > bmdbf.maxCritOps) bmdbf.maxCritOps = numCritOps;
        if (microSecs  < bmdbf.minTime   ) bmdbf.minTime    = microSecs;
        if (microSecs  > bmdbf.maxTime   ) bmdbf.maxTime    = microSecs;

        bmdbf.avgCritOps += (float)numCritOps / NumIterations;
        bmdbf.avgTime    += (float)microSecs  / NumIterations;

        /*--------------------------------------------------------------------*/
        /* Quickhull                                                          */
        /*--------------------------------------------------------------------*/

        StopwatchStart();
        numCritOps = Quickhull(ps, &hull);
        microSecs  = StopwatchStop();

        if (numCritOps < bmdqh.minCritOps) bmdqh.minCritOps = numCritOps;
        if (numCritOps > bmdqh.maxCritOps) bmdqh.maxCritOps = numCritOps;
        if (microSecs  < bmdqh.minTime   ) bmdqh.minTime    = microSecs;
        if (microSecs  > bmdqh.maxTime   ) bmdqh.maxTime    = microSecs;

        bmdqh.avgCritOps += (float)numCritOps / NumIterations;
        bmdqh.avgTime    += (float)microSecs  / NumIterations;

        // H�r ser vi till att skriva ut hur l�ngt i benchmarket vi kommit,
        // procentuellt sett, en g�ng varje sekund. S� att anv�ndaren inte tror
        // att programmet h�ngt sig.
        int benchmarkTime = 1000 * (clock() - start) / CLOCKS_PER_SEC;
        if (benchmarkTime >= 1000) {
            printf("%2.1f%%...\n", 100.0f * (float)i / NumIterations);
            start = clock();
        }
    }

    printf("100.0%%. Done!\n\n");

    PrintStatistics("Bruteforce", &bmdbf);
    PrintStatistics("Quickhull", &bmdqh);

    FreeHull(hull);
    FreePoints(ps);

    printf("\nPress ENTER to exit...");
    getchar();
}
