/*------------------------------------------------------------------------------
 * File: benchmark.h
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

#ifndef Benchmark_h
#define Benchmark_h

typedef struct {
    int   minCritOps;
    int   maxCritOps;
    float avgCritOps;

    int   minTime;
    int   maxTime;
    float avgTime;
} benchmarkdataT;

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: RunBenchmark()
 * Parameters:
 *   numPoints  Antal punkter att anv�nda vid utr�kning av det konvexa
 *              h�ljet.
 *
 * Description:
 *   K�r programmet i benchmark-l�ge.
 *------------------------------------*/
void RunBenchmark(int numPoints);

#endif // Benchmark_h
