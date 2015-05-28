/*------------------------------------------------------------------------------
 * File: benchmark.h
 * Created: May 16, 2015
 * Last changed: May 21, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Benchmark-funktioner f�r algoritmerna.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef _benchmark_h_
#define _benchmark_h_

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

#endif // _benchmark_h_
