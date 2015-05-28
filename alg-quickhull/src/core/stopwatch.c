/*------------------------------------------------------------------------------
 * File: stopwatch.c
 * Created: May 28, 2015
 * Last changed: May 28, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Win32-implementation av tidtagarur.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "stopwatch.h"

/*------------------------------------------------
 * EXTERNALS
 *----------------------------------------------*/

extern void   ResetStopwatchImpl(int id);
extern int  StopwatchElapsedImpl(int id);

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: ResetStopwatch()
 * Parameters:
 *   id  "Namnet" p� tidtagaruret. Anv�nds f�r att k�ra flera, parallella ur.
 *
 * Description:
 *   Nollst�ller det specificerade tidtagaruret.
 *------------------------------------*/
void ResetStopwatch(int id) {
    ResetStopwatchImpl(id);
}

/*--------------------------------------
 * Function: StopwatchElapsed()
 * Parameters:
 *   id  "Namnet" p� tidtagaruret. Anv�nds f�r att k�ra flera, parallella ur.
 *
 * Description:
 *   Returnerar det antal mikrosekunder som passerat sedan StopwatchReset()-
 *   funktionen anropades f�r samma tidtagarur.
 *------------------------------------*/
int StopwatchElapsed(int id) {
    return StopwatchElapsedImpl(id);
}
