/*------------------------------------------------------------------------------
 * File: stopwatch.h
 * Created: May 28, 2015
 * Last changed: May 28, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Erbjuder funktioner f�r tidsm�tning.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef _stopwatch_h_
#define _stopwatch_h_

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
void ResetStopwatch(int id);

/*--------------------------------------
 * Function: StopwatchElapsed()
 * Parameters:
 *   id  "Namnet" p� tidtagaruret. Anv�nds f�r att k�ra flera, parallella ur.
 *
 * Description:
 *   Returnerar det antal mikrosekunder som passerat sedan StopwatchReset()-
 *   funktionen anropades f�r samma tidtagarur.
 *------------------------------------*/
int StopwatchElapsed(int id);

#endif
