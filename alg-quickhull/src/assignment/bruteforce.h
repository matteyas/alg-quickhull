/*------------------------------------------------------------------------------
 * File: bruteforce.h
 * Created: May 21, 2015
 * Last changed: May 21, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Inneh�ller bruteforce-algoritmen f�r att l�sa det konvexa h�ljet.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef _bruteforce_h_
#define _bruteforce_h_

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "core/math.h"

#include "assignment/algorithmdata.h"

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: BruteforceHull()
 * Parameters:
 *   ps    Punktupps�ttningen f�r vilken ett h�lje ska genereras.
 *   hull  En pekare till h�ljet.
 *
 * Description:
 *   Genererar att konvext h�lje f�r punktupps�ttningen genom utt�mmande
 *   s�kning. Returnerar information om algoritmens arbete.
 *------------------------------------*/
algorithmDataT BruteforceHull(pointsetT ps, hullT *hull);

#endif // _bruteforce_h_
