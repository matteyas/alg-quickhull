/*------------------------------------------------------------------------------
 * File: quickhull.h
 * Created: May 21, 2015
 * Last changed: May 21, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Inneh�ller algoritmen quickhull f�r att l�sa det konvexa h�ljet.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef _quickhull_h_
#define _quickhull_h_

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "core/math.h"

#include "assignment/algorithmdata.h"

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: Quickhull()
 * Parameters:
 *   ps    Punktupps�ttningen f�r vilken ett h�lje ska genereras.
 *   hull  En pekare till h�ljet.
 *
 * Description:
 *   Genererar att konvext h�lje f�r punktupps�ttningen med hj�lp av algoritmen
 *   quickhull. Returnerar data om algoritmens arbete.
 *------------------------------------*/
algorithmDataT Quickhull(pointsetT ps, hullT *hull);

#endif // _quickhull_h_
