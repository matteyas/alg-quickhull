/*------------------------------------------------------------------------------
 * File: math.c
 * Created: May 15, 2015
 * Last changed: May 21, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Matematikbibliotek f�r l�sning av olika matematiska problem.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "array.h"
#include "math.h"

#include <math.h>

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: CreatePoints()
 * Parameters:
 *   n  Antalet punkter att skapa i upps�ttningen.
 *
 * Description:
 *   Skapar en upps�ttning punkter.
 *------------------------------------*/
pointsetT CreatePoints(int n) {
    pointsetT ps;

    ps.points    = malloc(sizeof(pointT) * n);
    ps.numPoints = n;

    return ps;
}

/*--------------------------------------
 * Function: FreePoints()
 * Parameters:
 *   ps  Den upps�ttning punkter som ska avallokeras.
 *
 * Description:
 *   Avallokerar en upps�ttning punkter,
 *------------------------------------*/
void FreePoints(pointsetT ps) {
    free(ps.points);
    ps.numPoints = 0;
}

/*--------------------------------------
 * Function: InitHull()
 * Parameters:
 *   ps  Den upps�ttning punkter som ska anv�ndas till att initiera h�ljet.
 *
 * Description:
 *   Initierar ett nytt h�lje.
 *------------------------------------*/
hullT InitHull(pointsetT ps) {
    hullT hull;

    hull.lines    = malloc(sizeof(lineT) * ps.numPoints);
    hull.numLines = 0;
    hull.maxLines = ps.numPoints;

    return hull;
}

/*--------------------------------------
 * Function: FreeHull()
 * Parameters:
 *   hull Det h�lje som ska avallokeras.
 *
 * Description:
 *   Avallokerar ett h�lje.
 *------------------------------------*/
void FreeHull(hullT hull) {
    free(hull.lines);
    hull.maxLines = 0;
    hull.numLines = 0;
}

/*--------------------------------------
 * Function: RandomizePoints()
 * Parameters:
 *   ps    Punktupps�ttningen vars positioner ska slumpas.
 *
 * Description:
 *   Slumpar positionerna f�r punkterna i den specificerade punktupps�ttningen.
 *------------------------------------*/
void RandomizePoints(pointsetT ps) {
    for (int i = 0; i < ps.numPoints; i++) {
        ps.points[i].x = (float)rand()/RAND_MAX - 0.5f;
        ps.points[i].y = (float)rand()/RAND_MAX - 0.5f;
    }
}

/*--------------------------------------
 * Function: Reflect()
 * Parameters:
 *   d  Riktningsvektorn.
 *   n  Normalvektorn som riktningsvektorn ska reflekteras mot.
 *
 * Description:
 *   Reflekterar en riktningsvektor mot en normalvektor och returnerar
 *   resultatet.
 *------------------------------------*/
pointT Reflect(pointT d, pointT n) {
    float dp = 2.0f * (d.x*n.x + d.y*n.y);
    
    pointT p = {
        d.x - dp * n.x,
        d.y - dp * n.y
    };

    return p;
}
