/*------------------------------------------------------------------------------
 * File: math.c
 * Created: May 15, 2015
 * Last changed: May 15, 2015
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

#include "common.h"
#include "math.h"

#include <stdlib.h>
#include <time.h>

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
 * Function: BruteforceHull()
 * Parameters:
 *   ps    Punktupps�ttningen f�r vilken ett h�lje ska genereras.
 *   hull  En pekare till h�ljet.
 *
 * Description:
 *   Genererar att konvext h�lje f�r punktupps�ttningen genom utt�mmande
 *   s�kning.
 *------------------------------------*/
void BruteforceHull(pointsetT ps, hullT *hull) {
    // For-looparna med i och j anv�nds f�r att konstruera alla t�nkbara
    // kombinationer av par bland punkterna.

    hull->numLines = 0;
    for (int i = 0; i < ps.numPoints; i++) {
        for (int j = (i+1); j < (i+ps.numPoints); j++) {
            pointT *a, *b, *c;

            // For-loopen med k anv�nds f�r att kolla att alla punkter ligger p�
            // r�tt sida av linjen mellan punkt i och j.

            bool outside = FALSE;
            for (int k = 0; k < ps.numPoints; k++) {
                // Linjen a---b �r (potentiellt) ett segment i h�ljet. Vi avg�r
                // det genom att se vilken sida om linjen som punkten c ligger
                // p�.
                a = &ps.points[i];
                b = &ps.points[j % ps.numPoints];
                c = &ps.points[k];

                // Nedan avg�r vi vilken sida om linjen punkten ligger p�.
                // d<0.0 : Punkten �r utanf�r.
                // d=0.0 : Punkten ligger p� linjen.
                // d>0.0 : Punkten �r innanf�r.
                float d = (b->x - a->x) * (c->y - a->y)
                        - (b->y - a->y) * (c->x - a->x);

                if (d < 0.0f) {
                    outside = TRUE;
                    break;
                }
            }

            // Om alla punkter visade sig vara innanf�r linjen, s� beh�ller vi
            // den som ett segment i det konvexa h�ljet.
            if (!outside) {
                hull->lines[hull->numLines].a = a;
                hull->lines[hull->numLines].b = b;
                hull->numLines++;

                if (hull->numLines >= hull->maxLines)
                    return;
            }
        }
    }

    // -
}

/*--------------------------------------
 * Function: Quickhull()
 * Parameters:
 *   ps    Punktupps�ttningen f�r vilken ett h�lje ska genereras.
 *   hull  En pekare till h�ljet.
 *
 * Description:
 *   Genererar att konvext h�lje f�r punktupps�ttningen med hj�lp av algoritmen
 *   Quickhull.
 *------------------------------------*/
void Quickhull(pointsetT ps, hullT *hull) {
    BruteforceHull(ps, hull);
    /*pointT *a = &ps.points[0],
           *b = a;

    for (int i = 1; i < ps.numPoints; i++) {
        pointT *p = &ps.points[i];

        if (p->x < a->x) a = p;
        if (p->x > b->x) b = p;
    }*/
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
        // Positionskomponenterna slumpas i intervallet [-0.5, 0.5).
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
