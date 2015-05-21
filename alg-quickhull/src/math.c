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

#include "array.h"
#include "benchmark.h"
#include "common.h"
#include "debug.h"
#include "math.h"

#include <float.h>
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
 *   s�kning. Returnerar det totala antalet kritiska operationer som utf�rdes.
 *------------------------------------*/
int BruteforceHull(pointsetT ps, hullT *hull) {
    int numCritOps = 0;

    // For-looparna med i och j anv�nds f�r att konstruera alla t�nkbara
    // kombinationer av par bland punkterna. Vi konstruerar dem �t b�da h�ll,
    // d.v.s. b�de (a,b) och (b,a).

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
                //
                // d<0.0 : Punkten �r utanf�r.
                // d=0.0 : Punkten ligger p� linjen.
                // d>0.0 : Punkten �r innanf�r.
                //
                // Dessutom har jag l�tit detta utg�ra algoritmens kritiska
                // operation, vilket jag hoppas �r adekvat.
                numCritOps++;
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
                if (hull->numLines >= hull->maxLines) {
                    // Segmentet f�r inte plats i h�ljets datastruktur. Om det
                    // h�r sker s� har vi antingen allokerat ett f�r litet h�lje
                    // eller, om maxLines �r detsamma som ps.numPoints, s� vet vi
                    // att h�ljet egentligen borde f� plats. Det inneb�r att vi
                    // h�ller p� att l�gga in en dubblett. Det inneb�r att vi
                    // �nd� �r klara, s� vi kan returnera h�r.
                    return numCritOps;
                }

                hull->lines[hull->numLines].a = a;
                hull->lines[hull->numLines].b = b;
                hull->numLines++;
                break;
            }
        }
    }

    return numCritOps;
}

int Quickhull2(arrayT* hull, pointT *a, pointT *b, arrayT *subset) {

    int n = ArrayLength(subset);

    if (n == 0)
        return 0;

    if (n == 1) {
        ArrayAdd(hull, ArrayGet(subset, 0));
        return 1;
        for (int i = 0; i < ArrayLength(hull); i++) {
            if (*(pointT **)ArrayGet(hull, i) == a) {
                ArrayInsert(hull, i, ArrayGet(subset, 0));
                return i;
            }
        }

        Fail();
    }

    int numOps = 0;

    int   index = 0;
    float max   = -FLT_MAX;

    for (int i = 0; i < n; i++) {
        pointT *p = *(pointT **)ArrayGet(subset, i);

        float d = (b->x - a->x) * (p->y - a->y)
                - (b->y - a->y) * (p->x - a->x);

        if (d < 0.0f)
            d = -d;

        if (d > max) {
            index = i;
            max   = d;
        }

        numOps++;
    }

    pointT *p = *(pointT **)ArrayGet(subset, index);

    arrayT subsetA, subsetB;
    InitArray(&subsetA, sizeof(pointT *));
    InitArray(&subsetB, sizeof(pointT *));

    for (int i = 0; i < n; i++) {
        pointT *q = *(pointT **)ArrayGet(subset, i);

        float d1 = (p->x - a->x) * (q->y - a->y)
                 - (p->y - a->y) * (q->x - a->x);

        if (d1 > 0.0f) {
            ArrayAdd(&subsetA, &q);
        }

        float d2 = (b->x - p->x) * (q->y - p->y)
                 - (b->y - p->y) * (q->x - p->x);

        if (d2 > 0.0f) {
            ArrayAdd(&subsetA, &q);
        }

        numOps++;
    }

    numOps += Quickhull2(hull, a, p, &subsetA)
           +  Quickhull2(hull, p, b, &subsetB);

    FreeArray(&subsetA);
    FreeArray(&subsetB);

    return numOps;
}

/*--------------------------------------
 * Function: Quickhull()
 * Parameters:
 *   ps    Punktupps�ttningen f�r vilken ett h�lje ska genereras.
 *   hull  En pekare till h�ljet.
 *
 * Description:
 *   Genererar att konvext h�lje f�r punktupps�ttningen med hj�lp av algoritmen
 *   quickhull. Returnerar det totala antalet kritiska operationer som utf�rdes.
 *------------------------------------*/
int Quickhull(pointsetT ps, hullT *hull) {
    /*--------------------------------------------------------------------------
     * 1. INITIERING AV VARIABLER OCH DATA
     *
     *   H�r allokerar vi minne och arrayer/vektorer som vi beh�ver f�r att
     *   lagra information om h�ljet och de 'subset' punkter (tv� "halvor" av
     *   h�ljet) som vi ska jobba med.
     *
     *------------------------------------------------------------------------*/

    // Vi initierar b�de lp (left point) och rp (right point) till den f�rsta
    // punkten i upps�ttningen, bara f�r att inte beh�va hantera null-pekare.
    pointT *lp = &ps.points[0],
           *rp = lp;

    // Efter att algoritmen �r klar kommer hullPoints inneh�lla alla punkter i
    // h�ljer i medurs ordning.
    arrayT hullPoints, subsetA, subsetB;
    InitArray(&hullPoints, sizeof(pointT *));
    InitArray(&subsetA   , sizeof(pointT *));
    InitArray(&subsetB   , sizeof(pointT *));

    /*--------------------------------------------------------------------------
     * 2. S�KNING EFTER EXTREMPUNKTER
     *
     *   Vi letar upp v�nstra och h�gra extrempunkterna. Dessa �r garanterat en
     *   del av det konvexa h�ljet, och l�ggs d�rf�r in i punktlistan omg�ende.
     *
     *------------------------------------------------------------------------*/

    for (int i = 1; i < ps.numPoints; i++) {
        pointT *p = &ps.points[i];

        if (p->x < lp->x) lp = p;
        if (p->x > rp->x) rp = p;
    }

    ArrayAdd(&hullPoints, &lp);
    ArrayAdd(&hullPoints, &rp);

    /*--------------------------------------------------------------------------
     * 3. UPPDELNING I TV� HALVOR
     *
     *   Vi drar en linje mellan extrempunkterna. Punkterna p� de tv� sidorna om
     *   linjen hamnar i varsina punktsamlingar; subsetA och subsetB.
     *
     *------------------------------------------------------------------------*/

    for (int i = 0; i < ps.numPoints; i++) {
        pointT *p = &ps.points[i];

        float d = (rp->x - lp->x) * (p->y - lp->y)
                - (rp->y - lp->y) * (p->x - lp->x);

        if (d < 0.0f) ArrayAdd(&subsetA, &p); // �vre "halva."
        else          ArrayAdd(&subsetB, &p); // Nedra "halva."
    }

    /*--------------------------------------------------------------------------
     * 4. REKURSION
     *
     *   H�r hanterar vi de tv� punktsamlingarna rekursivt, vilket �r quickhull-
     *   algoritmens k�rna.
     *
     *------------------------------------------------------------------------*/

    int numOps = ps.numPoints;
               + Quickhull2(&hullPoints, lp, rp, &subsetA)
               + Quickhull2(&hullPoints, rp, lp, &subsetB);

    /*--------------------------------------------------------------------------
     * 5. AVALLOKERING
     *
     *   H�r rensar vi upp minnet efter oss.
     *
     *------------------------------------------------------------------------*/

    hull->numLines = 0;
    for (int i = 0; i < ArrayLength(&hullPoints); i++) {
        if (hull->numLines >= hull->maxLines)
            break;

        hull->lines[i].a = *(pointT **)ArrayGet(&hullPoints, i);
        hull->lines[i].b = *(pointT **)ArrayGet(&hullPoints, (i+1) % ArrayLength(&hullPoints));

        hull->numLines++;
    }

    FreeArray(&hullPoints);
    FreeArray(&subsetA);
    FreeArray(&subsetB);

    return 0;
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
