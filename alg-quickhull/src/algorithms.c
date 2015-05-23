/*------------------------------------------------------------------------------
 * File: algorithms.c
 * Created: May 21, 2015
 * Last changed: May 21, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Inneh�ller algoritmerna f�r att l�sa det konvexa h�ljet.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "algorithms.h"
#include "array.h"
#include "benchmark.h"
#include "common.h"
#include "debug.h"
#include "math.h"
#include "queue.h"

#include <float.h>

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*****************************************************
 *  _                _        __                     *
 * | |__  _ __ _   _| |_ ___ / _| ___  _ __ ___ ___  *
 * | '_ \| '__| | | | __/ _ \ |_ / _ \| '__/ __/ _ \ *
 * | |_) | |  | |_| | ||  __/  _| (_) | | | (_|  __/ *
 * |_.__/|_|   \__,_|\__\___|_|  \___/|_|  \___\___| *
 *                                                   *
 *****************************************************/

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
algorithmdataT BruteforceHull(pointsetT ps, hullT *hull) {
    algorithmdataT algo = { 0 };

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
                algo.numOps++;
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
                if (hull->numLines >= hull->maxLines)
                    Fail();

                hull->lines[hull->numLines].a = a;
                hull->lines[hull->numLines].b = b;
                hull->numLines++;
                break;
            }
        }
    }

    return algo;
}

//------------------------------------------------------------------------------

/*********************************************
 *              _      _    _           _ _  *
 *   __ _ _   _(_) ___| | _| |__  _   _| | | *
 *  / _` | | | | |/ __| |/ / '_ \| | | | | | *
 * | (_| | |_| | | (__|   <| | | | |_| | | | *
 *  \__, |\__,_|_|\___|_|\_\_| |_|\__,_|_|_| *
 *     |_|                                   *
 *********************************************/

static queueADT arrayPool;

arrayADT GetArray() {
    if (!arrayPool)
        arrayPool = NewQueue(32);

    if (QueueIsEmpty(arrayPool))
        return NewArray(sizeof(pointT *));

    return Dequeue(arrayPool);
}

void ReleaseArray(arrayADT a) {
    ResetArray(a);

    Enqueue(arrayPool, a);
}

/*--------------------------------------
 * Function: QH()
 * Parameters:
 *   hull    Array med h�ljets punkter.
 *   a       A-punkten i triangeln.
 *   b       B-punkten i triangeln.
 *   subset  Punkterna vi ska hantera.
 *
 * Description:
 *   Hanterar punkter och l�ser det konvexa h�ljet med hj�lp av den rekursiva
 *   algoritmen quickhull.
 *------------------------------------*/
algorithmdataT QH(arrayADT hull, pointT *a, pointT *b, arrayADT subset) {
    algorithmdataT algo = { 0 };

    int numPoints = ArrayLength(subset);

    /*--------------------------------------------------------------------------
     * 1a. INGA PUNKTER
     *
     *   Om vi inte fick n�gra punkter s� �r linjen a---b ett segment i h�ljet,
     *   s� vi beh�ver inte g�ra n�got mer h�r.
     *------------------------------------------------------------------------*/

    if (numPoints == 0)
        return algo;

    /*--------------------------------------------------------------------------
     * 1b. EN ENDA PUNKT
     *
     *   Om vi bara fick en punkt s� �r segmenten a---p---b segment i h�ljet,
     *   d�r p �r den enda punkten vi f�tt in. Vi l�gger in den mellan a och b,
     *   sedan �r vi klara.
     *------------------------------------------------------------------------*/

    if (numPoints == 1) {
        int n = ArrayLength(hull);
        for (int i = 0; i < n; i++) {
            pointT *point = *(pointT **)ArrayGet(hull, i);
            if (point == b) {
                // H�r stoppar vi in punkten mellan a och b.
                ArrayInsert(hull, i, ArrayGet(subset, 0));

                algo.numOps = i;
                return algo;
            }
        }

        // Hit kommer vi aldrig.
        Fail();
    }

    /*--------------------------------------------------------------------------
     * 1c. PUNKTEN L�NGST BORT
     *
     *   Vi bildar en triangel med a, b och p, d�r p �r den punkt som ligger
     *   l�ngst bort fr�n linjen a---b. Segmenten a---p---b �r garanterat delar
     *   av h�ljet. Punkterna inuti triangeln ignoreras. Sedan hanterar vi
     *   punkter p� varsin sida om triangeln rekursivt.
     *------------------------------------------------------------------------*/

    float dMax  = -FLT_MAX;
    int   index = -1;

    for (int i = 0; i < numPoints; i++) {
        pointT *point = *(pointT **)ArrayGet(subset, i);

        float d = (b->x - a->x) * (point->y - a->y)
                - (b->y - a->y) * (point->x - a->x);
        if (d < 0.0f) d = -d;

        if (d > dMax) {
            index = i;
            dMax  = d;
        }
    }

    algo.numOps += numPoints;

    // Punkten farPoint �r den punkt, av alla punkter i 'subset' som �r l�ngst
    // bort fr�n linjen a---b.
    pointT *farPoint = *(pointT **)ArrayGet(subset, index);

    int numHullPoints = ArrayLength(hull);
    for (int i = 0; i < numHullPoints; i++) {
        pointT *point = *(pointT **)ArrayGet(hull, i);
        if (point == b) {
            // H�r stoppar vi in farPoint mellan a och b.
            ArrayInsert(hull, i, &farPoint);
            algo.numOps += i;
            break;
        }
    }

    /*--------------------------------------------------------------------------
     * 2. ANDRA PUNKTER
     *
     *   Vi bildar en triangel med a, b och p, d�r p �r den punkt som ligger
     *   l�ngst bort fr�n linjen a---b. Segmenten a---p---b �r garanterat delar
     *   av h�ljet. Punkterna inuti triangeln ignoreras. Sedan hanterar vi
     *   punkter p� varsin sida om triangeln rekursivt.
     *------------------------------------------------------------------------*/

    arrayADT subsetA = GetArray(),//NewArray(sizeof(pointT *)),
             subsetB = GetArray();//NewArray(sizeof(pointT *));

    algo.numAllocs += 2;

    for (int i = 0; i < numPoints; i++) {
        pointT *point = *(pointT **)ArrayGet(subset, i);

        float dA = (farPoint->x - a->x) * (point->y - a->y)
                 - (farPoint->y - a->y) * (point->x - a->x);
        if (dA < 0.0f) ArrayAdd(subsetA, &point);

        float dB = (b->x - farPoint->x) * (point->y - farPoint->y)
                 - (b->y - farPoint->y) * (point->x - farPoint->x);
        if (dB < 0.0f) ArrayAdd(subsetB, &point);
    }

    algo.numOps += numPoints;

    /*--------------------------------------------------------------------------
     * 3. REKURSION
     *
     *   H�r hanterar vi de tv� delm�ngderna med punkter rekursivt. Detta �r
     *   algoritmens k�rna. Den ena delm�ngden inneh�ller punkter p� ena sidan
     *   om triangeln, och den andra inneh�ller punkterna p� andra sidan
     *   triangeln:
     *
     *              farPoint
     *                /\
     *               /  \
     *   (subsetA)  /    \  (subsetB)
     *             /      \
     *            /________\
     *          a            b
     *------------------------------------------------------------------------*/

    algorithmdataT algoA = QH(hull, a       , farPoint, subsetA);
    algorithmdataT algoB = QH(hull, farPoint, b       , subsetB);

    algo.numOps    += algoA.numOps        + algoB.numOps   ;
    algo.numAllocs += algoA.numAllocs     + algoB.numAllocs;
    algo.numBytes  += algoA.numBytes      + algoB.numBytes ;

    /*--------------------------------------------------------------------------
     * 4. DEALLOKERING
     *
     *   Dags att rensa upp minnet.
     *------------------------------------------------------------------------*/

    algo.numBytes += ArrayBytes(subsetA) + ArrayBytes(subsetB);

    ReleaseArray(subsetA);//FreeArray(subsetA);
    ReleaseArray(subsetB);//FreeArray(subsetB);

    return algo;
}

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
algorithmdataT Quickhull(pointsetT ps, hullT *hull) {
    algorithmdataT algo = { 0 };

    /*--------------------------------------------------------------------------
     * 1. INITIERING AV VARIABLER OCH DATA
     *
     *   H�r allokerar vi minne och arrayer/vektorer som vi beh�ver f�r att
     *   lagra information om h�ljet och de 'subset' punkter (tv� "halvor" av
     *   h�ljet) som vi ska jobba med.
     *
     *------------------------------------------------------------------------*/

    // Vi initierar b�de left point och right point till den f�rsta punkten i
    // upps�ttningen, bara f�r att inte beh�va hantera null-pekare.
    pointT *leftPoint  = &ps.points[0],
           *rightPoint = leftPoint;

    // Efter att algoritmen �r klar kommer hullPoints inneh�lla alla punkter i
    // h�ljer i medurs ordning.
    arrayADT hullPoints = GetArray(),//NewArray(sizeof(pointT *)),
             subsetA    = GetArray(),//NewArray(sizeof(pointT *)),
             subsetB    = GetArray();//NewArray(sizeof(pointT *));

    // Tre allokeringar ovan.
    algo.numOps += 3;

    /*--------------------------------------------------------------------------
     * 2. S�KNING EFTER EXTREMPUNKTER
     *
     *   Vi letar upp v�nstra och h�gra extrempunkterna. Dessa �r garanterat en
     *   del av det konvexa h�ljet, och l�ggs d�rf�r in i punktlistan omg�ende.
     *------------------------------------------------------------------------*/

    for (int i = 1; i < ps.numPoints; i++) {
        pointT *point = &ps.points[i];

        if (point->x < leftPoint ->x) leftPoint  = point;
        if (point->x > rightPoint->x) rightPoint = point;
    }

    ArrayAdd(hullPoints, &leftPoint );
    ArrayAdd(hullPoints, &rightPoint);

    algo.numOps += ps.numPoints - 1;

    /*--------------------------------------------------------------------------
     * 3. UPPDELNING I TV� HALVOR
     *
     *   Vi drar en linje mellan extrempunkterna. Punkterna p� de tv� sidorna om
     *   linjen hamnar i varsina punktsamlingar; subsetA och subsetB.
     *------------------------------------------------------------------------*/

    for (int i = 0; i < ps.numPoints; i++) {
        pointT *point = &ps.points[i];

        // V�nster och h�ger extrempunkter har vi redan hanterat.
        if (point==leftPoint || point==rightPoint)
            continue;

        float d = (rightPoint->x - leftPoint->x) * (point->y - leftPoint->y)
                - (rightPoint->y - leftPoint->y) * (point->x - leftPoint->x);

        if (d < 0.0f) ArrayAdd(subsetA, &point); // �vre "halva."
        else          ArrayAdd(subsetB, &point); // Nedra "halva."
    }

    algo.numOps += ps.numPoints;

    /*--------------------------------------------------------------------------
     * 4. REKURSION
     *
     *   H�r hanterar vi de tv� punktsamlingarna rekursivt, vilket �r quickhull-
     *   algoritmens rekursiva k�rna.
     *------------------------------------------------------------------------*/

    algorithmdataT algoA = QH(hullPoints, leftPoint , rightPoint, subsetA);
    algorithmdataT algoB = QH(hullPoints, rightPoint, leftPoint , subsetB);

    algo.numOps    += algoA.numOps    + algoB.numOps   ;
    algo.numAllocs += algoA.numAllocs + algoB.numAllocs;
    algo.numBytes  += algoA.numBytes  + algoB.numBytes ;

    /*--------------------------------------------------------------------------
     * 5. IHOPS�TTNING AV H�LJE
     *
     *   Vi "monterar ihop" h�ljet genom att l�nka ihop punkterna. De �r sparade
     *   i r�tt ordning i listan, s� vi kan bara linjer fr�n en punkt till n�sta
     *   tills vi g�tt ett helt varv.
     *------------------------------------------------------------------------*/

    hull->numLines = 0;
    int numHullPoints = ArrayLength(hullPoints);
    for (int i = 0; i < numHullPoints; i++) {
        if (hull->numLines >= hull->maxLines)
            Fail();

        hull->lines[i].a = *(pointT **)ArrayGet(hullPoints, i);
        hull->lines[i].b = *(pointT **)ArrayGet(hullPoints, (i+1) % numHullPoints);
        hull->numLines++;
    }

    algo.numOps += numHullPoints;

    /*--------------------------------------------------------------------------
     * 6. DEALLOKERING
     *
     *   H�r rensar vi upp minnet efter oss.
     *------------------------------------------------------------------------*/

    algo.numBytes += ArrayBytes(hullPoints)
                  +  ArrayBytes(subsetA   )
                  +  ArrayBytes(subsetB   );

    ReleaseArray(hullPoints);//FreeArray(hullPoints);
    ReleaseArray(subsetA);//FreeArray(subsetA);
    ReleaseArray(subsetB);//FreeArray(subsetB);

    return algo;
}

//------------------------------------------------------------------------------
