/*------------------------------------------------------------------------------
 * File: algorithms.c
 * Created: May 21, 2015
 * Last changed: May 21, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Innehåller algoritmerna för att lösa det konvexa höljet.
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
 *   ps    Punktuppsättningen för vilken ett hölje ska genereras.
 *   hull  En pekare till höljet.
 *
 * Description:
 *   Genererar att konvext hölje för punktuppsättningen genom uttömmande
 *   sökning. Returnerar information om algoritmens arbete.
 *------------------------------------*/
algorithmDataT BruteforceHull(pointsetT ps, hullT *hull) {
    algorithmDataT algo = { 0 };

    // For-looparna med i och j används för att konstruera alla tänkbara
    // kombinationer av par bland punkterna. Vi konstruerar dem åt båda håll,
    // d.v.s. både (a,b) och (b,a).

    hull->numLines = 0;
    for (int i = 0; i < ps.numPoints; i++) {
		pointT *a, *b, *c;
		a = &ps.points[i];
        for (int j = (i+1); j < (i+ps.numPoints); j++) {
			b = &ps.points[j % ps.numPoints];
            // For-loopen med k används för att kolla att alla punkter ligger på
            // rätt sida av linjen mellan punkt i och j.

            bool outside = FALSE;
            for (int k = 0; k < ps.numPoints; k++) {
                // Linjen a---b är (potentiellt) ett segment i höljet. Vi avgör
                // det genom att se vilken sida om linjen som punkten c ligger
                // på.                
                c = &ps.points[k];

                // Nedan avgör vi vilken sida om linjen punkten ligger på.
                //
                // d<0.0 : Punkten är utanför.
                // d=0.0 : Punkten ligger på linjen.
                // d>0.0 : Punkten är innanför.
                //
                // Dessutom har jag låtit detta utgöra algoritmens kritiska
                // operation, vilket jag hoppas är adekvat.
                algo.numOps++;
                float d = (b->x - a->x) * (c->y - a->y)
                        - (b->y - a->y) * (c->x - a->x);

                if (d < 0.0f) {
                    outside = TRUE;
                    break;
                }
            }

            // Om alla punkter visade sig vara innanför linjen, så behåller vi
            // den som ett segment i det konvexa höljet.
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

/*--------------------------------------
 * Variable: arrayPool
 *
 * Description:
 *   En kö där vi poolar arrayer för att slippa omallokeringar. Vi återanvänder
 *   helt enkelt minne istället.
 *------------------------------------*/
static queueADT arrayPool;
static queueADT intPool;

/*--------------------------------------
 * Function: GetPointArray()
 * Parameters:
 *
 * Description:
 *   Returnerar en array med pekare till punkter. Antingen skapar funktionen en
 *   ny array, eller så återanvänder den en array från poolen.
 *------------------------------------*/
static arrayADT GetPointArray() {
    if (!arrayPool)
        arrayPool = NewQueue(32);

    if (QueueIsEmpty(arrayPool))
        return NewArray(sizeof(pointT *));

    return Dequeue(arrayPool);
}
static arrayADT GetIntArray() {
	if (!intPool)
		intPool = NewQueue(32);

	if (QueueIsEmpty(intPool))
		return NewArray(sizeof(int));

	return Dequeue(intPool);
}

/*--------------------------------------
 * Function: ReleaseArray()
 * Parameters:
 *   a  Arrayen som ska släppas tillbaka till array-poolen.
 *
 * Description:
 *   Släpper tillbaka en array till poolen.
 *------------------------------------*/
static void ReleaseArray(arrayADT a) {
    if (QueueIsFull(arrayPool)) {
        FreeArray(a);
        printf("Warning: Array pool is not big enough.\n");
        return;
    }

    ResetArray(a);
    Enqueue(arrayPool, a);
}

/*--------------------------------------
 * Function: QH()
 * Parameters:
 *   hull    Array med höljets punkter.
 *   a       A-punkten i triangeln.
 *   b       B-punkten i triangeln.
 *   subset  Punkterna vi ska hantera.
 *
 * Description:
 *   Hanterar punkter och löser det konvexa höljet med hjälp av den rekursiva
 *   algoritmen quickhull.
 *------------------------------------*/
static algorithmDataT QH(arrayADT hull, pointT *a, pointT *b, arrayADT subset, pointT *r, bool clockwise) {
    algorithmDataT algo = { 0 };

    int numPoints = ArrayLength(subset);

    /*--------------------------------------------------------------------------
     * 1a. INGA PUNKTER
     *
     *   Om vi inte fick några punkter så är linjen a---b ett segment i höljet,
     *   så vi behöver inte göra något mer här.
     *------------------------------------------------------------------------*/

    if (numPoints == 0)
        return algo;

    /*--------------------------------------------------------------------------
     * 1b. EN ENDA PUNKT
     *
     *   Om vi bara fick en punkt så är segmenten a---p---b segment i höljet,
     *   där p är den enda punkten vi fått in. Vi lägger in den mellan a och b,
     *   sedan är vi klara.
     *------------------------------------------------------------------------*/

    if (numPoints == 1) {
        int n = ArrayLength(hull);
        for (int i = 0; i < n; i++) {
            pointT *point = *(pointT **)ArrayGet(hull, i);
            if (point == b) {
                // Här stoppar vi in punkten mellan a och b.
                ArrayInsert(hull, i, &r);

                algo.numOps = i;
                return algo;
            }
        }

        // Hit kommer vi aldrig.
        Fail();
    }

    /*--------------------------------------------------------------------------
     * 1c. PUNKTEN LÄNGST BORT
     *
     *   Vi bildar en triangel med a, b och p, där p är den punkt som ligger
     *   längst bort från linjen a---b. Segmenten a---p---b är garanterat delar
     *   av höljet. Punkterna inuti triangeln ignoreras. Sedan hanterar vi
     *   punkter på varsin sida om triangeln rekursivt.
     *------------------------------------------------------------------------*/

 //   float dMax  = -FLT_MAX;
 //   int   index = -1;

	//float dx = (b->x - a->x);
	//float dy = (b->y - a->y);
 //   for (int i = 0; i < numPoints; i++) {
 //       pointT *point = *(pointT **)ArrayGet(subset, i);

 //       float d = dx * (point->y - a->y)
 //               - dy * (point->x - a->x);
 //       if (d < 0.0f) d = -d;

 //       if (d > dMax) {
 //           index = i;
 //           dMax  = d;
 //       }
 //   }

 //   algo.numOps += numPoints;

 //   // Punkten farPoint är den punkt, av alla punkter i 'subset' som är längst
 //   // bort från linjen a---b.
 //   pointT *farPoint = *(pointT **)ArrayGet(subset, index);
		
    int numHullPoints = ArrayLength(hull);
    for (int i = 0; i < numHullPoints; i++) {
        pointT *point = *(pointT **)ArrayGet(hull, i);
        if (point == b) {
            // Här stoppar vi in farPoint mellan a och b.
            ArrayInsert(hull, i, &r);
            algo.numOps += i;
            break;
        }
    }

    /*--------------------------------------------------------------------------
     * 2. ANDRA PUNKTER
     *
     *   Vi bildar en triangel med a, b och p, där p är den punkt som ligger
     *   längst bort från linjen a---b. Segmenten a---p---b är garanterat delar
     *   av höljet. Punkterna inuti triangeln ignoreras. Sedan hanterar vi
     *   punkter på varsin sida om triangeln rekursivt.
     *------------------------------------------------------------------------*/

    arrayADT subsetA = GetPointArray(),
             subsetB = GetPointArray();

    algo.numAllocs += 2;

	float maxA = 0, maxB = 0;
	pointT *rA, *rB;
	pointT *comp1, *comp2;
	int ctr = 0;
	float dAx = (r->x - a->x); float dAy = (r->y - a->y); float dAgamma = dAx * a->y - dAy * a->x;
	float dBx = (b->x - r->x); float dBy = (b->y - r->y); float dBgamma = dBx * r->y - dBy * r->x;
    for (int i = 0; i < numPoints; i++) {
        pointT *point = *(pointT **)ArrayGet(subset, i);

		if (clockwise) {
			comp1 = point;
			comp2 = r;
		}
		else {
			comp2 = point;
			comp1 = r;
		}

		if (comp1->x < comp2->x) {  // hack 2: man behöver inte kolla punkter som är på fel sida så att säga
			float dA = dAy * point->x - dAx * point->y + dAgamma;
			if (dA < 0.0f) {
				ArrayAdd(subsetA, &point);
				if (dA < maxA) { // hack 1: hämta maxpunkten samtidigt som mängden byggs upp
					maxA = dA;
					rA = point;
				}
			}
		}

		if (comp1->x > comp2->x) {  // hack 2: man behöver inte kolla punkter som är på fel sida så att säga
			float dB = dBy * point->x - dBx * point->y + dBgamma;
			if (dB < 0.0f) {
				ArrayAdd(subsetB, &point);
				if (dB < maxB) { // hack 1: hämta maxpunkten samtidigt som mängden byggs upp
					maxB = dB;
					rB = point;
				}
			}
		}
    }

    algo.numOps += numPoints;

    /*--------------------------------------------------------------------------
     * 3. REKURSION
     *
     *   Här hanterar vi de två delmängderna med punkter rekursivt. Detta är
     *   algoritmens kärna. Den ena delmängden innehåller punkter på ena sidan
     *   om triangeln, och den andra innehåller punkterna på andra sidan
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

    algorithmDataT algoA = QH(hull, a, r, subsetA, rA, clockwise);
	algorithmDataT algoB = QH(hull, r, b, subsetB, rB, clockwise);

    algo.numOps    += algoA.numOps        + algoB.numOps   ;
    algo.numAllocs += algoA.numAllocs     + algoB.numAllocs;
    algo.numBytes  += algoA.numBytes      + algoB.numBytes ;

    /*--------------------------------------------------------------------------
     * 4. DEALLOKERING
     *
     *   Dags att rensa upp minnet.
     *------------------------------------------------------------------------*/

    algo.numBytes += ArrayBytes(subsetA) + ArrayBytes(subsetB);

    ReleaseArray(subsetA);
    ReleaseArray(subsetB);

    return algo;
}

/*--------------------------------------
 * Function: Quickhull()
 * Parameters:
 *   ps    Punktuppsättningen för vilken ett hölje ska genereras.
 *   hull  En pekare till höljet.
 *
 * Description:
 *   Genererar att konvext hölje för punktuppsättningen med hjälp av algoritmen
 *   quickhull. Returnerar data om algoritmens arbete.
 *------------------------------------*/
algorithmDataT Quickhull(pointsetT ps, hullT *hull) {
    algorithmDataT algo = { 0 };

    /*--------------------------------------------------------------------------
     * 1. INITIERING AV VARIABLER OCH DATA
     *
     *   Här allokerar vi minne och arrayer/vektorer som vi behöver för att
     *   lagra information om höljet och de 'subset' punkter (två "halvor" av
     *   höljet) som vi ska jobba med.
     *
     *------------------------------------------------------------------------*/

    // Vi initierar både left point och right point till den första punkten i
    // uppsättningen, bara för att inte behöva hantera null-pekare.
    pointT *leftPoint  = &ps.points[0],
           *rightPoint = leftPoint;

    // Efter att algoritmen är klar kommer hullPoints innehålla alla punkter i
    // höljer i medurs ordning.
    arrayADT hullPoints = GetPointArray(),
             subsetA    = GetPointArray(),
             subsetB    = GetPointArray();;

    // Tre allokeringar ovan.
    algo.numAllocs += 3;

    /*--------------------------------------------------------------------------
     * 2. SÖKNING EFTER EXTREMPUNKTER
     *
     *   Vi letar upp vänstra och högra extrempunkterna. Dessa är garanterat en
     *   del av det konvexa höljet, och läggs därför in i punktlistan omgående.
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
     * 3. UPPDELNING I TVÅ HALVOR
     *
     *   Vi drar en linje mellan extrempunkterna. Punkterna på de två sidorna om
     *   linjen hamnar i varsina punktsamlingar; subsetA och subsetB.
     *------------------------------------------------------------------------*/
	pointT *rA, *rB;
	float maxA = 0, maxB = 0;
	float dx = (rightPoint->x - leftPoint->x);
	float dy = (rightPoint->y - leftPoint->y);
	float dGamma = dx*leftPoint->y - dy * leftPoint->x;
    for (int i = 0; i < ps.numPoints; i++) {
        pointT *point = &ps.points[i];

        // Vänster och höger extrempunkter har vi redan hanterat.
        if (point==leftPoint || point==rightPoint)
            continue;

		float d = dy * point->x - dx * point->y + dGamma;

		if (d < 0.0f) {
			ArrayAdd(subsetA, &point); // Övre "halva."
			if (d < maxA) { // hack 1: hämta maxpunkten samtidigt som mängden byggs upp
				maxA = d;
				rA = point;
			}
		}
		else {
			ArrayAdd(subsetB, &point); // Nedra "halva."
			if (d > maxB) { // hack 1: hämta maxpunkten samtidigt som mängden byggs upp
				maxB = d;
				rB = point;
			}
		}
    }

    algo.numOps += ps.numPoints;

    /*--------------------------------------------------------------------------
     * 4. REKURSION
     *
     *   Här hanterar vi de två punktsamlingarna rekursivt, vilket är quickhull-
     *   algoritmens rekursiva kärna.
     *------------------------------------------------------------------------*/

	algorithmDataT algoA = QH(hullPoints, leftPoint, rightPoint, subsetA, rA, 1);
	algorithmDataT algoB = QH(hullPoints, rightPoint, leftPoint, subsetB, rB, 0);

    algo.numOps    += algoA.numOps    + algoB.numOps   ;
    algo.numAllocs += algoA.numAllocs + algoB.numAllocs;
    algo.numBytes  += algoA.numBytes  + algoB.numBytes ;

    /*--------------------------------------------------------------------------
     * 5. IHOPSÄTTNING AV HÖLJE
     *
     *   Vi "monterar ihop" höljet genom att länka ihop punkterna. De är sparade
     *   i rätt ordning i listan, så vi kan bara linjer från en punkt till nästa
     *   tills vi gått ett helt varv.
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
     *   Här rensar vi upp minnet efter oss.
     *------------------------------------------------------------------------*/

    algo.numBytes += ArrayBytes(hullPoints)
                  +  ArrayBytes(subsetA   )
                  +  ArrayBytes(subsetB   );

    ReleaseArray(hullPoints);
    ReleaseArray(subsetA   );
    ReleaseArray(subsetB   );

    return algo;
}

//------------------------------------------------------------------------------

#define ui unsigned int
const ui kHist = 2048;
const ui bytes_kHist = 2048 * 3 * 4;

__inline ui FloatFlip(ui f)
{
	ui mask = -(long)(f >> 31) | 0x80000000;
	return f ^ mask;
}
__inline void FloatFlipX(ui *f)
{
	ui mask = -(long)(*f >> 31) | 0x80000000;
	*f ^= mask;
}
__inline ui IFloatFlip(ui f)
{
	ui mask = ((f >> 31) - 1) | 0x80000000;
	return f ^ mask;
}
__inline void IFloatFlipX(ui *f)
{
	ui mask = ((*f >> 31) - 1) | 0x80000000;
	*f ^= mask;
}

// kollar om linjen AP ligger till höger om AB (perspektiv från A till B)
__inline char Right(pointT *A, pointT *B, pointT *P) {
	float d = (B->x - A->x) * (P->y - A->y)
		- (B->y - A->y) * (P->x - A->x);

	if (d < 0.0f)
		return 1;
	else
		return 0;
}

// kollar om linjen AP ligger till vänster om AB (perspektiv från A till B)
__inline char Left(pointT *A, pointT *B, pointT *P) {
	float d = (B->x - A->x) * (P->y - A->y)
		- (B->y - A->y) * (P->x - A->x);

	if (d > 0.0f)
		return 1;
	else
		return 0;
}

// sort it like a bitch, 4 passes tot + 2048 shits
// notera att sorted och tempsort ska vara malloc(size(int)*ps.numpoints)
#define _0(x)   (x & 0x7FF)
#define _1(x)   (x >> 11 & 0x7FF)
#define _2(x)   (x >> 22 )
algorithmDataT RaddeSort(pointsetT ps, pointT **sorted, pointT **tempsort) {
	algorithmDataT algo = { 0 };
	int ops;

	// konstruera histogram (det är typ counters på antal tal i olika storleksspann)
	
	ui *b0 = calloc(3 * kHist, 4);
	ui *b1 = b0 + kHist;
	ui *b2 = b1 + kHist;

	algo.numAllocs++;
	algo.numBytes = bytes_kHist;

	// pass 0, beräkna histograms
	for (ops = 0; ops < ps.numPoints; ops++) {
		ui fi = FloatFlip(*(ui*)&ps.points[ops].x);

		b0[_0(fi)] ++;
		b1[_1(fi)] ++;
		b2[_2(fi)] ++;
	}

	// gör nån supersummering, shrug... offsets typ
	ui sum0 = 0, sum1 = 0, sum2 = 0;
	ui tsum;
	for (int i = 0; i < kHist; i++) {
		tsum = b0[i] + sum0;
		b0[i] = sum0 - 1;
		sum0 = tsum;

		tsum = b1[i] + sum1;
		b1[i] = sum1 - 1;
		sum1 = tsum;

		tsum = b2[i] + sum2;
		b2[i] = sum2 - 1;
		sum2 = tsum;

		ops++;
	}

	// pass 1, sortera på första 11 bits
	for (int i = 0; i < ps.numPoints; i++) {
		ui *fi = (ui*)&ps.points[i].x;
		FloatFlipX(fi);
		ui pos = _0(*fi);

		sorted[++b0[pos]] = &ps.points[i];
		ops++;
	}

	// pass 2, sortera på nästa 11 bits
	for (int i = 0; i < ps.numPoints; i++) {
		ui si = *(ui*)&sorted[i]->x;
		ui pos = _1(si);

		tempsort[++b1[pos]] = sorted[i];
		ops++;
	}

	// pass 3, sortera på sista 10 bits
	for (int i = 0; i < ps.numPoints; i++) {
		ui *ai = (ui*)&tempsort[i]->x;
		ui pos = _2(*ai);

		IFloatFlipX(ai);

		sorted[++b2[pos]] = tempsort[i];
		ops++;
	}

	free(b0);

	algo.numOps = ops;
	return algo;
}

// shrug, blev snabbare att ha en funktion för det här lol
void InitHullPoints(pointT *p1, pointT *p2) {
	p1->next = p2;
	p1->prev = p2;

	p2->next = p1;
	p2->prev = p1;
}

// incremental hull
algorithmDataT Incrementhull(pointsetT ps, hullT *hull) {
	algorithmDataT algo = { 0 }, algoRadde;
	
	int msize = 4 * ps.numPoints;

	pointT **sortset = malloc(msize);
	pointT **tempsort = malloc(msize);
	
	algo.numAllocs = 2;
	algo.numBytes = msize << 1;

	algoRadde = RaddeSort(ps, sortset, tempsort);
	algo.numAllocs += algoRadde.numAllocs;
	algo.numBytes += algoRadde.numBytes;
	
	int ops = algoRadde.numOps;

	pointT
		*p1 = sortset[0],
		*p2 = sortset[1],
		*LastHP = p2, *P;
	
	InitHullPoints(p1, p2);
	// ta in nästa punkt, den är med i hull (den är längst till höger), ta reda på vilka tidigare hullpoints som ska bort
	for (int i = 2; i < ps.numPoints; i++) {
		pointT *A = sortset[i];

		// stega från senaste hullpoint åt ena hållet (det här är bilden med blå linjer)
		pointT *B = LastHP;
		pointT *P = LastHP->prev;
		while (Right(A, B, P)) {
			B = P;
			P = P->prev;
			ops++;
		}

		// stega från senaste hullpoint åt andra hållet (bilden med gröna linjer, fast B=C då hum)
		pointT *C = LastHP;
		P = LastHP->next;
		while (Left(A, C, P)) {
			C = P;
			P = P->next;
			ops++;
		}

		// uppdatera neighbours
		B->next = A;
		A->prev = B;

		C->prev = A;
		A->next = C;

		LastHP = A;
	}

	// konstruera hullfan
	P = LastHP;

	hull->numLines = 0;
	while (1) {
		hull->lines[hull->numLines].a = LastHP;
		hull->lines[hull->numLines++].b = LastHP->prev;

		LastHP = LastHP->prev;
		if (LastHP == P)
			break;
	}

	free(sortset);
	free(tempsort);

	algo.numOps = ops;

	return algo;
}