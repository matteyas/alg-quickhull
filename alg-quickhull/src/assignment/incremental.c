/*------------------------------------------------------------------------------
 * File: quickhull.c
 * Created: May 28, 2015
 * Last changed: May 28, 2015
 *
 * Author(s): "Philip Arvidsson" (philip@philiparvidsson.com)
 *
 * Description:
 *   Innehåller algoritmen quickhull för att lösa det konvexa höljet.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

/*------------------------------------------------
 * DEFINES
 *----------------------------------------------*/

/*--------------------------------------
 * Define: ArrayPooling
 *
 * Description:
 *   Kommentera ut raden nedan för att inte använda array-poolen.
 *------------------------------------*/
#define ArrayPooling

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "incremental.h"

#include "core/common.h"
#include "core/debug.h"
#include "core/math.h"
#include "core/collections/array.h"
#include "core/collections/queue.h"

#include "assignment/algorithmdata.h"

#include <float.h>

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Variable: arrayPool
 *
 * Description:
 *   En kö där vi poolar arrayer för att slippa omallokeringar. Vi återanvänder
 *   helt enkelt minne istället.
 *------------------------------------*/
static queueADT arrayPool;

/*--------------------------------------
 * Function: GetPointArray()
 * Parameters:
 *
 * Description:
 *   Returnerar en array med pekare till punkter. Antingen skapar funktionen en
 *   ny array, eller så återanvänder den en array från poolen.
 *------------------------------------*/
static arrayADT GetPointArray() {
#   ifdef ArrayPooling
        if (!arrayPool)
            arrayPool = NewQueue(32);

        if (QueueIsEmpty(arrayPool))
            return NewArray(sizeof(pointT *));

        return Dequeue(arrayPool);
#   else
        return NewArray(sizeof(pointT *));
#   endif
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
#   ifdef ArrayPooling
        if (QueueIsFull(arrayPool)) {
            FreeArray(a);
            printf("Warning: Array pool is not big enough.\n");
            return;
        }

        ResetArray(a);
        Enqueue(arrayPool, a);
   #else
        FreeArray(a);
#   endif
}

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

#define _0(x)   (x & 0x7FF)
#define _1(x)   (x >> 11 & 0x7FF)
#define _2(x)   (x >> 22 )
algorithmDataT RaddeSort(pointsetT ps, arrayADT sorted) {
	algorithmDataT algo = { 0 };
	int ops;

	arrayADT tempsort = GetPointArray();
	algo.numAllocs++;

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

		pointT *p = &ps.points[i];
		ArraySet(sorted, ++b0[pos], &p);
		//sorted[++b0[pos]] = &ps.points[i];

		ops++;
	}

	// pass 2, sortera på nästa 11 bits
	for (int i = 0; i < ps.numPoints; i++) {
		pointT *p = *(pointT **)ArrayGet(sorted, i);
		ui si = *(ui*)&p->x;
		ui pos = _1(si);

		ArraySet(tempsort, ++b1[pos], &p);
		//tempsort[++b1[pos]] = sorted[i];
		ops++;
	}

	// pass 3, sortera på sista 10 bits
	for (int i = 0; i < ps.numPoints; i++) {
		pointT *p = *(pointT **)ArrayGet(tempsort, i);
		ui *ai = (ui*)&p->x;
		ui pos = _2(*ai);

		IFloatFlipX(ai);

		ArraySet(sorted, ++b2[pos], &p);
		//sorted[++b2[pos]] = tempsort[i];
		ops++;
	}

	free(b0);

	algo.numBytes += ArrayBytes(tempsort);
	ReleaseArray(tempsort);

	algo.numOps = ops;
	return algo;
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

void InitHullPoints(pointT *p1, pointT *p2) {
	p1->next = p2;
	p1->prev = p2;

	p2->next = p1;
	p2->prev = p1;
}

/*--------------------------------------
 * Function: Incremental()
 * Parameters:
 *   ps    Punktuppsättningen för vilken ett hölje ska genereras.
 *   hull  En pekare till höljet.
 *
 * Description:
 *   Genererar att konvext hölje för punktuppsättningen med hjälp av algoritmen
 *   incremental convex hull. Returnerar data om algoritmens arbete.
 *------------------------------------*/
algorithmDataT Incremental(pointsetT ps, hullT *hull) {
    algorithmDataT algo = { 0 };
	
	arrayADT sorted = GetPointArray();
	algo.numAllocs++;
	
	algorithmDataT algoRadde = RaddeSort(ps, sorted);
	algo.numAllocs += algoRadde.numAllocs;
	algo.numBytes += algoRadde.numBytes;
	
	int ops = algoRadde.numOps;
	
	pointT
		*p1 = *(pointT**)ArrayGet(sorted, 0),
		*p2 = *(pointT**)ArrayGet(sorted, 1),
		*LastHP = p2,
		*P;

	InitHullPoints(p1, p2);
	// ta in nästa punkt, den är med i hull (den är längst till höger), ta reda på vilka tidigare hullpoints som ska bort
	for (int i = 2; i < ps.numPoints; i++) {
		pointT *A = *(pointT**)ArrayGet(sorted, i);

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

	algo.numBytes += ArrayBytes(sorted);
	ReleaseArray(sorted);

	algo.numOps = ops;

    return algo;
}
