/*------------------------------------------------------------------------------
 * File: quickhull.c
 * Created: May 28, 2015
 * Last changed: May 28, 2015
 *
 * Author(s): "Philip Arvidsson" (philip@philiparvidsson.com)
 *
 * Description:
 *   Inneh�ller algoritmen quickhull f�r att l�sa det konvexa h�ljet.
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
 *   Kommentera ut raden nedan f�r att inte anv�nda array-poolen.
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
 *   En k� d�r vi poolar arrayer f�r att slippa omallokeringar. Vi �teranv�nder
 *   helt enkelt minne ist�llet.
 *------------------------------------*/
static queueADT arrayPool;

/*--------------------------------------
 * Function: GetPointArray()
 * Parameters:
 *
 * Description:
 *   Returnerar en array med pekare till punkter. Antingen skapar funktionen en
 *   ny array, eller s� �teranv�nder den en array fr�n poolen.
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
 *   a  Arrayen som ska sl�ppas tillbaka till array-poolen.
 *
 * Description:
 *   Sl�pper tillbaka en array till poolen.
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

	// konstruera histogram (det �r typ counters p� antal tal i olika storleksspann)

	ui *b0 = calloc(3 * kHist, 4);
	ui *b1 = b0 + kHist;
	ui *b2 = b1 + kHist;

	algo.numAllocs++;
	algo.numBytes = bytes_kHist;

	// pass 0, ber�kna histograms
	for (ops = 0; ops < ps.numPoints; ops++) {
		ui fi = FloatFlip(*(ui*)&ps.points[ops].x);

		b0[_0(fi)] ++;
		b1[_1(fi)] ++;
		b2[_2(fi)] ++;
	}

	// g�r n�n supersummering, shrug... offsets typ
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

	// pass 1, sortera p� f�rsta 11 bits
	for (int i = 0; i < ps.numPoints; i++) {
		ui *fi = (ui*)&ps.points[i].x;
		FloatFlipX(fi);
		ui pos = _0(*fi);

		int j = ++b0[pos];
		int len = ArrayLength(sorted);
		printf("%d\n", len);
		ArraySet(sorted, j, &ps.points[i]);
		//sorted[++b0[pos]] = &ps.points[i];
		
		len = ArrayLength(sorted);
		printf("%d\n", len);
		printf("idx: %d\n", j);
		
		pointT *p = *(pointT **)ArrayGet(sorted, j);
		printf("%f\n", p->x);
		ops++;
	}

	// pass 2, sortera p� n�sta 11 bits
	for (int i = 0; i < ps.numPoints; i++) {
		pointT *p = *(pointT **)ArrayGet(sorted, i);
		printf("%f\n", p->x);
		ui si = *(ui*)&p->x;
		ui pos = _1(si);

		ArraySet(tempsort, ++b1[pos], &p);
		//tempsort[++b1[pos]] = sorted[i];
		ops++;
	}

	// pass 3, sortera p� sista 10 bits
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

// kollar om linjen AP ligger till h�ger om AB (perspektiv fr�n A till B)
__inline char Right(pointT *A, pointT *B, pointT *P) {
	float d = (B->x - A->x) * (P->y - A->y)
		- (B->y - A->y) * (P->x - A->x);

	if (d < 0.0f)
		return 1;
	else
		return 0;
}

// kollar om linjen AP ligger till v�nster om AB (perspektiv fr�n A till B)
__inline char Left(pointT *A, pointT *B, pointT *P) {
	float d = (B->x - A->x) * (P->y - A->y)
		- (B->y - A->y) * (P->x - A->x);

	if (d > 0.0f)
		return 1;
	else
		return 0;
}

/*--------------------------------------
 * Function: Incremental()
 * Parameters:
 *   ps    Punktupps�ttningen f�r vilken ett h�lje ska genereras.
 *   hull  En pekare till h�ljet.
 *
 * Description:
 *   Genererar att konvext h�lje f�r punktupps�ttningen med hj�lp av algoritmen
 *   incremental convex hull. Returnerar data om algoritmens arbete.
 *------------------------------------*/
algorithmDataT Incremental(pointsetT ps, hullT *hull) {
    algorithmDataT algo = { 0 };
	arrayADT sorted = GetPointArray();
	RaddeSort(ps, sorted);
	
	for (int i = 0; i < ps.numPoints; i += ps.numPoints / 10) {
		pointT *p = ArrayGet(sorted, i);
		printf("%f\n", p->x);
	}

    return algo;
}
