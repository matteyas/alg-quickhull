/*------------------------------------------------------------------------------
 * File: array.c
 * Created: January 2, 2015
 * Last changed: January 9, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Erbjuder en dynamisk array som v�xer automatiskt med antalet objekt som
 *   l�ggs in.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "array.h"
#include "debug.h"

#include <stdlib.h>
#include <string.h> // memcpy()

/*------------------------------------------------
 * CONSTANTS
 *----------------------------------------------*/

/*--------------------------------------
 * Type: arrayCDT
 *
 * Description:
 *   Representerar en dynamisk array med objekt i.
 *------------------------------------*/
typedef struct arrayCDT {
    void  *elements;
    int    numElements;
    int    maxElements;
    size_t elementSize;
} arrayCDT;

/*------------------------------------------------
 * CONSTANTS
 *----------------------------------------------*/

/*--------------------------------------
 * Constant: InitialMaxElements
 *
 * Description:
 *   Det initiala antalet element som minne ska allokeras f�r i Array_Init().
 *------------------------------------*/
#define InitialMaxElements 8 // 8 �r nog lagom.

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

static void DoubleArrayCapacity(arrayADT a) {
    // Arrayen �r full, s� vi dubblar kapaciteten och kopierar �ver de gamla
    // elementen till den nya minnesplatsen, sen sl�pper vi den gamla
    // arrayen ur minnet.

    void *oldElements    = a->elements;
    int   oldMaxElements = a->maxElements;

    a->maxElements *= 2;
    a->elements = malloc(a->maxElements * a->elementSize);
    memcpy(a->elements, oldElements, oldMaxElements * a->elementSize);
    free(oldElements);
}

/*--------------------------------------
 * Function: ArrayAdd()
 * Parameters:
 *   a      Den array till vilken vi ska l�gga ett element.
 *   value  Elementet som ska l�ggas till i arrayen.
 *
 * Description:
 *   L�gger till ett element i en array. Returnerar minnesadressen d�r noden
 *   lades in.
 *------------------------------------*/
void *ArrayAdd(arrayADT a, const void *value) {
    if (a->numElements >= a->maxElements)
        DoubleArrayCapacity(a);

    void *dest = (char *)a->elements + (a->numElements * a->elementSize);
    memcpy(dest, value, a->elementSize);

    a->numElements++;
    return dest;
}

/*--------------------------------------
 * Function: ArrayBytes()
 * Parameters:
 *   a  Arrayen vars minnesanv�ndning ska r�knas ut.
 *
 * Description:
 *   Returnerar den specificerade arrayens minnesanv�ndning, i antal bytes.
 *------------------------------------*/
int ArrayBytes(arrayADT a) {
    return sizeof(arrayCDT) + a->maxElements*a->elementSize;
}

/*--------------------------------------
 * Function: ArrayGet()
 * Parameters:
 *   a      Den array fr�n vilken vi ska l�sa ett element.
 *   i      Det index i arrayen fr�n vilket vi ska l�sa elementet.
 *
 * Description:
 *   L�ser ut och returnerar en pekare till det specificerade elementet i
 *   arrayen.
 *------------------------------------*/
void *ArrayGet(arrayADT a, int index) {
    Assert(0 <= index && index < a->numElements);
    return (char *)a->elements + (index * a->elementSize);
}

/*--------------------------------------
 * Function: ArrayInsert()
 * Parameters:
 *   a      Den array till vilken vi ska l�gga ett element.
 *   index  Det index i arrayen d�r elementet ska s�ttas in.
 *   value  Elementet som ska l�ggas till i arrayen.
 *
 * Description:
 *   L�gger in ett element i en array vid det specificerade indexet. Returnerar
 *   minnesadressen d�r noden lades in.
 *------------------------------------*/
void *ArrayInsert(arrayADT a, int index, const void *value) {
    Assert(0 <= index && index <= a->numElements);

    if (index == a->numElements)
        return ArrayAdd(a, value);

    if (a->numElements >= a->maxElements)
        DoubleArrayCapacity(a);

    int elementSize = a->elementSize;

    char *base = (char *)a->elements;
    for (int i = a->numElements; i > index; i--) {
        void *src  = base + ((i-1) * elementSize);
        void *dest = base + ( i    * elementSize);
        memcpy(dest, src, elementSize);
    }

    void *dest = base + (index*elementSize);
    memcpy(dest, value, elementSize);

    a->numElements++;
    return dest;
}

/*--------------------------------------
 * Function: ArrayLength()
 * Parameters:
 *   a  Arrayen vars l�ngd ska l�sas ut.
 *
 * Description:
 *   Returnerar den specificerade arrayens l�ngd.
 *------------------------------------*/
int ArrayLength(arrayADT a) {
    return a->numElements;
}

/*--------------------------------------
 * Function: FreeArray()
 * Parameters:
 *   a  Den array som ska avallokeras.
 *
 * Description:
 *   Sl�pper en array ur minnet.
 *------------------------------------*/
void FreeArray(arrayADT a) {
    Assert(a->elements != NULL);

    free(a->elements);
    free(a);
}

/*--------------------------------------
 * Function: NewArray()
 * Parameters:
 *   elementSize  Storleken, i bytes, p� arrayens element.
 *
 * Description:
 *   Initierar och allokerar en array.
 *------------------------------------*/
arrayADT NewArray(size_t elementSize) {
    arrayADT array = malloc(sizeof(arrayCDT));

    array->numElements = 0;
    array->maxElements = InitialMaxElements;
    array->elementSize = elementSize;

    array->elements = malloc(array->maxElements * array->elementSize);

    return array;
}
