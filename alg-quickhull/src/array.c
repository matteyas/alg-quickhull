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
 * Constant: InitialMaxElements
 *
 * Description:
 *   Det initiala antalet element som minne ska allokeras f�r i Array_Init().
 *------------------------------------*/
#define InitialMaxElements 8 // 8 �r nog lagom.

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

static void DoubleArrayCapacity(arrayT *array) {
    // Arrayen �r full, s� vi dubblar kapaciteten och kopierar �ver de gamla
    // elementen till den nya minnesplatsen, sen sl�pper vi den gamla
    // arrayen ur minnet.

    void *oldElements    = array->elements;
    int   oldMaxElements = array->maxElements;

    array->maxElements *= 2;
    array->elements = malloc(array->maxElements * array->elementSize);
    memcpy(array->elements, oldElements, oldMaxElements * array->elementSize);
    free(oldElements);
}

/*--------------------------------------
 * Function: ArrayAdd()
 * Parameters:
 *   array  Den array till vilken vi ska l�gga ett element.
 *   value  Elementet som ska l�ggas till i arrayen.
 *
 * Description:
 *   L�gger till ett element i en array. Returnerar minnesadressen d�r noden
 *   lades in.
 *------------------------------------*/
void *ArrayAdd(arrayT *array, const void *value) {
    if (array->numElements >= array->maxElements)
        DoubleArrayCapacity(array);

    void *dest = (char *)array->elements + (array->numElements * array->elementSize);
    memcpy(dest, value, array->elementSize);

    array->numElements++;
    return dest;
}

/*--------------------------------------
 * Function: ArrayGet()
 * Parameters:
 *   array  Den array fr�n vilken vi ska l�sa ett element.
 *   i      Det index i arrayen fr�n vilket vi ska l�sa elementet.
 *
 * Description:
 *   L�ser ut och returnerar en pekare till det specificerade elementet i
 *   arrayen.
 *------------------------------------*/
void *ArrayGet(const arrayT *array, int index) {
    Assert(0 <= index && index < array->numElements);
    return (char *)array->elements + (index * array->elementSize);
}

/*--------------------------------------
 * Function: ArrayInsert()
 * Parameters:
 *   array  Den array till vilken vi ska l�gga ett element.
 *   index  Det index i arrayen d�r elementet ska s�ttas in.
 *   value  Elementet som ska l�ggas till i arrayen.
 *
 * Description:
 *   L�gger in ett element i en array vid det specificerade indexet. Returnerar
 *   minnesadressen d�r noden lades in.
 *------------------------------------*/
void *ArrayInsert(arrayT *array, int index, const void *value) {
    Assert(0 <= index && index <= array->numElements);

    if (index == array->numElements)
        return ArrayAdd(array, value);

    if (array->numElements >= array->maxElements)
        DoubleArrayCapacity(array);

    int elementSize = array->elementSize;

    char *base = (char *)array->elements;
    for (int i = array->numElements-1; i > index; i--) {
        void *src  = base + ((i-1) * elementSize);
        void *dest = base + ( i    * elementSize);
        memcpy(dest, src, elementSize);
    }

    void *dest = base + (index * elementSize);
    memcpy(dest, value, elementSize);

    return dest;
}

/*--------------------------------------
 * Function: ArrayLength()
 * Parameters:
 *   array  Arrayen vars l�ngd ska l�sas ut.
 *
 * Description:
 *   Returnerar den specificerade arrayens l�ngd.
 *------------------------------------*/
int ArrayLength(const arrayT *array) {
    return array->numElements;
}

/*--------------------------------------
 * Function: FreeArray()
 * Parameters:
 *   array  Den array som ska avallokeras.
 *
 * Description:
 *   Sl�pper en array ur minnet.
 *------------------------------------*/
void FreeArray(arrayT *array) {
    Assert(array->elements != NULL);

    free(array->elements);

    array->elements    = NULL;
    array->numElements = 0;
    array->maxElements = 0;
    array->elementSize = 0;
}

/*--------------------------------------
 * Function: InitArray()
 * Parameters:
 *   array  Den array som ska initieras.
 *
 * Description:
 *   Initialiserar och allokerar en array.
 *------------------------------------*/
void InitArray(arrayT *array, size_t elementSize) {
    array->numElements = 0;
    array->maxElements = InitialMaxElements;
    array->elementSize = elementSize;

    array->elements = malloc(array->maxElements * array->elementSize);
}
