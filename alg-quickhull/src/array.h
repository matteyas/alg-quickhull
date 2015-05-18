/*------------------------------------------------------------------------------
 * File: array.h
 * Created: January 2, 2015
 * Last changed: February 19, 2015
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

#ifndef Array_h
#define Array_h

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "common.h"
#include "debug.h"

#include <stdlib.h>

/*------------------------------------------------
 * TYPES
 *----------------------------------------------*/

/*--------------------------------------
 * Type: Array
 *
 * Description:
 *   Representerar en dynamisk array med objekt i.
 *------------------------------------*/
typedef struct {
    void  *elements;
    int    numElements;
    int    maxElements;
    size_t elementSize;
} arrayT;

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

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
void* ArrayAdd(arrayT *array, const void *value);

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
void *ArrayGet(const arrayT *array, int i);

/*--------------------------------------
 * Function: ArrayLength()
 * Parameters:
 *   array  Arrayen vars l�ngd ska l�sas ut.
 *
 * Description:
 *   Returnerar den specificerade arrayens l�ngd.
 *------------------------------------*/
int ArrayLength(const arrayT *array);

/*--------------------------------------
 * Function: FreeArray()
 * Parameters:
 *   array  Den array som ska avallokeras.
 *
 * Description:
 *   Sl�pper en array ur minnet.
 *------------------------------------*/
void FreeArray(arrayT *array);

/*--------------------------------------
 * Function: InitArray()
 * Parameters:
 *   array  Den array som ska initieras.
 *
 * Description:
 *   Initialiserar och allokerar en array.
 *------------------------------------*/
void InitArray(arrayT *array, size_t elementSize);

#endif // Array_h
