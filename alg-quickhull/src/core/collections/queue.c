/*------------------------------------------------------------------------------
 * File: queue.c
 * Created: May 23, 2015
 * Last changed: May 24, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Erbjuder en k�-struktur. Implementationen anv�nder sig av en cirkul�r
 *   vektor f�r att slippa omallokeringar.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "queue.h"

#include "core/common.h"
#include "core/debug.h"

#include <stdlib.h>

/*------------------------------------------------
 * TYPES
 *----------------------------------------------*/

/*--------------------------------------
 * Type: queueCDT
 *
 * Description:
 *   Representerar en k� med objekt i.
 *------------------------------------*/
typedef struct queueCDT {
    void **data;
    int    head, tail;
    int    size;
} queueCDT;

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: NewQueue()
 * Parameters:
 *   size  Det maximala antalet objekt som k�n kan inneh�lla.
 *
 * Description:
 *   Skapar en ny k� av den givna storleken.
 *------------------------------------*/
queueADT NewQueue(int size) {
    Assert(size > 0);

    // Vi �kar storleken med ett eftersom vi beh�ver ett tomt element i vektorn.
    // Anledningen �r att vi har en cirkelbuffert. Utan ett extra, tomt element
    // kan vi inte skilja p� full och tom k�.
    size++;

    queueADT queue = malloc(sizeof(queueCDT));

    queue->data = malloc(sizeof(void *) * size);
    queue->head = 0;
    queue->tail = 0;
    queue->size = size;

    return queue;
}

/*--------------------------------------
 * Function: FreeQueue()
 * Parameters:
 *   queue  K�n som ska frig�ras fr�n minnet.
 *
 * Description:
 *   Frig�r den specificerade k�n fr�n minnet.
 *------------------------------------*/
void FreeQueue(queueADT queue) {
    while (!QueueIsEmpty(queue))
        Dequeue(queue);

    free(queue);
}

/*--------------------------------------
 * Function: Enqueue()
 * Parameters:
 *   queue  K�n till vilket ett v�rde ska l�ggas.
 *   value  V�rdet som ska l�ggas till i k�n.
 *
 * Description:
 *   L�gger till ett v�rde i en k�.
 *------------------------------------*/
void Enqueue(queueADT queue, void *value) {
    Assert(!QueueIsFull(queue));

    queue->data[queue->tail] = value;

    queue->tail++;
    if (queue->tail >= queue->size)
        queue->tail = 0;
}

/*--------------------------------------
 * Function: Dequeue()
 * Parameters:
 *   queue  K�n fr�n vilket ett v�rde ska tas.
 *
 * Description:
 *   Tar ut det f�rsta v�rdet i k�n och returnerar det.
 *------------------------------------*/
void *Dequeue(queueADT queue) {
    Assert(!QueueIsEmpty(queue));

    void *value = queue->data[queue->head];

    queue->head++;
    if (queue->head >= queue->size)
        queue->head = 0;

    return value;
}


/*--------------------------------------
 * Function: QueueIsEmpty()
 * Parameters:
 *   queue  K�n som ska kontrolleras.
 *
 * Description:
 *   Returnerar sant om den specificerade k�n �r tom.
 *------------------------------------*/
bool QueueIsEmpty(queueADT queue) {
    return (queue->head == queue->tail);
}

/*--------------------------------------
 * Function: QueueIsFull()
 * Parameters:
 *   queue  K�n som ska kontrolleras.
 *
 * Description:
 *   Returnerar sant om den specificerade k�n �r full.
 *------------------------------------*/
bool QueueIsFull(queueADT queue) {
    return (queue->head == ((queue->tail+1) % queue->size));
}
