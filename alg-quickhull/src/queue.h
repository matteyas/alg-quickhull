/*------------------------------------------------------------------------------
 * File: queue.h
 * Created: May 23, 2015
 * Last changed: May 24, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Erbjuder en k�-struktur. Implementationen anv�nder sig av en cirkul�r
 *   buffer f�r att slippa omallokeringar.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef _queue_h_
#define _queue_h_

/*------------------------------------------------
 * TYPES
 *----------------------------------------------*/

/*--------------------------------------
 * Type: queueADT
 *
 * Description:
 *   Representerar en k� med objekt i.
 *------------------------------------*/
typedef struct queueCDT *queueADT;

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: NewQueue()
 * Parameters:
 *   size  Det maximala antalet objekt som k�n kan inneh�lla.
 *
 * Description:
 *   Sakapar en ny k� av den givna storleken och returnerar den.
 *------------------------------------*/
queueADT NewQueue(int size);

/*--------------------------------------
 * Function: FreeQueue()
 * Parameters:
 *   queue  K�n som ska frig�ras fr�n minnet.
 *
 * Description:
 *   Frig�r den specificerade k�n fr�n minnet.
 *------------------------------------*/
void FreeQueue(queueADT queue);

void Enqueue(queueADT queue, void *value);
void *Dequeue(queueADT queue);

bool QueueIsEmpty(queueADT queue);
bool QueueIsFull(queueADT queue);

#endif // _queue_h_
