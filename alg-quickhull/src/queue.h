/*------------------------------------------------------------------------------
 * File: queue.h
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
 *   Skapar en ny k� av den givna storleken och returnerar den.
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

/*--------------------------------------
 * Function: Enqueue()
 * Parameters:
 *   queue  K�n till vilket ett v�rde ska l�ggas.
 *   value  V�rdet som ska l�ggas till i k�n.
 *
 * Description:
 *   L�gger till ett v�rde i en k�.
 *------------------------------------*/
void Enqueue(queueADT queue, void *value);

/*--------------------------------------
 * Function: Dequeue()
 * Parameters:
 *   queue  K�n fr�n vilket ett v�rde ska tas.
 *
 * Description:
 *   Tar ut det f�rsta v�rdet i k�n och returnerar det.
 *------------------------------------*/
void *Dequeue(queueADT queue);

/*--------------------------------------
 * Function: QueueIsEmpty()
 * Parameters:
 *   queue  K�n som ska kontrolleras.
 *
 * Description:
 *   Returnerar sant om den specificerade k�n �r tom.
 *------------------------------------*/
bool QueueIsEmpty(queueADT queue);

/*--------------------------------------
 * Function: QueueIsFull()
 * Parameters:
 *   queue  K�n som ska kontrolleras.
 *
 * Description:
 *   Returnerar sant om den specificerade k�n �r full.
 *------------------------------------*/
bool QueueIsFull(queueADT queue);

#endif // _queue_h_
