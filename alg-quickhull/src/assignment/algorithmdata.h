/*------------------------------------------------------------------------------
 * File: algorithmdata.h
 * Created: May 21, 2015
 * Last changed: May 21, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Datatyp f�r att utv�rdera algoritmer.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef _algorithmdata_h_
#define _algorithmdata_h_

/*------------------------------------------------
 * TYPES
 *----------------------------------------------*/

/*
 * Type: algorithmDataT
 *
 * Description:
 *   Inneh�ller data om en algoritms arbete.
 */
typedef struct {
    // Antal "kritska operationer" som algoritmen utf�rt. Detta �r menat f�r att
    // f�rst� en algoritms komplexitet, inte f�r att j�mf�ra tv� olika
    // algoritmer.
    int numOps;

    // Antalet minnesallokeringar som algoritmen gjort.
    int numAllocs;

    // Den totala m�ngd minne som algoritmen anv�nt, i bytes.
    int numBytes;
} algorithmDataT;

#endif // _algorithmdata_h_
