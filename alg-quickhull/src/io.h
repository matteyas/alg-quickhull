/*------------------------------------------------------------------------------
 * File: io.h
 * Created: May 15, 2015
 * Last changed: May 15, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Input-funktioner f�r att l�sa in data fr�n anv�ndaren.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef IO_h
#define IO_h

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "common.h"

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: GetBoolFromUser()
 * Parameters:
 *   defaultVal  V�rdet som ska returneras om anv�ndaren inte svarar.
 *
 * Description:
 *   L�ter anv�ndaren skriva in ja eller nej.
 *------------------------------------*/
bool GetBoolFromUser(bool defaultVal);

/*--------------------------------------
 * Function: GetIntFromUser()
 * Parameters:
 *
 * Description:
 *   L�ter anv�ndaren skriva in ett heltal.
 *------------------------------------*/
int GetIntFromUser();

/*--------------------------------------
 * Function: GetStringFromUser()
 * Parameters:
 *
 * Description:
 *   L�ter anv�ndaren skriva in en str�ng. Gl�m inte anropa free() efter�t f�r
 *   att f�rhindra minnesl�ckage.
 *------------------------------------*/
string GetStringFromUser();

#endif // IO_h
