#ifndef IO_H_
#define IO_H_

//------------------------------------------------
// INCLUDES
//------------------------------------------------

#include "common.h"

//------------------------------------------------
// FUNCTIONS
//------------------------------------------------

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

#endif // IO_H_
