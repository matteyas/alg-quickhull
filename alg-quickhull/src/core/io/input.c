/*------------------------------------------------------------------------------
 * File: input.c
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

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "input.h"

#include "core/common.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
bool GetBoolFromUser(bool defaultVal) {
    string s = GetStringFromUser();
    char   c = s[0];

    free(s);

    if (c=='\0')
        return defaultVal;

    if (c=='Y' || c=='y')
        return TRUE;

    return FALSE;
}

/*--------------------------------------
 * Function: GetIntFromUser()
 * Parameters:
 *
 * Description:
 *   L�ter anv�ndaren skriva in ett heltal.
 *------------------------------------*/
int GetIntFromUser() {
    char buf[16];
    bool validInt = FALSE;

    while (!validInt) {
        validInt = TRUE;

        fgets(buf, sizeof(buf), stdin);

        int len = strlen(buf);
        for (int i = 0; i < len; i++) {
            char c = buf[i];

            if (c=='\r' || c=='\n') {
                buf[i] = '\0';
                break;
            }

            if (!isdigit(c)) {
                validInt = FALSE;
                break;
            }
        }

        if (strlen(buf) == 0)
            validInt = FALSE;

        if (!validInt)
            printf("Invalid integer. Try again: ");
    }

    return atoi(buf);
}

/*--------------------------------------
 * Function: GetStringFromUser()
 * Parameters:
 *
 * Description:
 *   L�ter anv�ndaren skriva in en str�ng. Gl�m inte anropa free() efter�t f�r
 *   att f�rhindra minnesl�ckage.
 *------------------------------------*/
string GetStringFromUser() {
    char buf[1024];

    fgets(buf, sizeof(buf), stdin);

    int len = strlen(buf);
    for (int i = 0; i < len; i++) {
        char c = buf[i];
        if (c == '\r' || c == '\n') {
            buf[i] = '\0';
            break;
        }
    }

    return strdup(buf);
}
