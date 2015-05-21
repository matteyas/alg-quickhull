/*------------------------------------------------------------------------------
 * File: graphics.h
 * Created: May 15, 2015
 * Last changed: May 15, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Enkelt grafikbibliotek f�r uppritning av simpel geometri. H�r anv�nder vi
 *   OpenGL.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef _graphics_h_
#define _graphics_h_

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "common.h"
#include "math.h"

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Grundl�ggande funktioner f�r att initiera grafikf�nstret, st�lla in        */
/* egenskaper, med mera.                                                      */
/*----------------------------------------------------------------------------*/

/*--------------------------------------
 * Function: InitGraphics()
 * Parameters:
 *
 * Description:
 *   Initierar grafikbiblioteket.
 *------------------------------------*/
void InitGraphics();

/*--------------------------------------
 * Function: SetFrameRate()
 * Parameters:
 *   fps  Antalet bildrutor per sekund.
 *
 * Description:
 *   St�ller in hur m�nga bildrutor per sekund som ska ritas.
 *------------------------------------*/
void SetFrameRate(int fps);

/*--------------------------------------
 * Function: SetColor()
 * Parameters:
 *   red    R�d f�rgkomponent.
 *   green  Gr�n f�rgkomponent.
 *   blue   Bl� f�rgkomponent.
 *   alpha  Alphav�rde.
 *
 * Description:
 *   St�ller in f�rgen f�r n�sta ritoperation.
 *------------------------------------*/
void SetColor(float red, float green, float blue, float alpha);

/*--------------------------------------
 * Function: SetLineWidth()
 * Parameters:
 *   width  Bredden p� linjer, i pixlar.
 *
 * Description:
 *   St�ller in bredden p� linjer.
 *------------------------------------*/
void SetLineWidth(float width);

/*--------------------------------------
 * Function: WindowIsOpen()
 * Parameters:
 *
 * Description:
 *   Returnerar sant om grafikf�nstret �r �ppet.
 *------------------------------------*/
bool WindowIsOpen();

/*----------------------------------------------------------------------------*/
/* Ritfunktioner f�r att rensa ritytan, rita punkter, linjer, trianglar m.m.  */
/*----------------------------------------------------------------------------*/

/*--------------------------------------
 * Function: ClearCanvas()
 * Parameters:
 *   red    R�d f�rgkomponent.
 *   green  Gr�n f�rgkomponent.
 *   blue   Bl� f�rgkomponent.
 *
 * Description:
 *   Rensar ritytan till den specificerade f�rgen.
 *------------------------------------*/
void ClearCanvas(float red, float green, float blue);

/*--------------------------------------
 * Function: DrawHull()
 * Parameters:
 *   hull  Det h�lje som ska ritas upp.
 *
 * Description:
 *   Ritar upp ett h�lje.
 *------------------------------------*/
void DrawHull(hullT hull);

/*--------------------------------------
 * Function: DrawLine()
 * Parameters:
 *   line  Den linje som ska ritas upp.
 *
 * Description:
 *   Ritar upp en linje.
 *------------------------------------*/
void DrawLine(lineT line);

/*--------------------------------------
 * Function: DrawPoint()
 * Parameters:
 *   p  Punkten som ska ritas.
 *
 * Description:
 *   Ritar en punkt p� ritytan.
 *------------------------------------*/
void DrawPoint(pointT p);

/*--------------------------------------
 * Function: DrawPoints()
 * Parameters:
 *   ps  Den upps�ttning punkter som ska ritas upp.
 *
 * Description:
 *   Ritar en upps�ttning punkter.
 *------------------------------------*/
void DrawPoints(pointsetT ps);

/*--------------------------------------
 * Function: UpdateDisplay()
 * Parameters:
 *
 * Description:
 *   Uppdaterar ritytan.
 *------------------------------------*/
void UpdateDisplay();

/*----------------------------------------------------------------------------*/
/* Funktioner f�r input, t ex associera knappar med callback-funktioner m.m.  */
/*----------------------------------------------------------------------------*/

/*--------------------------------------
 * Function: OnKeyPress()
 * Parameters:
 *   c    Tecknet som ska associeras till en funktion.
 *   cb   Callback-funktionen som ska anropas n�r knappen med det specificerade
 *        tecknet trycks ned.
 *   arg  Det argument som ska skickas till callback-funktionen.
 *
 * Description:
 *   Registerar en callback-funktion f�r en viss knapp.
 *------------------------------------*/
void OnKeyPress(char c, actionT cb, void *arg);

#endif // _graphics_h_
