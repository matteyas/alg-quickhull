/*------------------------------------------------------------------------------
 * File: graphics.c
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

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "graphics.h"

#include "core/common.h"
#include "core/debug.h"
#include "core/stopwatch.h"

#include "extra/sandbox.h"

#ifdef _WIN32
#include <windows.h>
#endif
#include <GL/gl.h>

#pragma comment(lib, "opengl32.lib")

/*------------------------------------------------
 * EXTERNALS
 *----------------------------------------------*/

extern void  InitGraphicsImpl();
extern void    OnKeyPressImpl(char c, actionT cb, void *arg);
extern void UpdateDisplayImpl();
extern bool  WindowIsOpenImpl();

/*------------------------------------------------
 * CONSTANTS
 *----------------------------------------------*/

/*--------------------------------------
 * Constant: FrameRateStopwatchID
 *
 * Description:
 *   ID f�r tidtagaruret som anv�nds f�r att synkronisera antalet bildrutor som
 *   visas varje sekund.
 *------------------------------------*/
#define FrameRateStopwatchID 93

/*------------------------------------------------
 * GLOBALS
 *----------------------------------------------*/

/*--------------------------------------
 * Variable: initialized
 *
 * Description:
 *   Indikerar om grafiksystemet �r initierat.
 *------------------------------------*/
static bool initialized;

/*--------------------------------------
 * Variable: microSecsPerFrame
 *
 * Description:
 *   Antal mikrosekunder som varje bildruta ska visas.
 *------------------------------------*/
static int microSecsPerFrame;

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*--------------------------------------
 * Function: CheckInitGraphics()
 * Parameters:
 *
 * Description:
 *   Anropar Error()-funktionen om grafibiblioteket inte �r initierat.
 *------------------------------------*/
static void CheckInitGraphics() {
    if (!initialized)
        Error("InitGraphics() must be called first");
}

/*--------------------------------------
 * Function: InitOpenGL()
 * Parameters:
 *
 * Description:
 *   Initierar OpenGL.
 *------------------------------------*/
static void InitOpenGL() {
    // Utan GL_BLEND fungerar inte kantutj�mningen f�r linjer.
    glEnable   (GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Runda, fina prickar! :-)
    glEnable   (GL_POINT_SMOOTH);
    glHint     (GL_POINT_SMOOTH_HINT, GL_NICEST);
    glPointSize(8.0f);

    // Mjuka, fina linjer!
    glEnable   (GL_LINE_SMOOTH);
    glHint     (GL_LINE_SMOOTH_HINT, GL_NICEST);
    glLineWidth(2.0f);
}

/*----------------------------------------------------------------------------*/
/* Grundl�ggande funktioner f�r att initiera grafikf�nstret, st�lla in        */
/* egenskaper, med mera.                                                      */
/*----------------------------------------------------------------------------*

/*--------------------------------------
 * Function: InitGraphics()
 * Parameters:
 *
 * Description:
 *   Initierar grafikbiblioteket.
 *------------------------------------*/
void InitGraphics() {
    if (initialized)
        Error("Graphics already initialized");

    InitGraphicsImpl();
    InitOpenGL();

    initialized = TRUE;

    SetFrameRate(30);
}

/*--------------------------------------
 * Function: SetFrameRate()
 * Parameters:
 *   fps  Antalet bildrutor per sekund.
 *
 * Description:
 *   St�ller in hur m�nga bildrutor per sekund som ska ritas.
 *------------------------------------*/
void SetFrameRate(int fps) {
    CheckInitGraphics();

    microSecsPerFrame = 1000000 / fps;
    ResetStopwatch(FrameRateStopwatchID);
}

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
void SetColor(float red, float green, float blue, float alpha) {
    CheckInitGraphics();

    glColor4f(red, green, blue, alpha);
}

/*--------------------------------------
 * Function: SetLineWidth()
 * Parameters:
 *   width  Bredden p� linjer, i pixlar.
 *
 * Description:
 *   St�ller in bredden p� linjer.
 *------------------------------------*/
void SetLineWidth(float width) {
    CheckInitGraphics();

         if (width <  0.0f) width =  0.0f;
    else if (width > 10.0f) width = 10.0f;

    glLineWidth(width);
}

/*--------------------------------------
 * Function: WindowIsOpen()
 * Parameters:
 *
 * Description:
 *   Returnerar sant om grafikf�nstret �r �ppet.
 *------------------------------------*/
bool WindowIsOpen() {
    CheckInitGraphics();

    return WindowIsOpenImpl();
}

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
void ClearCanvas(float red, float green, float blue) {
    CheckInitGraphics();

    glClearColor(red, green, blue, 1.0f);
    glClear     (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

/*--------------------------------------
 * Function: DrawHull()
 * Parameters:
 *   hull  Det h�lje som ska ritas upp.
 *
 * Description:
 *   Ritar upp ett h�lje.
 *------------------------------------*/
void DrawHull(hullT hull) {
    CheckInitGraphics();

    for (int i = 0; i < hull.numLines; i++)
        DrawLine(hull.lines[i]);
}

/*--------------------------------------
 * Function: DrawLine()
 * Parameters:
 *   line  Den linje som ska ritas upp.
 *
 * Description:
 *   Ritar upp en linje.
 *------------------------------------*/
void DrawLine(lineT line) {
    CheckInitGraphics();

    glBegin(GL_LINES);
    glVertex2f(line.a->x, line.a->y);
    glVertex2f(line.b->x, line.b->y);
    glEnd     ();
}

/*--------------------------------------
 * Function: DrawPoint()
 * Parameters:
 *   p  Punkten som ska ritas.
 *
 * Description:
 *   Ritar en punkt p� ritytan.
 *------------------------------------*/
void DrawPoint(pointT p) {
    CheckInitGraphics();

    glBegin   (GL_POINTS);
    glVertex2f(p.x, p.y);
    glEnd     ();
}

/*--------------------------------------
 * Function: DrawPoints()
 * Parameters:
 *   ps  Den upps�ttning punkter som ska ritas upp.
 *
 * Description:
 *   Ritar en upps�ttning punkter.
 *------------------------------------*/
void DrawPoints(pointsetT ps) {
    CheckInitGraphics();

    for (int i = 0; i < ps.numPoints; i++) {
        pointT *p = &ps.points[i];

        glBegin   (GL_POINTS);
        glVertex2f(p->x, p->y);
        glEnd     ();
    }
}

/*--------------------------------------
 * Function: UpdateDisplay()
 * Parameters:
 *
 * Description:
 *   Uppdaterar ritytan.
 *------------------------------------*/
void UpdateDisplay() {
    CheckInitGraphics();

    UpdateDisplayImpl();

    // H�r pausar vi tillr�ckligt l�nge f�r att synkronisera hastigheten som
    // bildrutorna visas i.
    while (StopwatchElapsed(FrameRateStopwatchID) < microSecsPerFrame) {
        // Wiiie, h�r snurrar det p� rej�lt!
    }

    ResetStopwatch(FrameRateStopwatchID);
}

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
void OnKeyPress(char c, actionT cb, void *arg) {
    CheckInitGraphics();

    OnKeyPressImpl(c, cb, arg);
}
