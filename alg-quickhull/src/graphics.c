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

#include "common.h"
#include "debug.h"
#include "graphics.h"
#include "sandbox.h"

#include <tchar.h>
#include <Windows.h>

#include <gl/GL.h>

#pragma comment(lib, "opengl32.lib")

/*------------------------------------------------
 * CONSTANTS
 *----------------------------------------------*/

/*--------------------------------------
 * Constant: ClassName
 *
 * Description:
 *   Namnet p� klassen som anv�nds f�r att skapa grafikf�nstret.
 *------------------------------------*/
#define ClassName _T("alg-quickhull")

/*--------------------------------------
 * Constant: WindowTitle
 *
 * Description:
 *   Grafikf�nstrets titel.
 *------------------------------------*/
#define WindowTitle _T("Quickhull Sandbox")

/*------------------------------------------------
 * GLOBALS
 *----------------------------------------------*/

/*--------------------------------------
 * Variable: frameInterval
 *
 * Description:
 *   Anv�nds f�r att lagra intervallet i vilket bildrutor ska visas.
 *------------------------------------*/
static int frameInterval;

/*--------------------------------------
 * Variable: hdc
 *
 * Description:
 *   Grafikf�nstrets DC (device context).
 *------------------------------------*/
static HDC hdc;

/*--------------------------------------
 * Variable: hwnd
 *
 * Description:
 *   Grafikf�nstrets f�nster-handtag (kommentarer p� svenska �ger).
 *------------------------------------*/
static HWND hwnd;

/*--------------------------------------
 * Variable: initialized
 *
 * Description:
 *   Indikerar om grafiksystemet �r initierat.
 *------------------------------------*/
static bool initialized;

/*--------------------------------------
 * Variable: keyPressCB
 *
 * Description:
 *   Callback-funktioner f�r olika knapptryck.
 *------------------------------------*/
static actionT keyPressCB[UCHAR_MAX];

/*--------------------------------------
 * Variable: keyPressArg
 *
 * Description:
 *   Argumentet till de olika callback-funktionerna.
 *------------------------------------*/
static void *keyPressArg[UCHAR_MAX];

/*--------------------------------------
 * Variable: lastUpdate
 *
 * Description:
 *   Tidsst�mpel d� senaste uppdateringen av ritytan gjordes.
 *------------------------------------*/
static LARGE_INTEGER lastUpdate;

/*--------------------------------------
 * Variable: windowOpen
 *
 * Description:
 *   Indikerar om grafikf�nstret �r �ppet.
 *------------------------------------*/
static bool windowOpen;

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* Privata funktioner som bara anv�nds internet av philgraph-biblioteket.     */
/*----------------------------------------------------------------------------*/

/*--------------------------------------
 * Function: WindowProc()
 * Parameters:
 *
 * Description:
 *   http://en.wikipedia.org/wiki/WindowProc
 *------------------------------------*/
LRESULT CALLBACK WindowProc(_In_ HWND   hwnd,
                            _In_ UINT   uMsg,
                            _In_ WPARAM wParam,
                            _In_ LPARAM lParam)
{
    if (uMsg == WM_CLOSE)
        windowOpen = FALSE;

    if (uMsg == WM_CHAR) {
        bool repeated = lParam & 0x40000000;
        if (!repeated) {
            char    c  = tolower(wParam);
            actionT cb = keyPressCB[c];
            if (cb)
                cb(keyPressArg[c]);
        }
    }

    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

static void CheckInitGraphics() {
#ifdef _DEBUG
    if (!initialized)
        Error("InitGraphics() must be called first");
#endif
}

static void InitOpenGL() {
    HGLRC hglrc = wglCreateContext(hdc);

    if (!hglrc)
        Error("wglCreateContext() failed");

    if (!wglMakeCurrent(hdc, hglrc))
        Error("wglMakeCurrent() failed");

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

static void InitPixelFormat() {
    PIXELFORMATDESCRIPTOR pfd;

    pfd.nSize      = sizeof(PIXELFORMATDESCRIPTOR);
    pfd.nVersion   = 1;
    pfd.dwFlags    = PFD_DRAW_TO_WINDOW | PFD_DOUBLEBUFFER | PFD_SUPPORT_OPENGL;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 24;
    pfd.cAlphaBits = 8;
    pfd.cDepthBits = 32;
    pfd.iLayerType = PFD_MAIN_PLANE;

    int pixelFormat = ChoosePixelFormat(hdc, &pfd);
    if (pixelFormat == 0)
        Error("ChoosePixelFormat() failed");

    if (!SetPixelFormat(hdc, pixelFormat, &pfd))
        Error("SetPixelFormat() failed");
}

static void InitWindow() {
    WNDCLASSEX wcx;

    wcx.cbSize        = sizeof(WNDCLASSEX);
    wcx.style         = CS_HREDRAW | CS_VREDRAW;
    wcx.lpfnWndProc   = WindowProc;
    wcx.cbClsExtra    = 0;
    wcx.cbWndExtra    = 0;
    wcx.hInstance     = GetModuleHandle(NULL);
    wcx.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    wcx.hCursor       = LoadCursor(NULL, IDC_ARROW);
    wcx.hbrBackground = (HBRUSH)(COLOR_WINDOW+1);
    wcx.lpszMenuName  = NULL;
    wcx.lpszClassName = ClassName;
    wcx.hIconSm       = NULL;


    if (!RegisterClassEx(&wcx))
        Error("RegisterClassEx() failed");


    DWORD dwStyle = WS_OVERLAPPEDWINDOW & ~WS_THICKFRAME & ~WS_MAXIMIZEBOX;

    RECT rect = { 0, 0, 600, 600 };

    if (!AdjustWindowRectEx(&rect, dwStyle, FALSE, WS_EX_LEFT))
        Error("AdjustWindowRectEx() failed");

    hwnd = CreateWindowEx(WS_EX_LEFT,
                          ClassName,
                          WindowTitle,
                          dwStyle,
                          CW_USEDEFAULT,
                          CW_USEDEFAULT,
                          rect.right - rect.left,
                          rect.bottom - rect.top,
                          HWND_DESKTOP,
                          NULL,
                          wcx.hInstance,
                          NULL);

    if (!hwnd)
        Error("CreateWindowEx() failed");

    ShowWindow(hwnd, SW_SHOW);
    SetFocus  (hwnd);

    hdc = GetDC(hwnd);
}

/*----------------------------------------------------------------------------*/
/* H�r nedanf�r �r publika funktioner som grafikbiblioteket exporterar till   */
/* andra program som vill anv�nda det. Allt ovan �r privata funktioner som    */
/* bara anv�nds internt av biblioteket.                                       */
/*----------------------------------------------------------------------------*/

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
    // TODO: B�r vi anropa Error() om man f�rs�ker k�ra InitGraphics() mer �n en
    //       g�ng?
    if (initialized)
        return;

    InitWindow     ();
    InitPixelFormat();
    InitOpenGL     ();

    initialized = TRUE;
    windowOpen  = TRUE;

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

    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);

    frameInterval = (int)(freq.QuadPart / fps);
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

    return windowOpen;
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
    glBegin   (GL_LINES);
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

    SwapBuffers(hdc);

    while (TRUE) {
        MSG msg;
        while (PeekMessage(&msg, hwnd, 0, 0, PM_REMOVE)) {
            TranslateMessage(&msg);
            DispatchMessage (&msg);
        }

        LARGE_INTEGER pc;
        QueryPerformanceCounter(&pc);
        pc.QuadPart -= lastUpdate.QuadPart;

        if (pc.QuadPart > frameInterval)
            break;
    }

    QueryPerformanceCounter(&lastUpdate);
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

    keyPressCB [c] = cb;
    keyPressArg[c] = arg;
}
