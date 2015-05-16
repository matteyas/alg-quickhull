/*------------------------------------------------------------------------------
 * File: alg-quickhull.c
 * Created: May 15, 2015
 * Last changed: May 15, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Huvudprogram f�r l�sning av det konvexa h�lja med hj�lpa av algoritmen
 *   quickhull.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

/*------------------------------------------------
 * INCLUDES
 *----------------------------------------------*/

#include "common.h"
#include "graphics.h"
#include "io.h"
#include "math.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*------------------------------------------------
 * CONSTANTS
 *----------------------------------------------*/

/*--------------------------------------
 * Constant: BenchmarkNumIterations
 *
 * Description:
 *   Antal iterationer som ska k�ras i benchmarkl�ge.
 *------------------------------------*/
#define BenchmarkNumIterations 1000

/*--------------------------------------
 * Constant: DemoSpeed
 *
 * Description:
 *   Hastighetsmultipel f�r demons k�rningshastighet vad g�ller simulering av
 *   punkterna. 0.5 = halv hastighet, 1.0 = normal, 2.0 = dubbel, etc.
 *------------------------------------*/
#define DemoSpeed 0.5f

/*--------------------------------------
 * Constant: FrameRate
 *
 * Description:
 *   Antal bildrutor per sekund.
 *------------------------------------*/
#define FrameRate 60

/*--------------------------------------
 * Constant: StepSize
 *
 * Description:
 *   Tidssteget vi simulerar i fysiksimuleringen varje bildruta.
 *------------------------------------*/
#define StepSize 1.0f / 120.0f

/*--------------------------------------
 * Constant: BottomEdge
 *
 * Description:
 *   Nedre kantens y-koordinat.
 *------------------------------------*/
#define BottomEdge -0.9f

/*--------------------------------------
 * Constant: LeftEdge
 *
 * Description:
 *   V�nsterkantens x-koordinat.
 *------------------------------------*/
#define LeftEdge -0.9f

/*--------------------------------------
 * Constant: RightEdge
 *
 * Description:
 *   H�gerkantens x-koordinat.
 *------------------------------------*/
#define RightEdge 0.9f

/*--------------------------------------
 * Constant: TopEdge
 *
 * Description:
 *   �vre kantens y-koordinat.
 *------------------------------------*/
#define TopEdge 0.9f

/*--------------------------------------
 * Constant: SpringCoefficient
 *
 * Description:
 *   Fj�dringskoefficienten. H�gre v�rde ger "stelare" gummiband runt prickarna.
 *------------------------------------*/
#define SpringCoefficient 30.0f

/*--------------------------------------
 * Constant: SpringDamping
 *
 * Description:
 *   Mjukar sammandragningen och expansionen av gummibandet.
 *------------------------------------*/
#define SpringDamping 0.5f

/*--------------------------------------
 * Constant: SpringLength
 *
 * Description:
 *   L�ngden som segmenten i gummibandet efterstr�var.
 *------------------------------------*/
#define SpringLength 0.1f

/*------------------------------------------------
 * GLOBALS
 *----------------------------------------------*/

bool blackHole   = FALSE;
bool damping     = FALSE;
bool drawHull    = TRUE;
bool drawPoints  = TRUE;
bool gravity     = FALSE;
bool rubberBand  = FALSE;
bool slopedFloor = FALSE;
bool wind        = FALSE;

/*------------------------------------------------
 * FUNCTIONS
 *----------------------------------------------*/

static void RandomizePS(void *arg) {
    pointsetT *ps = (pointsetT *)arg;
    RandomizePoints(*ps);
}

static void ToggleBlackHole(void *arg) {
    blackHole = !blackHole;

    if (blackHole) printf(":: Black hole enabled.\n");
    else           printf(":: Black hole disabled.\n");
}

static void ToggleDamping(void *arg) {
    damping = !damping;

    if (damping) printf(":: Damping enabled.\n");
    else         printf(":: Damping disabled.\n");
}

static void ToggleGravity(void *arg) {
    gravity = !gravity;

    if (gravity) printf(":: Gravity enabled.\n");
    else         printf(":: Gravity disabled.\n");
}

static void ToggleHull(void *arg) {
    drawHull = !drawHull;

    if (drawHull) printf(":: Hull rendering enabled.\n");
    else          printf(":: Hull rendering disabled.\n");
}

static void TogglePoints(void *arg) {
    drawPoints = !drawPoints;

    if (drawPoints) printf(":: Points rendering enabled.\n");
    else            printf(":: Points rendering disabled.\n");
}

static void ToggleRubberBand(void *arg) {
    rubberBand = !rubberBand;

    if (rubberBand) printf(":: Rubber band enabled.\n");
    else            printf(":: Rubber band disabled.\n");
}

static void ToggleWind(void *arg) {
    wind = !wind;

    if (wind) printf(":: Wind enabled.\n");
    else      printf(":: Wind disabled.\n");
}

static void ToggleSlopedFloor(void *arg) {
    slopedFloor = !slopedFloor;

    if (slopedFloor) printf(":: Sloped floor enabled.\n");
    else             printf(":: Sloped floor disabled.\n");
}

/*--------------------------------------
 * Function: PrintIntroMessage()
 * Parameters:
 *
 * Description:
 *   Skriver ut introduktionsmeddelandet.
 *------------------------------------*/
static void PrintIntroMessage() {
    printf("Quickhull Demo v%s by %s\n"
           "\n"
           " - The most awesome Quickhull demo ever made!\n\n",
           ProgramVersion, ProgramAuthors);
}

static void PrintInstructions() {
    printf("\nINSTRUCTIONS:\n\n"
           "KEYS\n------------\n"
           "Key    Effect\n\n"
           "  b    Toggle black hole (gravity towards center).\n"
           "  d    Toggle damping.\n"
           "  g    Toggle gravity.\n"
           "  h    Toggle hull rendering.\n"
           "  p    Toggle point rendering.\n"
           "  r    Toggle rubber band mode (hull becomes a rubber band).\n"
           "  s    Toggle sloped floor.\n"
           "  w    Toggle wind.\n"
           "  x    Randomize point velocities.\n"
           "  z    Randomize point positions.\n"
           "------------\n\n");
}

/*--------------------------------------
 * Function: RunBenchmark()
 * Parameters:
 *   numPoints  Antal punkter att anv�nda vid utr�kning av det konvexa
 *              h�ljet.
 *
 * Description:
 *   K�r programmet i benchmark-l�ge.
 *------------------------------------*/
static void RunBenchmark(int numPoints) {
    printf("Creating random points...");
    pointsetT ps = CreatePoints(numPoints);
    printf(" done.\n");

    printf("Initializing hull...");
    hullT hull = InitHull(ps);
    printf(" done.\n");

    printf("Benchmark will now run.\n");
    printf("Benchmarking...\n");

    for (int i = 1; i <= BenchmarkNumIterations; i++) {
        RandomizePoints(ps);

        // H�r r�knar vi ut det konvexa h�ljet.
        BruteforceHull(ps, &hull);

        if ((i % 10)==0)
            printf("%d/%d iterations done...\n", i, BenchmarkNumIterations);
    }

    FreeHull  (hull);
    FreePoints(ps);
}

static void UpdatePoints(pointsetT aps, pointsetT ps, pointsetT vps, hullT hull,
                         float dt)
{
    if (rubberBand) {
        // H�r nedan till�mpar vi Hookes lag p� segmenten i h�ljet, s� att
        // h�ljet drar sig samman runt punkterna likt ett gummiband.
        //
        // Se http://en.wikipedia.org/wiki/Hooke%27s_law f�r mer information.
        for (int k = 0; k < hull.numLines; k++) {
            lineT line = hull.lines[k];

            // Nedan anv�nder vi lite ful pekararitmetik f�r att lista
            // ut punkternas index, s� att vi kan arbeta med deras
            // hastigheter och accelerationer.
            int i = line.a - ps.points;
            int j = line.b - ps.points;

            pointT dp = { ps.points[j].x - ps.points[i].x,
                          ps.points[j].y - ps.points[i].y };

            pointT dv = { vps.points[j].x - vps.points[i].x,
                          vps.points[j].y - vps.points[i].y };

            float r = (float)sqrt(dp.x*dp.x + dp.y*dp.y);

            // Punkterna ligger f�r n�ra varandra f�r att vi ska kunna trycka
            // is�r dem utan avrundningsproblem.
            if (r < 0.000001f)
                continue;

            dp.x /= r;
            dp.y /= r;

            r -= SpringLength;

            pointT f = { dp.x*r*SpringCoefficient + dv.x*SpringDamping,
                         dp.y*r*SpringCoefficient + dv.y*SpringDamping };

            aps.points[i].x += f.x;
            aps.points[i].y += f.y;
            aps.points[j].x -= f.x;
            aps.points[j].y -= f.y;
        }
    }

    for (int i = 0; i < ps.numPoints; i++) {
        pointT *a = &aps.points[i],
               *p = &ps .points[i],
               *v = &vps.points[i];

        //-----------------------------------
        // H�r till�mpar vi milj�m�ssiga effekter p� punkterna.

        if (blackHole) {
            a->x -= p->x * 10.0f;
            a->y -= p->y * 10.0f;
        }

        if (gravity)
            a->y -= 9.81f;

        if (damping) {
            a->x -= v->x;
            a->y -= v->y;
        }

        if (wind)
            a->x += 7.0f;

        //-----------------------------------

        // Vi stegar punkterna fram�t linj�rt med Euler-integrering.
        v->x += a->x * dt; v->y += a->y * dt;
        p->x += v->x * dt; p->y += v->y * dt;

        // Efter att vi till�mpat accelerationerna s� nollst�ller vi dem.
        a->x = a->y = 0.0f;

        // H�r ser vi till att punkterna inte �ker f�r l�ngt ut mot kanterna,
        // utan att de i s�dana fall studsar mot kanten och �ker tillbaka.
             if (p->x <  LeftEdge) { p->x =  LeftEdge; v->x *= -1.0f; }
        else if (p->x > RightEdge) { p->x = RightEdge; v->x *= -1.0f; }
             if (p->y >   TopEdge) { p->y =   TopEdge; v->y *= -1.0f; }

        // Ett lutande golv beh�ver lite specialbehandling.
        if (slopedFloor) {
            float x     = (p->x-LeftEdge) / (RightEdge-LeftEdge);
            float floor = -0.9f * (1.0f-x) - 0.3f * x;

            if (p->y < floor) {
                p->y  = floor;

                // Vektorn n �r golvets perp-vektor normaliserad. Vi r�knar ut
                // den s�h�r:
                //
                //   g = (0.9 - (-0.9), (-0.3) - (-0.9))
                //     = (1.8, 0.6)
                //
                //   g_perp = (-0.6, 1.8)
                //
                //   n^ = g_perp / ||g_perp|| 
                //      = (-0.6, 1.8) / ||(-0.6, 1.8)||
                //      = (-0.6, 1.8) / 1.89737
                //      = (-0.316228, 0.948683)
                //
                pointT n = { -0.316228f, 0.948683f };
                *v = Reflect(*v, n);
            }
        }
        else if (p->y < BottomEdge) {
            p->y  = BottomEdge;
            v->y *= -1.0f;
        }
    }
}

/*--------------------------------------
 * Function: RunDemo()
 * Parameters:
 *   numPoints  Antal punkter att anv�nda vid utr�kning av det konvexa h�ljet.
 *
 * Description:
 *   K�r programmet i demo-l�ge.
 *------------------------------------*/
static void RunDemo(int numPoints) {
    // Vi skapar upps�ttningen med punkter som ska visas p� sk�rmen.
    printf("Creating random points...");
    pointsetT ps = CreatePoints(numPoints);
    RandomizePoints(ps);
    printf(" done.\n");

    // Vi skapar en till upps�ttning punkter - lika m�nga som i ps - och
    // sparar dem i v. Vi anv�nder i sj�lva verket dessa som hastigheter f�r
    // punkterna i ps, f�r att kunna l�ta dem studsa omkring p� sk�rmen!
    printf("Randomizing velocities...");
    pointsetT vps = CreatePoints(ps.numPoints);
    RandomizePoints(vps);
    printf(" done.\n");

    // Den tredje upps�ttningen anv�nder vi till accelerationer.
    printf("Creating acceleration set...");
    pointsetT aps = CreatePoints(ps.numPoints);
    for (int i = 0; i < aps.numPoints; i++) {
        aps.points[i].x = 0.0f;
        aps.points[i].y = 0.0f;
    }
    printf(" done.\n");

    // H�ljet skapar vi h�r f�r att slippa massa allokeringar inuti huvudloopen.
    printf("Initializing hull...");
    hullT hull = InitHull(ps);
    printf(" done.\n");

    // H�r skapar vi de fyra h�rnen och r�knar ut h�ljet(!) f�r dem, f�r att f�
    // en fin ruta som visar var kanterna g�r.
    printf("Initializing world hull...");
    pointsetT corners = CreatePoints(4);
    hullT     edges   = InitHull    (corners);

    corners.points[0].x =  LeftEdge; corners.points[0].y =    TopEdge;
    corners.points[1].x = RightEdge; corners.points[1].y =    TopEdge;
    corners.points[2].x = RightEdge; corners.points[2].y = BottomEdge;
    corners.points[3].x =  LeftEdge; corners.points[3].y = BottomEdge;

    BruteforceHull(corners, &edges);

    printf(" done.\n");

    PrintInstructions();

    printf("Press ENTER to start...\n");
    getchar();

    printf("Initializing graphics...");
    InitGraphics();
    SetFrameRate(FrameRate);

    // H�r st�ller vi in lite roliga knappar s� att demon blir lite mer
    // interaktiv och rolig.
    OnKeyPress('b', ToggleBlackHole  , NULL);
    OnKeyPress('d', ToggleDamping    , NULL);
    OnKeyPress('f', ToggleSlopedFloor, NULL);
    OnKeyPress('g', ToggleGravity    , NULL);
    OnKeyPress('h', ToggleHull       , NULL);
    OnKeyPress('p', TogglePoints     , NULL);
    OnKeyPress('r', ToggleRubberBand , NULL);
    OnKeyPress('w', ToggleWind       , NULL);
    OnKeyPress('x', RandomizePS      , &vps);
    OnKeyPress('z', RandomizePS      , &ps);

    printf(" done.\n");
    printf("Enjoy! :-)\n\n");

    float dt = 0.0f;
    while (WindowIsOpen()) {
        dt += DemoSpeed / FrameRate;
        while (dt > StepSize) {
            UpdatePoints(aps, ps, vps, hull, StepSize);
            dt -= StepSize;
        }

        // Dags att rita upp allting!
        ClearCanvas(0.992f, 0.964f, 0.890f);

        // F�rst ritar vi kanterna och h�ljet.
        SetColor(0.000f, 0.000f, 0.000f, 1.000f);
        DrawHull(edges);

        if (drawHull) {
            BruteforceHull(ps, &hull);
            DrawHull(hull);
        }

        // Sedan milj�grejer.
        if (slopedFloor) {
            pointT a = {  LeftEdge, -0.9f },
                   b = { RightEdge, -0.3f };

            lineT floor = { &a, &b };

            DrawLine(floor);
        }

        if (drawPoints) {
            // Sedan prickarna (s� att de ligger ovanp� h�ljet).
            SetColor  (0.866f, 0.000f, 0.000f, 1.000f);
            DrawPoints(ps);
        }

        // Fram med allt p� sk�rmen!
        UpdateDisplay();
    }

    printf("Exiting...\n");

    FreeHull  (edges);
    FreePoints(corners);

    FreeHull  (hull);
    FreePoints(vps);
    FreePoints(ps);
}

/*--------------------------------------
 * Function: main()
 * Parameters:
 *
 * Description:
 *   Programmets huvudfunktion.
 *------------------------------------*/
main() {
    PrintIntroMessage();

    printf("Number of points to use? ");
    int numPoints = GetIntFromUser();
    if (numPoints < 1)    numPoints = 1;
    if (numPoints > 1000) numPoints = 1000;

    printf("Do you want to run the benchmark? (y/N) ");
    bool benchmark = GetBoolFromUser(FALSE);

    printf("\n");

#ifndef _DEBUG
    // Vi slumpar inte "p� riktigt" i debugl�ge.
    srand((unsigned int)time(NULL));
#endif

    if (benchmark) RunBenchmark(numPoints);
    else           RunDemo     (numPoints);
}
