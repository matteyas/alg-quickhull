//------------------------------------------------
// INCLUDES
//------------------------------------------------

#include "common.h"
#include "graphics.h"
#include "io.h"
#include "math.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//------------------------------------------------
// CONSTANTS
//------------------------------------------------

/*--------------------------------------
 * Constant: BenchmarkNumIterations
 *
 * Description:
 *   Antal iterationer som ska k�ras i benchmarkl�ge.
 *------------------------------------*/
#define BenchmarkNumIterations 1000

#define LeftEdge -0.9f
#define RightEdge 0.9f
#define BottomEdge -0.9f
#define TopEdge 0.9f

//------------------------------------------------
// FUNCTIONS
//------------------------------------------------

/*--------------------------------------
 * Function: PrintIntroMessage()
 * Parameters:
 *
 * Description:
 *   Skriver ut introduktionsmeddelandet.
 *------------------------------------*/
static void PrintIntroMessage() {
    printf("Quickhull Demo v%s by %s\n\n", VERSION, AUTHOR);
}

static void RunBenchmark(int numPoints) {
    printf("Creating random points...\n");
    pointsetT ps = CreatePoints(numPoints);

    printf("Initializing hull...\n");
    hullT hull = InitHull(ps);

    printf("Done. Benchmark will now run.\n");
    printf("Benchmarking...\n");

    for (int i = 1; i <= BenchmarkNumIterations; i++) {
        RandomizePoints(ps);

        // H�r r�knar vi ut det konvexa h�ljet.
        BruteforceHull(ps, &hull);

        if ((i % 10)==0)
            printf("%d/%d iterations done...\n", i, BenchmarkNumIterations);
    }

    FreeHull(hull);
    FreePoints(ps);
}

static void RunDemo(int numPoints) {
    // Vi skapar upps�ttningen med punkter som ska visas p� sk�rmen.
    printf("Creating random points...\n");
    pointsetT ps = CreatePoints(numPoints);

    // Vi skapar en till upps�ttning punkter - lika m�nga som i ps - och
    // sparar dem i v. Vi anv�nder i sj�lva verket dessa som hastigheter f�r
    // punkterna i ps, f�r att kunna l�ta dem studsa omkring p� sk�rmen!
    printf("Randomizing velocities...\n");
    pointsetT vps = CreatePoints(ps.numPoints);

    // H�ljet skapar vi h�r f�r att slippa massa allokeringar inuti huvudloopen.
    printf("Initializing hull...\n");
    hullT hull = InitHull(ps);

    printf("Initializing graphics...\n");
    InitGraphics();
    SetFrameRate(60);
    printf("Done. Demo will now run.\n");

    pointsetT corners = CreatePoints(4);
    corners.points[0].x =  LeftEdge; corners.points[0].y =    TopEdge;
    corners.points[1].x = RightEdge; corners.points[1].y =    TopEdge;
    corners.points[2].x = RightEdge; corners.points[2].y = BottomEdge;
    corners.points[3].x =  LeftEdge; corners.points[3].y = BottomEdge;

    hullT edges = InitHull(corners);
    BruteforceHull(corners, &edges);

    while (WindowIsOpen()) {
        for (int i = 0; i < ps.numPoints; i++) {
            pointT *p = &ps .points[i],
                   *v = &vps.points[i];

            // Vi stegar punkterna fram�t linj�rt (d� vi inte har n�gon
            // acceleration).
            p->x += v->x / 100.0f;
            p->y += v->y / 100.0f;

            // H�r ser vi till att punkterna inte �ker f�r l�ngt ut mot kanterna,
            // utan att de i s�dana fall studsar mot kanten och �ker tillbaka.
                 if (p->x < LeftEdge  ) { p->x = LeftEdge  ; v->x *= -1.0f; }
            else if (p->x > RightEdge ) { p->x = RightEdge ; v->x *= -1.0f; }
                 if (p->y < BottomEdge) { p->y = BottomEdge; v->y *= -1.0f; }
            else if (p->y > TopEdge   ) { p->y = TopEdge   ; v->y *= -1.0f; }
        }

        // H�r r�knar vi ut det konvexa h�ljet.
        BruteforceHull(ps, &hull);

        // Dags att rita upp allting!
        ClearDisplay(0.992f, 0.964f, 0.890f);

        // F�rst ritar vi kanterna och h�ljet.
        SetColor(0.000f, 0.000f, 0.000f, 1.000f);
        DrawHull(edges);
        DrawHull(hull);

        // Sedan ritar vi prickarna (s� att de ligger ovanp� h�ljet).
        SetColor(0.866f, 0.000f, 0.000f, 1.000f);
        DrawPoints(ps);

        UpdateDisplay();
    }

    printf("Exiting...\n");

    FreePoints(vps);
    FreeHull(hull);
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
    string str       = GetStringFromUser();
    bool   benchmark = FALSE;
    if (str[0]=='Y' || str[0]=='y')
        benchmark = TRUE;
    free(str);

#ifndef _DEBUG
    // Vi slumpar inte "p� riktigt" i debugl�ge.
    srand((unsigned int)time(NULL));
#endif


    if (benchmark)
        RunBenchmark(numPoints);
    else
        RunDemo(numPoints);
}
