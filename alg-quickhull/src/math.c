//------------------------------------------------
// INCLUDES
//------------------------------------------------

#include "common.h"
#include "math.h"

#include <stdlib.h>
#include <time.h>

//------------------------------------------------
// FUNCTIONS
//------------------------------------------------

pointsetT CreatePoints(int n) {
    pointsetT ps;

    ps.points    = malloc(sizeof(pointT) * n);
    ps.numPoints = n;

    RandomizePoints(ps);

    return ps;
}

void FreePoints(pointsetT ps) {
    free(ps.points);
    ps.numPoints = 0;
}

hullT InitHull(pointsetT ps) {
    hullT hull;

    hull.lines    = malloc(sizeof(lineT) * ps.numPoints);
    hull.numLines = 0;
    hull.maxLines = ps.numPoints;

    return hull;
}

void FreeHull(hullT hull) {
    free(hull.lines);
    hull.maxLines = 0;
    hull.numLines = 0;
}

void BruteforceHull(pointsetT ps, hullT *hull) {
    // For-looparna med i och j anv�nds f�r att konstruera alla t�nkbara
    // kombinationer av par bland punkterna.

    hull->numLines = 0;
    for (int i = 0; i < ps.numPoints; i++) {
        for (int j = (i+1); j < (i+ps.numPoints); j++) {
            pointT *a, *b, *c;

            // For-loopen med k anv�nds f�r att kolla att alla punkter ligger p�
            // r�tt sida av linjen mellan punkt i och j.

            bool outside = FALSE;
            for (int k = 0; k < ps.numPoints; k++) {
                // Linjen a---b �r (potentiellt) ett segment i h�ljet. Vi avg�r
                // det genom att se vilken sida om linjen som punkten c ligger
                // p�.
                a = &ps.points[i];
                b = &ps.points[j % ps.numPoints];
                c = &ps.points[k];

                // Nedan avg�r vi vilken sida om linjen punkten ligger p�.
                // d<0.0 : Punkten �r utanf�r.
                // d=0.0 : Punkten ligger p� linjen.
                // d>0.0 : Punkten �r innanf�r.
                float d = (b->x - a->x) * (c->y - a->y)
                        - (b->y - a->y) * (c->x - a->x);

                if (d < 0.0f) {
                    outside = TRUE;
                    break;
                }
            }

            // Om alla punkter visade sig vara innanf�r linjen, s� beh�ller vi
            // den som ett segment i det konvexa h�ljet.
            if (!outside) {
                hull->lines[hull->numLines].a = a;
                hull->lines[hull->numLines].b = b;
                hull->numLines++;

                if (hull->numLines >= hull->maxLines)
                    return;
            }
        }
    }

    // -
}

void Quickhull(pointsetT ps, hullT *hull) {
    pointT *a = &ps.points[0], *b = a;

    for (int i = 1; i < ps.numPoints; i++) {
        pointT *p = &ps.points[i];

        if (p->x < a->x) a = p;
        if (p->x > b->x) b = p;
    }
}

void RandomizePoints(pointsetT ps) {
    for (int i = 0; i < ps.numPoints; i++) {
        ps.points[i].x = 1.6f * ((float)rand() / (float)RAND_MAX) - 0.8f;
        ps.points[i].y = 1.6f * ((float)rand() / (float)RAND_MAX) - 0.8f;
    }
}
