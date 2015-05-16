/*------------------------------------------------------------------------------
 * File: math.h
 * Created: May 15, 2015
 * Last changed: May 15, 2015
 *
 * Author(s): Philip Arvidsson (philip@philiparvidsson.com)
 *
 * Description:
 *   Matematikbibliotek f�r l�sning av olika matematiska problem.
 *
 * Changes:
 *
 *----------------------------------------------------------------------------*/

#ifndef Math_h
#define Math_h

/*------------------------------------------------
 * TYPES
 *----------------------------------------------*/

/*
 * Type: pointT
 *
 * Description:
 *   Representerar en punkt i euklidisk 2-dimensionell rymd.
 */
typedef struct {
    float x, y;
} pointT;

/*
 * Type: pointsetT
 *
 * Description:
 *   Representerar en upps�ttning punkter.
 */
typedef struct {
    pointT* points;
    int     numPoints;
} pointsetT;

/*
 * Type: lineT
 *
 * Description:
 *   Representerar en linje mellan tv� punkter.
 */
typedef struct {
    pointT *a, *b;
} lineT;

/*
 * Type: hullT
 *
 * Description:
 *   Representerar ett konvext h�lje.
 */
typedef struct {
    lineT *lines;
    int    numLines;
    int    maxLines;
} hullT;

/*--------------------------------------
 * Function: CreatePoints()
 * Parameters:
 *   n  Antalet punkter att skapa i upps�ttningen.
 *
 * Description:
 *   Skapar en upps�ttning punkter.
 *------------------------------------*/
pointsetT CreatePoints(int n);

/*--------------------------------------
 * Function: FreePoints()
 * Parameters:
 *   ps  Den upps�ttning punkter som ska avallokeras.
 *
 * Description:
 *   Avallokerar en upps�ttning punkter,
 *------------------------------------*/
void FreePoints(pointsetT ps);

/*--------------------------------------
 * Function: InitHull()
 * Parameters:
 *   ps  Den upps�ttning punkter som ska anv�ndas till att initiera h�ljet.
 *
 * Description:
 *   Initierar ett nytt h�lje.
 *------------------------------------*/
hullT InitHull(pointsetT ps);

/*--------------------------------------
 * Function: FreeHull()
 * Parameters:
 *   hull Det h�lje som ska avallokeras.
 *
 * Description:
 *   Avallokerar ett h�lje.
 *------------------------------------*/
void FreeHull(hullT hull);

/*--------------------------------------
 * Function: BruteforceHull()
 * Parameters:
 *   ps    Punktupps�ttningen f�r vilken ett h�lje ska genereras.
 *   hull  En pekare till h�ljet.
 *
 * Description:
 *   Genererar att konvext h�lje f�r punktupps�ttningen genom utt�mmande
 *   s�kning.
 *------------------------------------*/
void BruteforceHull(pointsetT ps, hullT *hull);

/*--------------------------------------
 * Function: RandomizePoints()
 * Parameters:
 *   ps    Punktupps�ttningen vars positioner ska slumpas.
 *
 * Description:
 *   Slumpar positionerna f�r punkterna i den specificerade punktupps�ttningen.
 *------------------------------------*/
void RandomizePoints(pointsetT ps);

/*--------------------------------------
 * Function: Reflect()
 * Parameters:
 *   d  Riktningsvektorn.
 *   n  Normalvektorn som riktningsvektorn ska reflekteras mot.
 *
 * Description:
 *   Reflekterar en riktningsvektor mot en normalvektor och returnerar
 *   resultatet.
 *------------------------------------*/
pointT Reflect(pointT d, pointT n);

#endif // Math_h
