#ifndef NuNurbsH
#define NuNurbsH

#include <Nu/Nu.h>

#ifndef NuNurbsPRIVATE
typedef struct NuNurbsCurveStruct { void *object; } *NuNurbsCurve;
typedef struct NuNurbsSurfStruct { void *object; } *NuNurbsSurf;
#endif

NuInt NuPROC NuNurbsCurveEvaluate(NuNurbsCurve nurbsCurve, NuDouble parameter,
	NuDouble point[3]);

void NuPROC NuNurbsCurveFree(NuNurbsCurve nurbsCurve);

void NuPROC NuNurbsCurveGetDegree(NuNurbsCurve nurbsCurve, NuInt *degree);

NuInt NuPROC NuNurbsCurveGetKnot(NuNurbsCurve nurbsCurve, NuDouble knot[],
	NuInt knotCount);

void NuPROC NuNurbsCurveGetKnotCount(NuNurbsCurve nurbsCurve,
	NuInt *knotCount);

NuInt NuPROC NuNurbsCurveGetNormal(NuNurbsCurve nurbsCurve,
	NuDouble parameter, NuDouble point[3]);

NuInt NuPROC NuNurbsCurveGetNormals(NuNurbsCurve nurbsCurve,
	NuDouble normals[][3], NuInt normalsCount);

NuInt NuPROC NuNurbsCurveGetPoints(NuNurbsCurve nurbsCurve,
	NuDouble points[][4], NuInt pointsCount);

void NuPROC NuNurbsCurveGetPointsCount(NuNurbsCurve nurbsCurve,
	NuInt *pointsCount);

NuNurbsCurve NuPROC NuNurbsCurveNew(NuInt degree, NuDouble knot[],
	NuInt knotCount, NuDouble points[][4], NuInt pointsCount);

NuNurbsCurve NuPROC NuNurbsCurveNewArc(NuDouble center[3],
	NuDouble normal[3], NuDouble startPoint[3], NuDouble degrees);

NuNurbsCurve NuPROC NuNurbsCurveNewCircle(NuDouble center[3],
	NuDouble normal[3], NuDouble radius);

NuNurbsCurve NuPROC NuNurbsCurveNewInsertKnots(NuNurbsCurve nurbsCurve,
	NuDouble knot[], NuInt knotCount);

NuNurbsCurve NuPROC NuNurbsCurveNewPolyline(NuDouble points[][3],
	NuInt pointsCount);

NuNurbsCurve NuPROC NuNurbsCurveNewRefined(NuNurbsCurve nurbsCurve,
	NuInt factor);

NuInt NuPROC NuNurbsCurveSetKnot(NuNurbsCurve nurbsCurve, NuDouble knot[],
	NuInt knotCount);

NuInt NuPROC NuNurbsCurveSetPoints(NuNurbsCurve nurbsCurve,
	NuDouble points[][4], NuInt pointsCount);

NuInt NuPROC NuNurbsSurfEvaluate(NuNurbsSurf nurbsSurf,
	NuDouble uParameter, NuDouble vParameter, NuDouble point[3]);

void NuPROC NuNurbsSurfFree(NuNurbsSurf nurbsSurf);

void NuPROC NuNurbsSurfGetDegrees(NuNurbsSurf nurbsSurf, NuInt *uDegree,
	NuInt *vDegree);

NuInt NuPROC NuNurbsSurfGetKnots(NuNurbsSurf nurbsSurf,	NuDouble uKnot[],
	NuDouble vKnot[], NuInt uKnotCount, NuInt vKnotCount);

void NuPROC NuNurbsSurfGetKnotCounts(NuNurbsSurf nurbSurf,
	NuInt *uKnotCount, NuInt *vKnotCount);

NuInt NuPROC NuNurbsSurfGetNormal(NuNurbsSurf nurbsSurf, NuDouble uParameter,
	NuDouble vParameter, NuDouble point[3]);

NuInt NuPROC NuNurbSurfGetNormals(NuNurbsSurf nurbsSurf, NuDouble normals[][3],
	NuInt uPointsCount, NuInt vPointsCount);

NuInt NuPROC NuNurbsSurfGetPoints(NuNurbsSurf nurbsSurf,
	NuDouble points[][4], NuInt uPointsCount, NuInt vPointsCount);

void NuPROC NuNurbsSurfGetPointsCounts(NuNurbsSurf nurbsSurf,
	NuInt *uPointsCount, NuInt *vPointsCount);

NuNurbsSurf NuPROC NuNurbsSurfNew(NuInt uDegree, NuInt vDegree,
	NuDouble uKnot[], NuDouble vKnot[], NuInt uKnotCount, NuInt vKnotCount,
	NuDouble points[][4], NuInt uPointsCount, NuInt vPointsCount);

NuNurbsSurf NuPROC NuNurbsSurfNewCylinder(NuDouble center[3],
	NuDouble normal[3], NuDouble radius, NuDouble distance);

NuNurbsSurf NuPROC NuNurbsSurfNewCone(NuDouble center[3],
	NuDouble normal[3], NuDouble radius1, NuDouble radius2,
	NuDouble distance);

NuNurbsSurf NuPROC NuNurbsSurfNewExtruded(NuNurbsCurve nurbsCurve,
	NuDouble vec[3], NuDouble distance);

NuNurbsSurf NuPROC NuNurbsSurfNewInsertUKnots(NuNurbsSurf nurbsSurf,
	NuDouble knot[], NuInt knotCount);

NuNurbsSurf NuPROC NuNurbsSurfNewInsertVKnots(NuNurbsSurf nurbsSurf,
	NuDouble knot[], NuInt knotCount);

NuNurbsSurf NuPROC NuNurbsSurfNewMesh(NuDouble points[][3],
	NuInt uPointsCount,	NuInt vPointsCount);

NuNurbsSurf NuPROC NuNurbsSurfNewPlane(NuDouble point1[3],
	NuDouble point2[3], NuDouble point3[3]);

NuNurbsSurf NuPROC NuNurbsSurfNewRefined(NuNurbsSurf nurbsSurf, NuInt factor);

NuNurbsSurf NuPROC NuNurbsSurfNewRev(NuNurbsCurve nurbsCurve,
	NuDouble point1[3], NuDouble point2[3], NuDouble degrees);

NuNurbsSurf NuPROC NuNurbsSurfNewRuled(NuNurbsCurve nurbsCurve1,
	NuNurbsCurve nurbsCurve2);

NuNurbsSurf NuPROC NuNurbsSurfNewSphere(NuDouble center[3],
	NuDouble radius);

NuInt NuPROC NuNurbsSurfSetPoints(NuNurbsSurf nurbsSurf,
	NuDouble points[][4], NuInt uPointsCount, NuInt vPointsCount);

NuInt NuPROC NuNurbsSurfSetKnots(NuNurbsSurf nurbsSurf,	NuDouble uKnot[],
	NuDouble vKnot[], NuInt uKnotCount, NuInt vKnotCount);


#endif
