/* Some NURBS code from 1993. The only really interesting piece of code is
   the insert knots code that uses Boehms multiple knot insertion algorithm
   as it is a very fast and efficient way to insert multiple knots into
   NURBS curves and surfaces

   Eric Swildens 1993 (updated 2016)
*/

#include <Nu/Nu.h>

typedef struct NuNurbsCurveStruct
	{
	NuMem mem;
	NuInt degree;			/* degree */
	NuDouble *knot;			/* knot vector */
	NuInt knotCount;		/* num of knots: m=n+p+1*/
	NuDouble (*points)[4];	/* control points: x,y,z,w */
	NuInt pointsCount;		/* num of control points */
	} *NuNurbsCurve;

typedef struct NuNurbsSurfStruct
	{
	NuMem mem;
	NuInt uDegree;
	NuInt vDegree;			/* degree */
	NuDouble *uKnot;
	NuDouble *vKnot;		/* u,v knot vectors */
	NuInt uKnotCount;
	NuInt vKnotCount;		/* num of knots: m=n+p+1*/
	NuDouble (*points)[4];	/* m x n matrix of control points: x,y,z,w */
	NuInt uPointsCount;
	NuInt vPointsCount;		/* num of control points */
	} *NuNurbsSurf;

static void NuNurbsCurveInsertKnot(NuNurbsCurve newNurbsCurve,
	NuNurbsCurve nurbsCurve, NuDouble knot[], NuInt i, NuInt j,
	NuInt oldKnotCount);

static void NuNurbsEvaluate(NuDouble ePoints[][4], NuDouble knot[],
	NuDouble parameter, NuInt degree, NuInt knotPos, NuDouble point[3]);

static void NuNurbsSurfInsertUKnot(NuNurbsSurf newNurbsSurf,
	NuNurbsSurf nurbsSurf, NuDouble knot[], NuInt i, NuInt j,
	NuInt oldUKnotCount);

static void NuNurbsSurfInsertVKnot(NuNurbsSurf newNurbsSurf,
	NuNurbsSurf nurbsSurf, NuDouble knot[], NuInt i, NuInt j,
	NuInt oldVKnotCount);

#define NuNurbsPRIVATE 1
#include <Nu/Nurbs.h>

/*
NOTES:
	for evaluation:
	((t >= n->knot[n->uPointsCount-1]) || (t < n->knot[n->p])) should be false
	parameter [0..1)

	for curves:
	n > p

	for surfs:
	m > p and n > q

	weights are (0..1000]
*/

NuInt NuPROC NuNurbsCurveEvaluate(NuNurbsCurve nurbsCurve,
	NuDouble parameter, NuDouble point[3])
/* parameter [0..1) */
/* 	An Introduction to Splines for Use in Computer Graphics
	and Geometric Modeling, pg. 392
*/	{
	NuDouble (*ePoints)[4];		/* evaluation points */
	NuInt i, j, knotPos, size;

	knotPos = nurbsCurve->degree;
	while (parameter >= nurbsCurve->knot[knotPos++])
		;
	knotPos--;
	size = sizeof(NuDouble) * (nurbsCurve->degree + 1) * 4;
	ePoints = (NuDouble (*)[4])NuMem_allocate(nurbsCurve->mem, size);
	i = knotPos - 1;
	for (j = 0; j < nurbsCurve->degree + 1; j++)
		{
		ePoints[j][3] = nurbsCurve->points[i][3];
		ePoints[j][0] = nurbsCurve->points[i][0] * ePoints[j][3];
		ePoints[j][1] = nurbsCurve->points[i][1] * ePoints[j][3];
		ePoints[j][2] = nurbsCurve->points[i][2] * ePoints[j][3];
		i--;
		}
	NuNurbsEvaluate(ePoints, nurbsCurve->knot, parameter,
		nurbsCurve->degree, knotPos, point);
	NuMem_deallocate(nurbsCurve->mem, (char *)ePoints);
	return 0;
	}

void NuPROC NuNurbsCurveFree(NuNurbsCurve nurbsCurve)
	{
	NuMem_free(nurbsCurve->mem);
	}

void NuPROC NuNurbsCurveGetDegree(NuNurbsCurve nurbsCurve, NuInt *degree)
	{
	*degree = nurbsCurve->degree;
	}

NuInt NuPROC NuNurbsCurveGetKnot(NuNurbsCurve nurbsCurve, NuDouble knot[],
	NuInt knotCount)
	{
	NuInt i;

	if (knotCount != nurbsCurve->knotCount)
		return -1;
	for (i = 0; i <knotCount; i++)
		knot[i] = nurbsCurve->knot[i];
	return 0;
	}

void NuPROC NuNurbsCurveGetKnotCount(NuNurbsCurve nurbsCurve,
	NuInt *knotCount)
	{
	*knotCount = nurbsCurve->knotCount;
	}

NuInt NuPROC NuNurbsCurveGetNormal(NuNurbsCurve nurbsCurve,
        NuDouble parameter, NuDouble point[3])
	{
	return -1;
	}

NuInt NuPROC NuNurbsCurveGetNormals(NuNurbsCurve nurbsCurve,
        NuDouble normals[][3], NuInt normalsCount)
/*	An Introduction to Splines for use in Computer Graphics
	and Geometric Modelsing, p 397
*/	{
    return -1;

#ifdef NEVER
	NuDouble (*c)[4],temp;
	NuInt i,j,k,r,s,s1,size;

	k = 0;
	while (t >= n->knot[k++])
		;
	k--;
	size = sizeof(NuDouble)*(n->p+1)*4;
	c = (NuDouble (*)[4])NuMem_allocate(n->mem,size);
	for (s=0; s<n->p+1; s++)
		{
		c[s][3] = n->points[k-s-1][3];
		c[s][0] = n->points[k-s-1][0];
		c[s][1] = n->points[k-s-1][1];
		c[s][2] = n->points[k-s-1][2];
		}
	for (r=1; r<d+1; r++)
		{
		i = k - 1;
		j = n->p - r + 1;
		s1 = 1;
		for (s=0; s<n->p-r+1; s++)
			{
			temp = n->knot[i+j] - n->knot[i];
			if (temp > 0.0)
				{
				c[s][0] = (j*(c[s][0] - c[s1][0]))/temp;
				c[s][1] = (j*(c[s][1] - c[s1][1]))/temp;
				c[s][2] = (j*(c[s][2] - c[s1][2]))/temp;
				c[s][3] = (j*(c[s][3] - c[s1][3]))/temp;
				}
			else
				{
				c[s][0] = 0.0;
				c[s][1] = 0.0;
				c[s][2] = 0.0;
				c[s][3] = 0.0;
				}
			s1++;
			i--;
			}
		}
	NuNurbs_eval(c,n->knot,t,n->p-d,k,p);
	NuMem_deallocate(n->mem,(char *)c);
	return(1);
#endif
	}

NuInt NuPROC NuNurbsCurveGetPoints(NuNurbsCurve nurbsCurve,
	NuDouble points[][4], NuInt pointsCount)
	{
	NuInt i;

	if (pointsCount != nurbsCurve->pointsCount)
		return -1;
	for (i = 0; i < pointsCount; i++)
		{
		points[i][0] = nurbsCurve->points[i][0];
		points[i][1] = nurbsCurve->points[i][1];
		points[i][2] = nurbsCurve->points[i][2];
		points[i][3] = nurbsCurve->points[i][3];
		}
	return 0;
	}

void NuPROC NuNurbsCurveGetPointsCount(NuNurbsCurve nurbsCurve,
	NuInt *pointsCount)
	{
	*pointsCount = nurbsCurve->pointsCount;
	}

static void NuNurbsCurveInsertKnot(NuNurbsCurve newNurbsCurve,
	NuNurbsCurve nurbsCurve, NuDouble knot[], NuInt i, NuInt j,
	NuInt oldKnotCount)
	{
	NuInt l, l1, l2;
	NuDouble B, Bw1, Bw2, Bwsum;

	for (l = 0; l < nurbsCurve->degree; l++)
		{
		B = newNurbsCurve->knot[j + l + 1] - knot[i];
		l1 = j - nurbsCurve->degree + l;
		l2 = l1 + 1;
		if (B <= 0.0)
			{
			newNurbsCurve->points[l1][0] = newNurbsCurve->points[l2][0];
			newNurbsCurve->points[l1][1] = newNurbsCurve->points[l2][1];
			newNurbsCurve->points[l1][2] = newNurbsCurve->points[l2][2];
			newNurbsCurve->points[l1][3] = newNurbsCurve->points[l2][3];
			}
		else
			{
			B /= (newNurbsCurve->knot[j + l + 1] -
				nurbsCurve->knot[oldKnotCount - nurbsCurve->degree + l + 1]);
			Bw1 = B * newNurbsCurve->points[l1][3];
			Bw2 = (1.0 - B) * newNurbsCurve->points[l2][3];
			Bwsum = Bw1 + Bw2;
			newNurbsCurve->points[l1][0] = (Bw1 * newNurbsCurve->points[l1][0] +
				Bw2 * newNurbsCurve->points[l2][0]) / Bwsum;
			newNurbsCurve->points[l1][1] = (Bw1 * newNurbsCurve->points[l1][1] +
				Bw2 * newNurbsCurve->points[l2][1]) / Bwsum;
			newNurbsCurve->points[l1][2] = (Bw1 * newNurbsCurve->points[l1][2] +
				Bw2 * newNurbsCurve->points[l2][2]) / Bwsum;
			newNurbsCurve->points[l1][3] = Bwsum;
			}
		}
	}

NuNurbsCurve NuPROC NuNurbsCurveNew(NuInt degree, NuDouble knot[],
	NuInt knotCount, NuDouble points[][4], NuInt pointsCount)
	{
	NuMem mem;
	NuNurbsCurve nurbsCurve;
	NuInt i, size;

	if (knotCount != pointsCount + degree + 1)
		return 0;
	mem = NuMem_new("NuNurbsCurve");
	size = sizeof(struct NuNurbsCurveStruct);
	nurbsCurve = (NuNurbsCurve)NuMem_allocate(mem, size);
	nurbsCurve->mem = mem;
	nurbsCurve->degree = degree;
	nurbsCurve->knotCount = knotCount;
	size = sizeof(NuDouble) * knotCount;
	nurbsCurve->knot = (NuDouble *)NuMem_allocate(mem, size);
	for (i = 0; i < knotCount; i++)
		nurbsCurve->knot[i] = knot[i];
	nurbsCurve->pointsCount = pointsCount;
	size = sizeof(NuDouble) * pointsCount * 4;
	nurbsCurve->points = (NuDouble (*)[4])NuMem_allocate(mem, size);
	for (i = 0; i < pointsCount; i++)
		{
		nurbsCurve->points[i][0] = points[i][0];
		nurbsCurve->points[i][1] = points[i][1];
		nurbsCurve->points[i][2] = points[i][2];
		nurbsCurve->points[i][3] = points[i][3];
		}
	return(nurbsCurve);
	}

NuNurbsCurve NuPROC NuNurbsCurveNewArc(NuDouble center[3],
	NuDouble normal[3], NuDouble startPoint[3], NuDouble degrees)
	{
	return 0;

#ifdef NEVER
	r = 4.0;
	radius = NuDouble_sqrt(r);
	points[0][0] = p1[0];
	points[0][1] = p1[1];
	points[0][2] = p1[2];
	points[0][3] = 1.0;
	radius = NuDouble_sqrt ((p1[0] - center[0]) * (p1[0] - center[0]) +
		(p1[1] - center[1]) * (p1[1] - center[1]) +
		(p1[2] - center[2]) * (p1[2] - center[2]));
	r = NuDouble_sqrt ((p2[0] - center[0]) * (p2[0] - center[0]) +
		(p2[1] - center[1]) * (p2[1] - center[1]) +
		(p2[2] - center[2]) * (p2[2] - center[2]));
	midp[0] = (p1[0] + p2[0]) / 2.0;
	midp[1] = (p1[1] + p2[1]) / 2.0;
	midp[2] = (p1[2] + p2[2]) / 2.0;
	midp_len = (midp[0] - center[0]) * (midp[0] - center[0]) +
		(midp[1] - center[1]) * (midp[1] - center[1]) +
		(midp[2] - center[2]) * (midp[2] - center[2]);
	midp_len = NuDouble_sqrt(midp_len);
	l1 = radius - midp_len;
	l2 = NuDouble_sqrt((radius * radius) - (midp_len * midp_len));
	l3 = (2.0 * l1) / (1.0 - ((l1 / l2) * (l1 / l2)));
	fac = (l3 + midp_len) / midp_len;
	points[1][0] = (midp[0] - center[0]) * fac + center[0];
	points[1][1] = (midp[1] - center[1]) * fac + center[1];
	points[1][2] = (midp[2] - center[2]) * fac + center[2];
	points[1][3] = l2 / NuDouble_sqrt((l3 * l3) + (l2 * l2));
	points[2][0] = p2[0];
	points[2][1] = p2[1];
	points[2][2] = p2[2];
	points[2][3] = 1.0;
	nc = NuNurbsCurve_new((NuInt)2, knot, points, (NuInt)3);
	return(nc);
#endif
	}

NuNurbsCurve NuPROC NuNurbsCurveNewCircle(NuDouble center[3],
	NuDouble normal[3], NuDouble radius)
	{
	return 0;

#ifdef NEVER
	NuNurbsCurve nc;
	NuDouble points[7][4];
	static NuDouble knot[10]=
		{0.0, 0.0, 0.0, 1./4., 1./2., 1./2., 3./4., 1.0, 1.0, 1.0};
	static NuDouble points_scale[7][4] =
		{
		{ 1.0,   0.0,   0.0,   1.0},
		{ 1.0,   1.0,   0.0,    .5},
		{-1.0,   1.0,   0.0,    .5},
		{-1.0,   0.0,   0.0,   1.0},
		{-1.0,  -1.0,   0.0,    .5},
		{ 1.0,  -1.0,   0.0,    .5},
		{ 1.0,   0.0,   0.0,   1.0}
		};
	NuInt loop;

	for (loop=0; loop<7; loop++)
		{
		points[loop][0] = points_scale[loop][0] * radius;
		points[loop][1] = points_scale[loop][1] * radius;
		points[loop][2] = points_scale[loop][2] * radius;
		points[loop][3] = points_scale[loop][3];
		}
	nc = NuNurbsCurve_new(2, knot, points, 7);
	return(nc);
#endif
	}

NuNurbsCurve NuPROC NuNurbsCurveNewInsertKnots(NuNurbsCurve nurbsCurve,
	NuDouble knot[], NuInt knotCount)
/*	Boehms algorithm.  See "Inserting new knots into B-Spline Curves",
	Computer Aided Design Vol 12 No 4 (1980), pg 199-202.  This
	algorithm is the modified Boehm algorithm from "The Insertion
	Algorithm", Computer Aided Design Vol 17 (1985), pg 58-59.
*/	{
	NuMem mem;
	NuNurbsCurve newNurbsCurve;
	NuInt i, j;
	NuInt l1, l2;
	NuInt size, oldKnotCount, insert;

	mem = NuMem_new("NuNurbsCurve");
	size = sizeof(struct NuNurbsCurveStruct);
	newNurbsCurve = (NuNurbsCurve)NuMem_allocate(mem, size);
	newNurbsCurve->mem = mem;
	newNurbsCurve->degree = nurbsCurve->degree;
	newNurbsCurve->knotCount = nurbsCurve->knotCount + knotCount;
	size = newNurbsCurve->knotCount * sizeof(NuDouble);
	newNurbsCurve->knot = (NuDouble *)NuMem_allocate(mem, size);
	newNurbsCurve->pointsCount = nurbsCurve->pointsCount + knotCount;
	size = sizeof(NuDouble) * newNurbsCurve->pointsCount * 4;
	newNurbsCurve->points = (NuDouble (*)[4])NuMem_allocate(mem, size);

	i = newNurbsCurve->knotCount - nurbsCurve->knotCount - 1;
	oldKnotCount = nurbsCurve->knotCount - 1;
	for (j = newNurbsCurve->knotCount - 1; j >= 0; j--)
		{
		if ((j != newNurbsCurve->knotCount - 1) &&
			(j - newNurbsCurve->degree >= 0))
			{
			l1 = j - newNurbsCurve->degree;
			l2 = oldKnotCount - nurbsCurve->degree;
			newNurbsCurve->points[l1][0] = nurbsCurve->points[l2][0];
			newNurbsCurve->points[l1][1] = nurbsCurve->points[l2][1];
			newNurbsCurve->points[l1][2] = nurbsCurve->points[l2][2];
			newNurbsCurve->points[l1][3] = nurbsCurve->points[l2][3];
			}
		insert = 0;
		if (i >= 0)
			if (knot[i] >= nurbsCurve->knot[oldKnotCount])
				insert = 1;
		if (insert)
			{
			NuNurbsCurveInsertKnot(newNurbsCurve, nurbsCurve, knot,
				i, j, oldKnotCount);
			newNurbsCurve->knot[j] = knot[i--];
			}
		else
			newNurbsCurve->knot[j] = nurbsCurve->knot[oldKnotCount--];
		}
	return(newNurbsCurve);
	}

NuNurbsCurve NuPROC NuNurbsCurveNewPolyline(NuDouble points[][3],
	NuInt pointsCount)
	{
	return 0;
	}

NuNurbsCurve NuPROC NuNurbsCurveNewRefined(NuNurbsCurve nurbsCurve,
	NuInt factor)
	{
	NuNurbsCurve newNurbsCurve;
	NuInt knotCount, size;
	NuInt i, j;
	NuDouble step, *knot;

	size = (nurbsCurve->knotCount - 1) * factor * sizeof(NuDouble);
	knot = (NuDouble *)NuMem_allocate(nurbsCurve->mem, size);
	knotCount = 0;
	for (j = nurbsCurve->degree; j < nurbsCurve->pointsCount; j++)
		{
		step = (nurbsCurve->knot[j + 1] - nurbsCurve->knot[j]) /
			((NuDouble)factor + 1.0);
		if (step > 0.0)
			for (i = 0; i < factor; i++)
				knot[knotCount++] = nurbsCurve->knot[j] +
					step * (NuDouble)(i + 1);
		}
	newNurbsCurve = NuNurbsCurveNewInsertKnots(nurbsCurve, knot, knotCount);
	NuMem_deallocate(nurbsCurve->mem, (char *)knot);
	return(newNurbsCurve);
	}

NuInt NuPROC NuNurbsCurveSetKnot(NuNurbsCurve nurbsCurve, NuDouble knot[],
	NuInt knotCount)
	{
	NuInt i;

	if (knotCount != nurbsCurve->knotCount)
		return -1;
	for (i = 0; i < knotCount; i++)
		nurbsCurve->knot[i] = knot[i];
	}

NuInt NuPROC NuNurbsCurveSetPoints(NuNurbsCurve nurbsCurve,
        NuDouble points[][4], NuInt pointsCount)
	{
	NuInt i;

	if (pointsCount != nurbsCurve->pointsCount)
		return -1;
	for (i = 0; i < pointsCount; i++)
		{
		nurbsCurve->points[i][0] = points[i][0];
		nurbsCurve->points[i][1] = points[i][1];
		nurbsCurve->points[i][2] = points[i][2];
		nurbsCurve->points[i][3] = points[i][3];
		}
	}

static void NuNurbsEvaluate(NuDouble ePoints[][4], NuDouble knot[],
	NuDouble parameter, NuInt degree, NuInt knotPos, NuDouble point[3])
	{
	NuInt i, j, s, s1;
	NuDouble o, o1;

	for (j = degree + 1; j > 1; j--)
		{
		i = knotPos - 1;
		s1 = 1;
		for (s=0; s < j - 1; s++)
			{
			o = (knot[i + j - 1] - knot[i]);
			if (o > 0.0)
				{
				o = (parameter - knot[i]) / o;
				o1 = 1.0 - o;
				ePoints[s][0] = o * ePoints[s][0] + o1 * ePoints[s1][0];
				ePoints[s][1] = o * ePoints[s][1] + o1 * ePoints[s1][1];
				ePoints[s][2] = o * ePoints[s][2] + o1 * ePoints[s1][2];
				ePoints[s][3] = o * ePoints[s][3] + o1 * ePoints[s1][3];
				}
			else
				{
				ePoints[s][0] = ePoints[s1][0];
				ePoints[s][1] = ePoints[s1][1];
				ePoints[s][2] = ePoints[s1][2];
				ePoints[s][3] = ePoints[s1][3];
				}
			i--;
			s1++;
			}
		}
	if (ePoints[0][3] == 0.0)
		{
#ifdef DEBUG
		printf ("major error in eval %f %f\n",t,ePoints[0][3]);
#endif
		return;
		}
	if (ePoints[0][3] < 0.0)
		ePoints[0][3] = -ePoints[0][3];
	point[0] = ePoints[0][0] / ePoints[0][3];
	point[1] = ePoints[0][1] / ePoints[0][3];
	point[2] = ePoints[0][2] / ePoints[0][3];
	}

NuInt NuPROC NuNurbsSurfEvaluate(NuNurbsSurf nurbsSurf,
	NuDouble uParameter, NuDouble vParameter, NuDouble point[3])
/*	Algorithm from "Curve and surface constuctions using rational
	B-splines", Computer Aided Design vol 19 num 3 Nov 1987,
	pg 489.
*/	{
	NuInt i, j, k, l, r, s, s1, size, t;
	NuDouble o, o1;
	NuDouble (*ePoints)[4];

	k = nurbsSurf->uDegree;
	while (uParameter >= nurbsSurf->uKnot[k++])
		;
	k--;
	l = nurbsSurf->vDegree;
	while (vParameter >= nurbsSurf->vKnot[l++])
		;
	l--;
	size = (nurbsSurf->uDegree + 1) * (nurbsSurf->vDegree + 1) *
		sizeof(NuDouble) * 4;
	ePoints = (NuDouble (*)[4])NuMem_allocate(nurbsSurf->mem, size);
	i = 0;
	j = (k - 1) * nurbsSurf->vPointsCount;
	for (r = 0; r < nurbsSurf->uDegree + 1; r++)
		{
		t = j + l - 1;
		for (s = 0; s < nurbsSurf->vDegree + 1; s++)
			{
			ePoints[i][3] = nurbsSurf->points[t][3];
			ePoints[i][0] = nurbsSurf->points[t][0] * ePoints[i][3];
			ePoints[i][1] = nurbsSurf->points[t][1] * ePoints[i][3];
			ePoints[i][2] = nurbsSurf->points[t][2] * ePoints[i][3];
			i++;
			t--;
			}
		j -= nurbsSurf->vPointsCount;
		}
	for (r = nurbsSurf->uDegree + 1; r > 1; r--)
		{
		i = k - 1;
		j = 0;
		s1 = nurbsSurf->vDegree + 1;
		for (s = 0; s < r - 1; s++)
			{
			o = (nurbsSurf->uKnot[i + r - 1] - nurbsSurf->uKnot[i]);
			if (o > 0.0)
				{
				o = (uParameter - nurbsSurf->uKnot[i]) / o;
				o1 = 1.0 - o;
				for (t = 0; t < nurbsSurf->vDegree + 1; t++)
					{
					ePoints[j][0] = o * ePoints[j][0] + o1 * ePoints[s1][0];
					ePoints[j][1] = o * ePoints[j][1] + o1 * ePoints[s1][1];
					ePoints[j][2] = o * ePoints[j][2] + o1 * ePoints[s1][2];
					ePoints[j][3] = o * ePoints[j][3] + o1 * ePoints[s1][3];
					j++;
					s1++;
					}
				}
			else
				for (t = 0; t < nurbsSurf->vDegree + 1; t++)
					{
					ePoints[j][0] = ePoints[s1][0];
					ePoints[j][1] = ePoints[s1][1];
					ePoints[j][2] = ePoints[s1][2];
					ePoints[j][3] = ePoints[s1][3];
					j++;
					s1++;
					}
			i--;
			}
		}
	NuNurbsEvaluate(ePoints, nurbsSurf->vKnot, vParameter,
		nurbsSurf->vDegree, l, point);
	NuMem_deallocate(nurbsSurf->mem, (char *)ePoints);
	return 0;
	}

void NuPROC NuNurbsSurfFree(NuNurbsSurf nurbsSurf)
	{
	NuMem_free(nurbsSurf->mem);
	}

void NuPROC NuNurbsSurfGetDegrees(NuNurbsSurf nurbsSurf, NuInt *uDegree,
	NuInt *vDegree)
	{
	*uDegree = nurbsSurf->uDegree;
	*vDegree = nurbsSurf->vDegree;
	}

NuInt NuPROC NuNurbsSurfGetKnots(NuNurbsSurf nurbsSurf, NuDouble uKnot[],
	NuDouble vKnot[], NuInt uKnotCount, NuInt vKnotCount)
	{
	NuInt i;

	if (uKnotCount != nurbsSurf->uKnotCount)
		return -1;
	if (vKnotCount != nurbsSurf->vKnotCount)
		return -1;
	for (i = 0; i < uKnotCount; i++)
		uKnot[i] = nurbsSurf->uKnot[i];
	for (i = 0; i < vKnotCount; i++)
		vKnot[i] = nurbsSurf->vKnot[i];
	}

void NuPROC NuNurbsSurfGetKnotCounts(NuNurbsSurf nurbsSurf,
	NuInt *uKnotCount, NuInt *vKnotCount)
	{
	*uKnotCount = nurbsSurf->uKnotCount;
	*vKnotCount = nurbsSurf->vKnotCount;
	}

NuInt NuPROC NuNurbsSurfGetNormal(NuNurbsSurf nurbsSurf, NuDouble uParameter,
	NuDouble vParameter, NuDouble point[3])
	{
	return -1;
#ifdef NEVER
	generate u deriv
	generate v deriv - use eval ?
	cross product them
#endif
	}

NuInt NuPROC NuNurbSurfGetNormals(NuNurbsSurf nurbsSurf, NuDouble normals[][3],
	NuInt uPointsCount, NuInt vPointsCount)
	{
	return -1;
	}

NuInt NuPROC NuNurbsSurfGetPoints(NuNurbsSurf nurbsSurf,
	NuDouble points[][4], NuInt uPointsCount, NuInt vPointsCount)
	{
	NuInt i;

	if (uPointsCount != nurbsSurf->uPointsCount)
		return -1;
	if (vPointsCount != nurbsSurf->vPointsCount)
		return -1;
	for (i = 0; i < uPointsCount * vPointsCount; i++)
		{
		points[i][0] = nurbsSurf->points[i][0];
		points[i][1] = nurbsSurf->points[i][1];
		points[i][2] = nurbsSurf->points[i][2];
		points[i][3] = nurbsSurf->points[i][3];
		}
	return 0;
	}

void NuPROC NuNurbsSurfGetPointsCounts(NuNurbsSurf nurbsSurf,
	NuInt *uPointsCount, NuInt *vPointsCount)
	{
	*uPointsCount = nurbsSurf->uPointsCount;
	*vPointsCount = nurbsSurf->vPointsCount;
	}

static void NuNurbsSurfInsertUKnot(NuNurbsSurf newNurbsSurf,
	NuNurbsSurf nurbsSurf, NuDouble knot[], NuInt i, NuInt j,
	NuInt oldUKnotCount)
	{
	NuInt k, l, l1, l2;
	NuDouble B, Bw1, Bw2, Bwsum;

	for (l = 0; l < nurbsSurf->uDegree; l++)
		{
		B = newNurbsSurf->uKnot[j + l + 1] - knot[i];
		l1 = (j - nurbsSurf->uDegree + l) * newNurbsSurf->vPointsCount;
		l2 = l1 + newNurbsSurf->vPointsCount;
		if (B <= 0.0)
			{
			for (k = 0; k < newNurbsSurf->vPointsCount; k++)
				{
				newNurbsSurf->points[l1][0] = newNurbsSurf->points[l2][0];
				newNurbsSurf->points[l1][1] = newNurbsSurf->points[l2][1];
				newNurbsSurf->points[l1][2] = newNurbsSurf->points[l2][2];
				newNurbsSurf->points[l1][3] = newNurbsSurf->points[l2][3];
				l1++;
				l2++;
				}
			}
		else
			{
			B /= (newNurbsSurf->uKnot[j + l + 1] -
				nurbsSurf->uKnot[oldUKnotCount - nurbsSurf->uDegree + l + 1]);
			for (k = 0; k < newNurbsSurf->vPointsCount; k++)
				{
				Bw1 = B * newNurbsSurf->points[l1][3];
				Bw2 = (1.0 - B) * newNurbsSurf->points[l2][3];
				Bwsum = Bw1 + Bw2;
				newNurbsSurf->points[l1][0] =
					(Bw1 * newNurbsSurf->points[l1][0] +
					 Bw2 * newNurbsSurf->points[l2][0]) / Bwsum;
				newNurbsSurf->points[l1][1] =
					(Bw1 * newNurbsSurf->points[l1][1] +
					 Bw2 * newNurbsSurf->points[l2][1]) / Bwsum;
				newNurbsSurf->points[l1][2] =
					(Bw1 * newNurbsSurf->points[l1][2]+
					Bw2 * newNurbsSurf->points[l2][2]) / Bwsum;
				newNurbsSurf->points[l1][3] = Bwsum;
				l1++;
				l2++;
				}
			}
		}
	}

static void NuNurbsSurfInsertVKnot(NuNurbsSurf newNurbsSurf,
	NuNurbsSurf nurbsSurf, NuDouble knot[], NuInt i, NuInt j,
	NuInt oldVKnotCount)
	{
	NuInt k, l, l1, l2;
	NuDouble B, Bw1, Bw2, Bwsum;

	for (l = 0; l < nurbsSurf->vDegree; l++)
		{
		B = newNurbsSurf->vKnot[j + l + 1] - knot[i];
		l1 = j - nurbsSurf->vDegree + l;
		l2 = l1 + 1;
		if (B <= 0.0)
			{
			for (k = 0; k < newNurbsSurf->uPointsCount; k++)
				{
				newNurbsSurf->points[l1][0] = newNurbsSurf->points[l2][0];
				newNurbsSurf->points[l1][1] = newNurbsSurf->points[l2][1];
				newNurbsSurf->points[l1][2] = newNurbsSurf->points[l2][2];
				newNurbsSurf->points[l1][3] = newNurbsSurf->points[l2][3];
				l1 += newNurbsSurf->vPointsCount;
				l2 += newNurbsSurf->vPointsCount;
				}
			}
		else
			{
			B /= (newNurbsSurf->vKnot[j + l + 1] -
				nurbsSurf->vKnot[oldVKnotCount - nurbsSurf->vDegree + l + 1]);
			for (k = 0; k < newNurbsSurf->uPointsCount; k++)
				{
				Bw1 = B * newNurbsSurf->points[l1][3];
				Bw2 = (1.0 - B) * newNurbsSurf->points[l2][3];
				Bwsum = Bw1 + Bw2;
				newNurbsSurf->points[l1][0] =
					(Bw1 * newNurbsSurf->points[l1][0] +
					Bw2 * newNurbsSurf->points[l2][0]) / Bwsum;
				newNurbsSurf->points[l1][1] =
					(Bw1 * newNurbsSurf->points[l1][1] +
					Bw2 * newNurbsSurf->points[l2][1]) / Bwsum;
				newNurbsSurf->points[l1][2] =
					(Bw1 * newNurbsSurf->points[l1][2] +
					Bw2 * newNurbsSurf->points[l2][2]) / Bwsum;
				newNurbsSurf->points[l1][3] = Bwsum;
				l1 += newNurbsSurf->vPointsCount;
				l2 += newNurbsSurf->vPointsCount;
				}
			}
		}
	}

NuNurbsSurf NuPROC NuNurbsSurfNew(NuInt uDegree, NuInt vDegree,
	NuDouble uKnot[], NuDouble vKnot[], NuInt uKnotCount, NuInt vKnotCount,
	NuDouble points[][4], NuInt uPointsCount, NuInt vPointsCount)
	{
	NuMem mem;
	NuNurbsSurf nurbsSurf;
	NuInt size, i;

	if (uKnotCount != uDegree + uPointsCount + 1)
		return 0;
	if (vKnotCount != vDegree + vPointsCount + 1)
		return 0;
	mem = NuMem_new("NuNurbsSurf");
	size = sizeof(struct NuNurbsSurfStruct);
	nurbsSurf = (NuNurbsSurf)NuMem_allocate(mem, size);
	nurbsSurf->mem = mem;
	nurbsSurf->uDegree = uDegree;
	nurbsSurf->vDegree = vDegree;
	nurbsSurf->uKnotCount = uKnotCount;
	size = uKnotCount * sizeof(NuDouble);
	nurbsSurf->uKnot = (NuDouble *)NuMem_allocate(mem, size);
	for (i = 0; i < uKnotCount; i++)
		nurbsSurf->uKnot[i] = uKnot[i];
	nurbsSurf->vKnotCount = vKnotCount;
	size = vKnotCount * sizeof(NuDouble);
	nurbsSurf->vKnot = (NuDouble *)NuMem_allocate(mem, size);
	for (i = 0; i < vKnotCount; i++)
		nurbsSurf->vKnot[i] = vKnot[i];
	nurbsSurf->uPointsCount = uPointsCount;
	nurbsSurf->vPointsCount = vPointsCount;
	size = uPointsCount * vPointsCount * sizeof(NuDouble) * 4;
	nurbsSurf->points = (NuDouble (*)[4])NuMem_allocate(mem, size);
	for (i = 0; i < uPointsCount * vPointsCount; i++)
		{
		nurbsSurf->points[i][0] = points[i][0];
		nurbsSurf->points[i][1] = points[i][1];
		nurbsSurf->points[i][2] = points[i][2];
		nurbsSurf->points[i][3] = points[i][3];
		}
	return(nurbsSurf);
	}

NuNurbsSurf NuPROC NuNurbsSurfNewCylinder(NuDouble center[3],
	NuDouble normal[3], NuDouble radius, NuDouble distance)
	{
	return 0;
	}

NuNurbsSurf NuPROC NuNurbsSurfNewCone(NuDouble center[3],
	NuDouble normal[3], NuDouble radius1, NuDouble radius2,
	NuDouble distance)
	{
	return 0;
	}

NuNurbsSurf NuPROC NuNurbsSurfNewExtruded(NuNurbsCurve nurbsCurve,
	NuDouble vec[3], NuDouble distance)
	{
	return 0;

#ifdef NEVER
	NuNurbsSurf ns;
	static NuDouble uknot[4]= {0.0, 0.0, 1.0, 1.0};
	NuDouble (*points)[4];
	NuInt loop, size;

	size = sizeof(NuDouble) * nc->vPointsCount * 4 * 2;
	points = (NuDouble (*)[4])NuMem_allocate(nc->mem, size);
	for (loop = 0; loop < nc->vPointsCount; loop++)
		{
		points[loop*2][0] = nc->points[loop][0];
		points[loop*2][1] = nc->points[loop][1];
		points[loop*2][2] = nc->points[loop][2] - .5;
		points[loop*2][3] = nc->points[loop][3];
		points[loop*2+1][0] = nc->points[loop][0];
		points[loop*2+1][1] = nc->points[loop][1];
		points[loop*2+1][2] = nc->points[loop][2] + .5;
		points[loop*2+1][3] = nc->points[loop][3];
		}
	ns = NuNurbsSurf_new(nc->p, (NuInt)1, nc->knot, uknot, points, nc->vPointsCount, (NuInt)2);
	NuMem_deallocate(nc->mem, (char *)points);
	return(ns);
#endif
	}

NuNurbsSurf NuPROC NuNurbsSurfNewInsertUKnots(NuNurbsSurf nurbsSurf,
	NuDouble knot[], NuInt knotCount)
/*	Boehms algorithm.  See "Inserting new knots into B-Spline Curves",
	Computer Aided Design Vol 12 No 4 (1980), pg 199-202.  This
	algorithm is the modified Boehm algorithm from "The Insertion
	Algorithm", Computer Aided Design Vol 17 (1985), pg 58-59.
*/	{
	NuMem mem;
	NuNurbsSurf newNurbsSurf;
	NuInt i, j, k, l1, l2;
	NuInt size, oldUKnotCount, insert;

	mem = NuMem_new("NuNurbsSurf");
	size = sizeof(struct NuNurbsSurfStruct);
	newNurbsSurf = (NuNurbsSurf)NuMem_allocate(mem, size);
	newNurbsSurf->mem = mem;
	newNurbsSurf->uDegree = nurbsSurf->uDegree;
	newNurbsSurf->vDegree = nurbsSurf->vDegree;
	newNurbsSurf->uKnotCount = knotCount + nurbsSurf->uKnotCount;
	size = newNurbsSurf->uKnotCount * sizeof(NuDouble);
	newNurbsSurf->uKnot = (NuDouble *)NuMem_allocate(mem, size);
	newNurbsSurf->vKnotCount = nurbsSurf->vKnotCount;
	size = newNurbsSurf->vKnotCount * sizeof(NuDouble);
	newNurbsSurf->vKnot = (NuDouble *)NuMem_allocate(mem, size);
	for (j=0; j < nurbsSurf->vKnotCount; j++)
		newNurbsSurf->vKnot[j] = nurbsSurf->vKnot[j];
	newNurbsSurf->uPointsCount = newNurbsSurf->uKnotCount -
		nurbsSurf->uDegree - 1;
	newNurbsSurf->vPointsCount = nurbsSurf->vPointsCount;
	size = sizeof(NuDouble) * newNurbsSurf->vPointsCount *
		newNurbsSurf->uPointsCount * 4;
	newNurbsSurf->points = (NuDouble (*)[4])NuMem_allocate(mem, size);

	i = newNurbsSurf->uKnotCount - nurbsSurf->uKnotCount - 1;
	oldUKnotCount = nurbsSurf->uKnotCount - 1;
	for (j = newNurbsSurf->uKnotCount - 1; j >= 0; j--)
		{
		if ((j != newNurbsSurf->uKnotCount - 1) &&
			(j - newNurbsSurf->uDegree >= 0))
			{
			l1 = (j - newNurbsSurf->uDegree) * newNurbsSurf->vPointsCount;
			l2 = (oldUKnotCount - nurbsSurf->uDegree) * nurbsSurf->vPointsCount;
			for (k = 0; k < newNurbsSurf->vPointsCount; k++)
				{
				newNurbsSurf->points[l1][0] = nurbsSurf->points[l2][0];
				newNurbsSurf->points[l1][1] = nurbsSurf->points[l2][1];
				newNurbsSurf->points[l1][2] = nurbsSurf->points[l2][2];
				newNurbsSurf->points[l1][3] = nurbsSurf->points[l2][3];
				l1++;
				l2++;
				}
			}
		insert = 0;
		if (i >= 0)
			if (knot[i] >= nurbsSurf->uKnot[oldUKnotCount])
				insert = 1;
		if (insert)
			{
			NuNurbsSurfInsertUKnot(newNurbsSurf, nurbsSurf, knot,
				i, j, oldUKnotCount);
			newNurbsSurf->uKnot[j] = knot[i--];
			}
		else
			newNurbsSurf->uKnot[j] = nurbsSurf->uKnot[oldUKnotCount--];
		}
	return(newNurbsSurf);
	}

NuNurbsSurf NuPROC NuNurbsSurfNewInsertVKnots(NuNurbsSurf nurbsSurf,
	NuDouble knot[], NuInt knotCount)
/*	Boehms algorithm.  See "Inserting new knots into B-Spline Curves",
	Computer Aided Design Vol 12 No 4 (1980), pg 199-202.  This
	algorithm is the modified Boehm algorithm from "The Insertion
	Algorithm", Computer Aided Design Vol 17 (1985), pg 58-59.
*/	{
	NuMem mem;
	NuNurbsSurf newNurbsSurf;
	NuInt i, j, k, l1, l2;
	NuInt size, oldVKnotCount, insert;

	mem = NuMem_new("NuNurbsSurf");
	size = sizeof(struct NuNurbsSurfStruct);
	newNurbsSurf = (NuNurbsSurf)NuMem_allocate(mem, size);
	newNurbsSurf->mem = mem;
	newNurbsSurf->uDegree = nurbsSurf->uDegree;
	newNurbsSurf->vDegree = nurbsSurf->vDegree;
	newNurbsSurf->uKnotCount = nurbsSurf->uKnotCount;
	size = newNurbsSurf->uKnotCount * sizeof(NuDouble);
	newNurbsSurf->uKnot = (NuDouble *)NuMem_allocate(mem, size);
	for (j=0; j < nurbsSurf->uKnotCount; j++)
		newNurbsSurf->uKnot[j] = nurbsSurf->uKnot[j];
	newNurbsSurf->vKnotCount = knotCount + nurbsSurf->vKnotCount;
	size = newNurbsSurf->vKnotCount * sizeof(NuDouble);
	newNurbsSurf->vKnot = (NuDouble *)NuMem_allocate(mem, size);
	newNurbsSurf->uPointsCount = nurbsSurf->uPointsCount;
	newNurbsSurf->vPointsCount = newNurbsSurf->vKnotCount - nurbsSurf->vDegree - 1;
	size = sizeof(NuDouble) * newNurbsSurf->vPointsCount * newNurbsSurf->uPointsCount * 4;
	newNurbsSurf->points = (NuDouble (*)[4])NuMem_allocate(mem, size);

	i = newNurbsSurf->vKnotCount - nurbsSurf->vKnotCount -1;
	oldVKnotCount = nurbsSurf->vKnotCount - 1;
	for (j = newNurbsSurf->vKnotCount - 1; j >= 0; j--)
		{
		if ((j != newNurbsSurf->vKnotCount - 1) &&
			(j - newNurbsSurf->vDegree >= 0))
			{
			l1 = j - newNurbsSurf->vDegree;
			l2 = oldVKnotCount - nurbsSurf->vDegree;
			for (k = 0; k < newNurbsSurf->uPointsCount; k++)
				{
				newNurbsSurf->points[l1][0] = nurbsSurf->points[l2][0];
				newNurbsSurf->points[l1][1] = nurbsSurf->points[l2][1];
				newNurbsSurf->points[l1][2] = nurbsSurf->points[l2][2];
				newNurbsSurf->points[l1][3] = nurbsSurf->points[l2][3];
				l1 += newNurbsSurf->vPointsCount;
				l2 += nurbsSurf->vPointsCount;
				}
			}
		insert = 0;
		if (i >= 0)
			if (knot[i] >= nurbsSurf->vKnot[oldVKnotCount])
				insert = 1;
		if (insert)
			{
			NuNurbsSurfInsertVKnot(newNurbsSurf, nurbsSurf, knot,
				i, j, oldVKnotCount);
			newNurbsSurf->vKnot[j] = knot[i--];
			}
		else
			newNurbsSurf->vKnot[j] = nurbsSurf->vKnot[oldVKnotCount--];
		}
	return(newNurbsSurf);
	}

NuNurbsSurf NuPROC NuNurbsSurfNewMesh(NuDouble points[][3],
	NuInt uPointsCount, NuInt vPointsCount)
	{
	return 0;
	}

NuNurbsSurf NuPROC NuNurbsSurfNewPlane(NuDouble point1[3],
	NuDouble point2[3], NuDouble point3[3])
	{
	return 0;
	}

NuNurbsSurf NuPROC NuNurbsSurfNewRefined(NuNurbsSurf nurbsSurf, NuInt factor)
	{
	NuNurbsSurf newNurbsSurf1, newNurbsSurf2;
	NuDouble *knot, step;
	NuInt count, size;
	NuInt i, j;

	/* generate insert u knot vector */
	size = (nurbsSurf->uKnotCount - 1) * factor * sizeof(NuDouble);
	knot = (NuDouble *)NuMem_allocate(nurbsSurf->mem, size);
	count = 0;
	for (j = nurbsSurf->uDegree; j < nurbsSurf->uPointsCount; j++)
		{
		step = (nurbsSurf->uKnot[j + 1] - nurbsSurf->uKnot[j])
			/ ((NuDouble)factor + 1.0);
		if (step > 0.0)
			for (i = 0; i < factor; i++)
				knot[count++] = nurbsSurf->uKnot[j] + step * (NuDouble)(i + 1);
		}
	/* insert u knots */
	newNurbsSurf1 = NuNurbsSurfNewInsertUKnots(nurbsSurf, knot, count);
	NuMem_deallocate(nurbsSurf->mem, (char *)knot);
	/* generate insert v knot vector */
	size = (nurbsSurf->vKnotCount - 1) * factor * sizeof(NuDouble);
	knot = (NuDouble *)NuMem_allocate(nurbsSurf->mem, size);
	count = 0;
	for (j = nurbsSurf->vDegree; j < nurbsSurf->vPointsCount; j++)
		{
		step = (nurbsSurf->vKnot[j + 1] - nurbsSurf->vKnot[j])
			/ ((NuDouble)factor + 1.0);
		if (step > 0.0)
			for (i = 0; i < factor; i++)
				knot[count++] = nurbsSurf->vKnot[j] + step * (NuDouble)(i + 1);
		}
	/* insert v knots */
	newNurbsSurf2 = NuNurbsSurfNewInsertVKnots(newNurbsSurf1, knot, count);
	NuMem_deallocate(nurbsSurf->mem, (char *)knot);
	NuNurbsSurfFree(newNurbsSurf1);
	return(newNurbsSurf2);
	}

NuNurbsSurf NuPROC NuNurbsSurfNewRev(NuNurbsCurve nurbsCurve,
	NuDouble point1[3], NuDouble point2[3], NuDouble degrees)
	{
	NuNurbsSurf nurbsSurf;
	NuDouble (*points)[4];
	NuInt loop, loopc, size, i;
	static NuDouble vKnot[10]=
		{0.0, 0.0, 0.0, .25, .5, .5, .75, 1.0, 1.0, 1.0};
	static NuDouble scalePoints[7][4] =
		{
		{ 1.0,   0.0,   0.0,   1.0},
		{ 1.0,   0.0,   1.0,    .5},
		{-1.0,   0.0,   1.0,    .5},
		{-1.0,   0.0,   0.0,   1.0},
		{-1.0,   0.0,  -1.0,    .5},
		{ 1.0,   0.0,  -1.0,    .5},
		{ 1.0,   0.0,   0.0,   1.0}
		};

	size = sizeof(NuDouble) * nurbsCurve->pointsCount * 4 * 7;
	points = (NuDouble (*)[4])NuMem_allocate(nurbsCurve->mem, size);
	i = 0;
	for (loop = 0; loop < nurbsCurve->pointsCount; loop++)
		for (loopc = 0; loopc < 7; loopc++)
			{
			points[i][0] = scalePoints[loopc][0] * nurbsCurve->points[loop][0];
			points[i][1] = nurbsCurve->points[loop][1];
			points[i][2] = scalePoints[loopc][2] * nurbsCurve->points[loop][0];
			points[i][3] = scalePoints[loopc][3] * nurbsCurve->points[loop][3];
			i++;
			}
	nurbsSurf = NuNurbsSurfNew(nurbsCurve->degree, 2, nurbsCurve->knot, vKnot,
		nurbsCurve->knotCount, 10, points, nurbsCurve->pointsCount, 7);
	NuMem_deallocate(nurbsCurve->mem, (char *)points);
	return(nurbsSurf);
	}

NuNurbsSurf NuPROC NuNurbsSurfNewRuled(NuNurbsCurve nurbsCurve1,
	NuNurbsCurve nurbsCurve2)
	{
	return 0;
	}

NuNurbsSurf NuPROC NuNurbsSurfNewSphere(NuDouble center[3], NuDouble radius)
	{
	NuInt loop;
	NuNurbsCurve nurbsCurve;
	NuNurbsSurf nurbsSurf;
	static NuDouble knot[7]=
		{0.0, 0.0, 0.0, .5, 1.0, 1.0, 1.0};
	static NuDouble points[4][4] =
		{
		{ 0.0, -1.0,  0.0,  1.0},
		{ 1.0, -1.0,  0.0,   .5},
		{ 1.0,  1.0,  0.0,   .5},
		{ 0.0,  1.0,  0.0,  1.0}
		};
	static NuDouble point1[3], point2[3];

	for (loop=0; loop<4; loop++)
		{
		points[loop][0] *= radius;
		points[loop][1] *= radius;
		}
	nurbsCurve = NuNurbsCurveNew(2, knot, 7, points, 4);
	nurbsSurf = NuNurbsSurfNewRev(nurbsCurve, point1, point2, 360.0);
	NuNurbsCurveFree(nurbsCurve);
	return(nurbsSurf);
	}

NuInt NuPROC NuNurbsSurfSetPoints(NuNurbsSurf nurbsSurf,
	NuDouble points[][4], NuInt uPointsCount, NuInt vPointsCount)
	{
	NuInt i;

	if (uPointsCount != nurbsSurf->uPointsCount)
		return -1;
	if (vPointsCount != nurbsSurf->vPointsCount)
		return -1;
	for (i = 0; i < uPointsCount * vPointsCount; i++)
		{
		nurbsSurf->points[i][0] = points[i][0];
		nurbsSurf->points[i][1] = points[i][1];
		nurbsSurf->points[i][2] = points[i][2];
		nurbsSurf->points[i][3] = points[i][3];
		}
	return 0;
	}

NuInt NuPROC NuNurbsSurfSetKnots(NuNurbsSurf nurbsSurf, NuDouble uKnot[],
        NuDouble vKnot[], NuInt uKnotCount, NuInt vKnotCount)
	{
	NuInt i;

	if (uKnotCount != nurbsSurf->uKnotCount)
		return -1;
	if (vKnotCount != nurbsSurf->vKnotCount)
		return -1;
	for (i = 0; i < uKnotCount; i++)
		nurbsSurf->uKnot[i] = uKnot[i];
	for (i = 0; i < vKnotCount; i++)
		nurbsSurf->vKnot[i] = vKnot[i];
	}
