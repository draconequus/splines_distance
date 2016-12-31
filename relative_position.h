#ifndef RELATIVE_POSITION_H_INCLUDED
#define RELATIVE_POSITION_H_INCLUDED
#define MESH_MAX 100

#define INTERSECTIONS_EPSILON (double)1e-09

#include<cmath>
#include "cubic_splines.h"


typedef struct ParamPointCouple {
    parametricpoint PointA;
    parametricpoint PointB;
} parametricpointcouple;

class RelativePositionDescriptor {
public:
    RelativePositionDescriptor(CubicSpline* Spline1, CubicSpline* Spline2);
    void InitializeIntersectionsList();
    std::vector<parametricpointcouple>* GetIntersectionsList();



private:
    CubicSpline *SplineA, *SplineB;
    bool IntersectionsFound;
    bool MinimalDistanceFound;
    bool DoRectanglesIntersect(double t1min, double t1max, double t2min, double t2max);
    std::vector<parametricpointcouple> IntersectionsList;
    void FindIntersections(double t1min, double t1max, double t2min, double t2max);



    };

std::vector<parametricpointcouple>* RelativePositionDescriptor::GetIntersectionsList() {
    return &IntersectionsList;
}

void RelativePositionDescriptor::InitializeIntersectionsList() {
    IntersectionsList.clear();
    FindIntersections(SplineA->Tmin(), SplineA->Tmax(), SplineB->Tmin(), SplineB->Tmax());
}

void RelativePositionDescriptor::FindIntersections(double t1min, double t1max, double t2min, double t2max) {
    if (!DoRectanglesIntersect(t1min, t1max, t2min, t2max)) return;

    double Precision1X = abs(SplineA->X(t1min) - SplineA->X(t1max));
    double Precision1Y = abs(SplineA->Y(t1min) - SplineA->Y(t1max));
    double Precision2X = abs(SplineB->X(t2min) - SplineB->X(t2max));
    double Precision2Y = abs(SplineB->Y(t2min) - SplineB->Y(t2max));

    if ((Precision1X < INTERSECTIONS_EPSILON) && (Precision1Y < INTERSECTIONS_EPSILON) && (Precision2X < INTERSECTIONS_EPSILON) && (Precision2Y < INTERSECTIONS_EPSILON)) {
        parametricpointcouple NewIntersection;
        NewIntersection.PointA.t = (t1max - t1min)/2;
        NewIntersection.PointA.x = SplineA->X(NewIntersection.PointA.t);
        NewIntersection.PointA.y = SplineA->Y(NewIntersection.PointA.t);
        NewIntersection.PointB.t = (t1max - t1min)/2;
        NewIntersection.PointB.x = SplineB->X(NewIntersection.PointB.t);
        NewIntersection.PointB.y = SplineB->Y(NewIntersection.PointB.t);
        IntersectionsList.push_back(NewIntersection);
    } else {
        if (Precision1X + Precision1Y > Precision2X + Precision2Y) {
            FindIntersections(t1min, (t1min+t1max)/2, t2min, t2max);
            FindIntersections((t1min+t1max)/2, t1max, t2min, t2max);
        } else {
            FindIntersections(t1min, t1max, t2min, (t2min+t2max)/2);
            FindIntersections(t1min, t1max, (t2min+t2max)/2, t2max);
        }

    }
    IntersectionsFound = true;
}


RelativePositionDescriptor::RelativePositionDescriptor(CubicSpline* Spline1, CubicSpline* Spline2) {
    SplineA = Spline1;
    SplineB = Spline2;
    IntersectionsFound = false;
    MinimalDistanceFound = false;

}


bool RelativePositionDescriptor::DoRectanglesIntersect(double t1min, double t1max, double t2min, double t2max) {
    double LeftBorder1 = SplineA->MinXInRange(t1min, t1max).value;
    double RightBorder1 = SplineA->MaxXInRange(t1min, t1max).value;
    double BottomBorder1 = SplineA->MinYInRange(t1min, t1max).value;
    double TopBorder1 = SplineA->MaxYInRange(t1min, t1max).value;
    double LeftBorder2 = SplineB->MinXInRange(t2min, t2max).value;
    double RightBorder2 = SplineB->MaxXInRange(t2min, t2max).value;
    double BottomBorder2 = SplineB->MinYInRange(t2min, t2max).value;
    double TopBorder2 = SplineB->MaxYInRange(t2min, t2max).value;

    bool Result;
    if ((RightBorder2 < LeftBorder1) || (RightBorder1 < LeftBorder2) || (TopBorder2 < BottomBorder1) || (TopBorder1 < BottomBorder2))
        Result = false;
    else
        Result = true;


    return Result;
}



#endif // RELATIVE_POSITION_H_INCLUDED
