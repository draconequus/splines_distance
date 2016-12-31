#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include "cubic_splines.h"
#include "relative_position.h"

int ReadPointsFromFile(std::vector<point2d> *ReferencePoints, std::string path);

int main()
{
    std::vector<point2d> RefPoints1, RefPoints2;
    std::string path1("ReferencePointsList1.txt");
    std::string path2("ReferencePointsList2.txt");
    ReadPointsFromFile(&RefPoints1, path1);
    if (RefPoints1.size() < 3) {
        std::cout << "Only " << RefPoints1.size() << " point(s) in file " << path1 << " while 3 or more required.\n";
        return 1;
        }
    ReadPointsFromFile(&RefPoints2, path2);
    if (RefPoints2.size() < 3) {
        std::cout << "Only " << RefPoints2.size() << " point(s) in file " << path2 << " while 3 or more required.\n";
        return 1;
        }


    CubicSpline SplineA(&RefPoints1);
    CubicSpline SplineB(&RefPoints2);

    std::ofstream OutputFile1, OutputFile2;
    OutputFile1.open("spline1.txt");
    OutputFile2.open("spline2.txt");
    double ShiftedT = 0;
    unsigned int IMAX = 200;
    for(unsigned int i=0; i<=IMAX; i++) {
        ShiftedT = SplineA.Tmin() + (SplineA.Tmax()-SplineA.Tmin())*(double)i/(double)IMAX;
        OutputFile1 /*<< ShiftedT << ' '*/ << SplineA.X(ShiftedT) << ' ' << SplineA.Y(ShiftedT) << '\n';
        OutputFile2 /*<< ShiftedT << ' '*/ << SplineB.X(ShiftedT) << ' ' << SplineB.Y(ShiftedT) << '\n';
    }
    OutputFile1.close();
    OutputFile2.close();
    RelativePositionDescriptor desc1(&SplineA, &SplineA);
    std::cout << "Intersection points found: " << desc1.GetIntersectionsList()->size() << '\n';
    desc1.InitializeIntersectionsList();
    std::cout << "Intersection points found: " << desc1.GetIntersectionsList()->size() << '\n';
    std::cout << "X: " << desc1.GetIntersectionsList()->at(0).PointA.x << " Y: " << desc1.GetIntersectionsList()->at(0).PointA.y << '\n';
    std::cout << "X: " << desc1.GetIntersectionsList()->at(0).PointB.x << " Y: " << desc1.GetIntersectionsList()->at(0).PointB.y << '\n';



//    parametricvalue temp_p;
//    temp_p = SplineA.MinXInRange(0,1);
//    std::cout << "Minimum X = " << temp_p.value << " at t = " << temp_p.t << '\n';
//    temp_p = SplineA.MinYInRange(0, 1);
//    std::cout << "Minimum Y = " << temp_p.value << " at t = " << temp_p.t << '\n';
//    temp_p = SplineA.MaxXInRange(0, 1);
//    std::cout << "Maximum X = " << temp_p.value << " at t = " << temp_p.t << '\n';
//    temp_p = SplineA.MaxYInRange(0,1);
//    std::cout << "Maximum Y = " << temp_p.value << " at t = " << temp_p.t << '\n';
//    //std::cout << "X = " << SplineA.X(0.148) << " at t = " << 0.148 << '\n';
//    //std::cout << "X = " << SplineA.X(0.149) << " at t = " << 0.149 << '\n';

    return 0;
}


int ReadPointsFromFile(std::vector<point2d> *ReferencePoints, std::string path) {
    point2d currentPoint;
    std::ifstream infile;
    infile.open(path.c_str());
    while (infile >> currentPoint.x >> currentPoint.y)
        ReferencePoints->push_back(currentPoint);
    infile.close();
    return 0;
}
