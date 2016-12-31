#ifndef CUBIC_SPLINES_H
#define CUBIC_SPLINES_H
#include <cassert>
#include <vector>
#include <math.h>

// ���������, ������������ � ����

#define DEFAULT_DENSITY 100
#define DEFAULT_TMIN (double)0
#define DEFAULT_TMAX (double)1
#define DEFAULT_BOUNDARY_CONDITIONS 0

#define IN_X 0
#define IN_Y 1
#define MIN_VALUE 0
#define MAX_VALUE 1



// ����������� �����


typedef struct PointIn2D {
    double x;
    double y;
} point2d;

typedef struct ParamValue {
    double value;
    double t;
    } parametricvalue;

typedef struct ParamPoint {
    double x;
    double y;
    double t;
} parametricpoint;

typedef struct CoefficientsOfaCubicPolynome {
    double A;
    double B;
    double C;
    double D;
} coefficients;

typedef struct tridiagmatr {
    std::vector<double> *LowerDiagonal;
    std::vector<double> *MainDiagonal;
    std::vector<double> *UpperDiagonal;
} TridiagonalMatrix;

// ����� ����������� �����


class CubicSpline
{
public:
    CubicSpline(std::vector<point2d>* ReferencePoints, double Tmin = DEFAULT_TMIN, double Tmax = DEFAULT_TMAX);
    double Tmin();                          // ����������� �������� ���������
    double Tmax();                          // ������������ �������� ���������
    std::vector<point2d> value(double t);   // ������� ����� ������, ��������������� �������� ���������
    double X(double t);
    double Y(double t);
    int test();
    int RangeNumber(double T);                                                             // ������� ������ ���������, � �������� ��������� ������ �������� t
    void SetCoefficients(int BoundaryConditions = DEFAULT_BOUNDARY_CONDITIONS);        // ������� ���������������� ������� ������������ 2-� ����������� � ������.
                                                                                            // �������� � ��������� ������� ��� �������� �����������.
    parametricvalue ExtremalValueInRange(int RangeType, int ExtremaType, double Tleft, double Tright);
    parametricvalue MinXInRange(double Tleft, double Tright);
    parametricvalue MaxXInRange(double Tleft, double Tright);
    parametricvalue MinYInRange(double Tleft, double Tright);
    parametricvalue MaxYInRange(double Tleft, double Tright);




private:
    double T_min;
    double T_max;
    int NumberOfPoints;
    std::vector<point2d>* RefPoints;          // ������, ��������� �� ������� ����� �������
    std::vector<double> *Ti, *XSecondDerivative, *YSecondDerivative;
    std::vector<coefficients> *XCoeff, *YCoeff;
    void SolveTridiagonalSystem(TridiagonalMatrix* matrix, std::vector<double>* var, std::vector<double>* RightSide);
    int renormalize(int density = DEFAULT_DENSITY);
    std::vector<double> PossibleExtremaInX(double Tleft, double Tright);
    std::vector<double> PossibleExtremaInY(double Tleft, double Tright);
    std::vector<double> PossibleExtremaInRange(double Tleft, double Tright, int type);
    void InitializeExtremaList(int type);


    std::vector<double> PossExtremaInX;
    std::vector<double> PossExtremaInY;
    bool PossExtremaInXSet;
    bool PossExtremaInYSet;
};


CubicSpline::CubicSpline(std::vector<point2d>* ReferencePoints, double Tmin /*= DEFAULT_TMIN*/, double Tmax /*= DEFAULT_TMAX*/) {
    T_min = Tmin;
    T_max = Tmax;
    NumberOfPoints = (int)ReferencePoints->size();
    RefPoints = new std::vector<point2d>(*ReferencePoints);
    Ti = new std::vector<double>(NumberOfPoints);
    XSecondDerivative = new std::vector<double>(NumberOfPoints);
    YSecondDerivative = new std::vector<double>(NumberOfPoints);
    XCoeff = new std::vector<coefficients>(NumberOfPoints);
    YCoeff = new std::vector<coefficients>(NumberOfPoints);
    PossExtremaInXSet = false;
    PossExtremaInYSet = false;

    // ��������� �������� T. � ������ ����������� ����� ���������� [T_i-1, T_i] ��������������� ������ ��������������� �������� ����� ������� �� ��������� Oxy.
    // �����, ����� ���������� ������������� ����������� �� ����������, t ����� ������� ����� "�����������" �������� renormalize()
    Ti->at(0) = 0;
    for(int i = 1; i < NumberOfPoints; ++i) {
        Ti->at(i) = Ti->at(i-1) + sqrt(
                                       ((RefPoints->at(i)).x - (RefPoints->at(i-1)).x)*((RefPoints->at(i)).x - (RefPoints->at(i-1)).x) +
                                       ((RefPoints->at(i)).y - (RefPoints->at(i-1)).y)*((RefPoints->at(i)).y - (RefPoints->at(i-1)).y)
                                       );
    };
    // ������� ����������� �������� �� ������� [0; TempLength] � �������� �� [Tmin; Tmax]
    double TempLength = Ti->at(NumberOfPoints - 1);
    for(int i=0; i < NumberOfPoints; i++) {
        Ti->at(i) = Ti->at(i)*(T_max-T_min)/TempLength + T_min;
    };
    //FIXIT! renormalize();
    // ��������� �������� ������������� ����������� �_i(t)  � y_i(t)
    SetCoefficients();
}

int CubicSpline::renormalize(int density/* = DEFAULT_DENSITY*/) {
    std::vector<double> SegmentLength(NumberOfPoints-1);
    double TotalLength = 0;
    for (int i = 0; i != NumberOfPoints; ++i) {
        SegmentLength[i] = 0;
        for (int j = 0; j != density; ++j) {
            SegmentLength[i] += sqrt(
                                    (X(Ti->at(i)+(Ti->at(i+1)-Ti->at(i))*(j+1)/density) - X(Ti->at(i)+(Ti->at(i+1)-Ti->at(i))*j/density))*
                                     (X(Ti->at(i)+(Ti->at(i+1)-Ti->at(i))*(j+1)/density) - X(Ti->at(i)+(Ti->at(i+1)-Ti->at(i))*j/density)) +
                                    (Y(Ti->at(i)+(Ti->at(i+1)-Ti->at(i))*(j+1)/density) - Y(Ti->at(i)+(Ti->at(i+1)-Ti->at(i))*j/density))*
                                     (Y(Ti->at(i)+(Ti->at(i+1)-Ti->at(i))*(j+1))/density - Y(Ti->at(i)+(Ti->at(i+1)-Ti->at(i))*j/density))
                                     );

        }
        TotalLength += SegmentLength[i];
    }



    return 0;
}

int CubicSpline::test() {

return 0;
}


void CubicSpline::SetCoefficients(int BoundaryConditions/* = DEFAULT_BOUNDARY_CONDITIONS*/) {
    TridiagonalMatrix TridiagM;
    TridiagM.LowerDiagonal = new std::vector<double>(NumberOfPoints);
    TridiagM.MainDiagonal = new std::vector<double>(NumberOfPoints);
    TridiagM.UpperDiagonal = new std::vector<double>(NumberOfPoints);
    std::vector<double>* RightSideForX = new std::vector<double>(NumberOfPoints);
    std::vector<double>* RightSideForY = new std::vector<double>(NumberOfPoints);

    // ������� �������� ������� � ������ ��������� ����� ���� �� ���������� �������, ��������� ������� �������� ������������ ���������� ���������
    TridiagM.LowerDiagonal->at(0) = TridiagM.LowerDiagonal->at(NumberOfPoints-1) = 0;
    TridiagM.UpperDiagonal->at(0) = TridiagM.UpperDiagonal->at(NumberOfPoints-1) = 0;

    // �������� BoundaryConditions �� ������ �������������� ��� ����� switch
    switch(BoundaryConditions) {
        case DEFAULT_BOUNDARY_CONDITIONS:
        default:
            {
            TridiagM.MainDiagonal->at(0) = TridiagM.MainDiagonal->at(NumberOfPoints-1) = 1;
            RightSideForX->at(0) = RightSideForX->at(NumberOfPoints-1) = 0;
            RightSideForY->at(0) = RightSideForY->at(NumberOfPoints-1) = 0;
            }
            break;
    }

    // ��������� �������� ������� � ������ ������ ���������, ��������� ���������
    for(int i=1; i < NumberOfPoints-1; ++i) {
        TridiagM.LowerDiagonal->at(i) = Ti->at(i) - Ti->at(i-1);
        TridiagM.UpperDiagonal->at(i) = Ti->at(i+1) - Ti->at(i);
        TridiagM.MainDiagonal->at(i) = (Ti->at(i+1) - Ti->at(i-1))*2;
        RightSideForX->at(i) = 6*(((RefPoints->at(i+1)).x - (RefPoints->at(i)).x)/(TridiagM.UpperDiagonal->at(i)) - ((RefPoints->at(i)).x - (RefPoints->at(i-1)).x)/(TridiagM.LowerDiagonal->at(i)));
        RightSideForY->at(i) = 6*(((RefPoints->at(i+1)).y - (RefPoints->at(i)).y)/(TridiagM.UpperDiagonal->at(i)) - ((RefPoints->at(i)).y - (RefPoints->at(i-1)).y)/(TridiagM.LowerDiagonal->at(i)));
    }

    // ������� � ������ ����� ���������, ������ ������� ������������ ������ ����������� X � Y
    SolveTridiagonalSystem(&TridiagM, XSecondDerivative, RightSideForX);
    SolveTridiagonalSystem(&TridiagM, YSecondDerivative, RightSideForY);

    // ������ ����������� �������, ��������� �������� ������������� ����������� �� ������ ���������
    double h;
    for(int i=0; i!=NumberOfPoints-1; ++i)
    {
        h = Ti->at(i+1) - Ti->at(i);

        (XCoeff->at(i+1)).A = (XSecondDerivative->at(i+1) - XSecondDerivative->at(i))/(6*h);
        (XCoeff->at(i+1)).B = (XSecondDerivative->at(i))/2;
        (XCoeff->at(i+1)).C = ((RefPoints->at(i+1)).x - (RefPoints->at(i)).x)/h - (XSecondDerivative->at(i+1))*h/6 - (XSecondDerivative->at(i))*h/3;
        (XCoeff->at(i+1)).D = ((RefPoints->at(i)).x);

        (YCoeff->at(i+1)).A = (YSecondDerivative->at(i+1) - YSecondDerivative->at(i))/(6*h);
        (YCoeff->at(i+1)).B = (YSecondDerivative->at(i))/2;
        (YCoeff->at(i+1)).C = ((RefPoints->at(i+1)).y - (RefPoints->at(i)).y)/h - (YSecondDerivative->at(i+1))*h/6 - (YSecondDerivative->at(i))*h/3;
        (YCoeff->at(i+1)).D = ((RefPoints->at(i)).y);

    }


}


void CubicSpline::SolveTridiagonalSystem(TridiagonalMatrix* matrix, std::vector<double>* var, std::vector<double>* RightSide) {

    // ����� ��� ��������: a[i]*var[i-1] + b[i]*var[i] + c[i]*var[i+1] = F[i], i = 0, ..., N-1

    std::vector<double> P(NumberOfPoints), Q(NumberOfPoints), gamma(NumberOfPoints);

    gamma.at(0) = matrix->MainDiagonal->at(0);
    P.at(0)=-(matrix->UpperDiagonal->at(0))/gamma.at(0);
    Q.at(0)=(RightSide->at(0))/(gamma.at(0));
    for(int i=1; i<NumberOfPoints; ++i)
    {
        gamma.at(i) = matrix->MainDiagonal->at(i) + (matrix->LowerDiagonal->at(i))*P.at(i-1);
        P.at(i) = -(matrix->UpperDiagonal->at(i))/gamma.at(i);
        Q.at(i) = (RightSide->at(i) - (matrix->LowerDiagonal->at(i))*Q.at(i-1))/gamma.at(i);
    }

    var->at(NumberOfPoints-1) = Q.at(NumberOfPoints-1);
    for(int i = NumberOfPoints-2; i>=0; --i)
        var->at(i) = P.at(i)*var->at(i+1) + Q.at(i);
}

int CubicSpline::RangeNumber(double T) {
    if (T == T_min) return 1;
    if (T == T_max) return NumberOfPoints-1;
    int i=1;
    while(Ti->at(i)<T) i++;
    return i;
}

double CubicSpline::X(double T) {
    int i = RangeNumber(T);
    double delta = T - Ti->at(i-1);
    return (XCoeff->at(i)).A * delta * delta *delta + (XCoeff->at(i)).B * delta * delta + (XCoeff->at(i)).C * delta + (XCoeff->at(i)).D;
}


double CubicSpline::Y(double T) {
    int i = RangeNumber(T);
    double delta = T - Ti->at(i-1);
    return (YCoeff->at(i)).A * delta * delta * delta + (YCoeff->at(i)).B * delta * delta + (YCoeff->at(i)).C * delta + (YCoeff->at(i)).D;
}

double CubicSpline::Tmin() {
    return T_min;
}

double CubicSpline::Tmax() {
    return T_max;
}


void CubicSpline::InitializeExtremaList(int type) {
    std::vector<double>* PossibleExtremaList;
    std::vector<coefficients>* CoefficientsList;

    switch (type) {
        case IN_X:
            PossibleExtremaList = &PossExtremaInX;
            CoefficientsList = XCoeff;
            break;
        case IN_Y:
            PossibleExtremaList = &PossExtremaInY;
            CoefficientsList = YCoeff;
            break;
    }

    const int imin = this->RangeNumber(T_min);
    const int imax = this->RangeNumber(T_max);
    double DD = 0;
    double CurrentT = 0, CurrentT2 = 0, SwapTemp = 0;
    int icurr = 0;

    PossibleExtremaList->clear();
    // �������� ��������� ������ "��������������" �����.
    PossibleExtremaList->push_back(T_min);


    // ������� ������������ �������� � T_min
    icurr = imin;
    DD = 4 * (CoefficientsList->at(icurr)).B * (CoefficientsList->at(icurr)).B - 12 * (CoefficientsList->at(icurr)).A * (CoefficientsList->at(icurr)).C;
    if (DD == 0) {
        CurrentT = -2*((CoefficientsList->at(icurr)).B)/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
        if ((T_min < CurrentT) && (CurrentT <= Ti->at(icurr)))
            PossibleExtremaList->push_back(CurrentT);
    }
    else if (DD > 0) {
        CurrentT = (-sqrt(DD)-2*((CoefficientsList->at(icurr)).B))/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
        CurrentT2 = (sqrt(DD)-2*((CoefficientsList->at(icurr)).B))/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
        if (CurrentT > CurrentT2) {
            SwapTemp = CurrentT;
            CurrentT = CurrentT2;
            CurrentT2 = SwapTemp;
        }
        if ((T_min < CurrentT) && (CurrentT <= Ti->at(icurr)))
            PossibleExtremaList->push_back(CurrentT);
        if ((T_min < CurrentT2) && (CurrentT2 <= Ti->at(icurr)))
            PossibleExtremaList->push_back(CurrentT2);
    }
    // ����� ������������� ���������
    for(icurr = imin + 1; icurr != imax; ++icurr) {
        DD = 4 * (CoefficientsList->at(icurr)).B * (CoefficientsList->at(icurr)).B - 12 * (CoefficientsList->at(icurr)).A * (CoefficientsList->at(icurr)).C;
        if (DD == 0) {
            CurrentT = -2*((CoefficientsList->at(icurr)).B)/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
            if ((Ti->at(icurr-1) < CurrentT) && (CurrentT <= Ti->at(icurr)))
                PossibleExtremaList->push_back(CurrentT);
            }
        else if (DD > 0) {
            CurrentT = (-sqrt(DD)-2*((CoefficientsList->at(icurr)).B))/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
            CurrentT2 = (sqrt(DD)-2*((CoefficientsList->at(icurr)).B))/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
            if (CurrentT > CurrentT2) {
            SwapTemp = CurrentT;
            CurrentT = CurrentT2;
            CurrentT2 = SwapTemp;
        }

            if ((Ti->at(icurr-1) < CurrentT) && (CurrentT <= Ti->at(icurr)))
                PossibleExtremaList->push_back(CurrentT);
            if ((Ti->at(icurr-1) < CurrentT2) && (CurrentT2 <= Ti->at(icurr)))
                PossibleExtremaList->push_back(CurrentT2);
        }
    }
    // � ������� - �������� � T_max
    icurr = imax;
    DD = 4 * (CoefficientsList->at(icurr)).B * (CoefficientsList->at(icurr)).B - 12 * (CoefficientsList->at(icurr)).A * (CoefficientsList->at(icurr)).C;
    if (DD == 0) {
        CurrentT = -2*((CoefficientsList->at(icurr)).B)/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
        if ((Ti->at(icurr-1) < CurrentT) && (CurrentT < T_max))
            PossibleExtremaList->push_back(CurrentT);
    }
    else if (DD > 0) {
        CurrentT = (-sqrt(DD)-2*((CoefficientsList->at(icurr)).B))/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
        CurrentT2 = (sqrt(DD)-2*((CoefficientsList->at(icurr)).B))/(6*(CoefficientsList->at(icurr)).A) + Ti->at(icurr-1);
        if (CurrentT > CurrentT2) {
            SwapTemp = CurrentT;
            CurrentT = CurrentT2;
            CurrentT2 = SwapTemp;
        }
        if ((Ti->at(icurr-1) < CurrentT) && (CurrentT < T_max))
            PossibleExtremaList->push_back(CurrentT);

        if ((Ti->at(icurr-1) < CurrentT2) && (CurrentT2 < T_max))
            PossibleExtremaList->push_back(CurrentT2);
    }

    // ����������� ��������� ������ "��������������" ����� ������ �������� ���������
    PossibleExtremaList->push_back(T_max);

    switch (type) {
    case IN_X:
        PossExtremaInXSet = true;
        break;
    case IN_Y:
        PossExtremaInYSet = true;
        break;
    }
}

parametricvalue CubicSpline::ExtremalValueInRange(int RangeType, int ExtremaType, double Tleft, double Tright) {
    parametricvalue Result;
    std::vector<double>* PossibleExtremaList;
    typedef double (CubicSpline::*ValueFunctionPointer)(double t);
    ValueFunctionPointer GetValue;

    switch (RangeType) {
        case IN_X:
            if (!PossExtremaInXSet)
                InitializeExtremaList(IN_X);
            PossibleExtremaList = &PossExtremaInX;
            GetValue = &(this->X);
            break;
        case IN_Y:
            if (!PossExtremaInYSet)
                InitializeExtremaList(IN_Y);
            PossibleExtremaList = &PossExtremaInY;
            GetValue = &(this->Y);
            break;
    }

    // ������������� �������� ��������� ��� � ������� ����� ���������, ��� � ����� �� �������������, ������� ������.
    // ��������� �� �� ������� (������, ������������ InitializeExtremaList, ������������ �� �����������)
    { // ��������� ������� ������������� LeftResult � RightResult
    double LeftResult = (*this.*GetValue)(Tleft);
    double RightResult = (*this.*GetValue)(Tright);
    if ((ExtremaType == MAX_VALUE && LeftResult > RightResult) || (ExtremaType == MIN_VALUE && LeftResult < RightResult)) {
        Result.t = Tleft;
        Result.value = LeftResult;

    } else {
        Result.t = Tright;
        Result.value = RightResult;
    }
    }

    int i;
    for(i = 0; PossibleExtremaList->at(i) <= Tleft; ++i);
    double CurrentT=0, CurrentValue=0;
    int ExtremaListSize = (int)PossibleExtremaList->size();
    for(;(PossibleExtremaList->at(i) < Tright) && (i != ExtremaListSize - 1); ++i) {
        CurrentT = PossibleExtremaList->at(i);
        CurrentValue = (*this.*GetValue)(CurrentT);
        if ((ExtremaType == MAX_VALUE && CurrentValue > Result.value) || (ExtremaType == MIN_VALUE && CurrentValue < Result.value)) {
            Result.t = CurrentT;
            Result.value = CurrentValue;
        }
    }
    return Result;
}


parametricvalue CubicSpline::MinXInRange(double Tleft, double Tright) {
    return ExtremalValueInRange(IN_X, MIN_VALUE, Tleft, Tright);
}

parametricvalue CubicSpline::MinYInRange(double Tleft, double Tright) {
    return ExtremalValueInRange(IN_Y, MIN_VALUE, Tleft, Tright);
}

parametricvalue CubicSpline::MaxXInRange(double Tleft, double Tright) {
    return ExtremalValueInRange(IN_X, MAX_VALUE, Tleft, Tright);
}

parametricvalue CubicSpline::MaxYInRange(double Tleft, double Tright) {
    return ExtremalValueInRange(IN_Y, MAX_VALUE, Tleft, Tright);
}


#endif // CUBIC_SPLINES_H
