#include <iomanip>
#include "common.h"
using namespace std;

#define DIG 7 // number of digits after the decimal point

void readEdges()
{
    std::ifstream input("edgeTypes.txt");
    for (i = 0; i < 4; i++)
        input >> edges[i];
    input.close();
}

void readSize()
{
    std::ifstream input("size.txt");
    input >> nX >> nY;
    input.close();

    N = nX * nY;
}

void readGrid()
{
    int i;

    ifstream inputX("nodesX.txt");
    for (i = 0; i < nX; i++)
        inputX >> gridX[i];
    inputX.close();

    ifstream inputY("nodesY.txt");
    for (i = 0; i < nY; i++)
        inputY >> gridY[i];
    inputY.close();
}

void printResult()
{
    int n = 0;
    double tmpNormErr = 0.0, tmpNormAn = 0.0;
    double xPoint, yPoint, analytical, error, relNorm;

    ofstream output("result.txt");

    output << setprecision(DIG) << setw(DIG + 10)
        << "Numerical" << setw(DIG + 20)
        << "Analytical" << setw(DIG + 20)
        << "Error" << setw(DIG + 17)
        << "x     y\n"
    << endl;

    for (k = 0; k < nY; k++)
        for (i = 0; i < nX; i++)
        {
            xPoint = gridX[i];
            yPoint = gridY[k];
            analytical = getAnalytical(xPoint, yPoint);
            error = std::abs(x[n] - analytical);

            output
                << setw(DIG + 10) << x[n]
                << setw(DIG + 20) << analytical
                << setw(DIG + 20) << error
                << setw(DIG + 10) << xPoint
                << setw(6) << yPoint
            << endl;

            tmpNormErr += error * error;
            tmpNormAn += analytical * analytical;

            n++;
        }

    relNorm = sqrt(tmpNormErr / tmpNormAn);

    output << setprecision(5) << "\n" << setw(DIG + 43) <<
        "||u - u*|| / ||u*|| = " << relNorm
    << endl;

    output.close();
}
