#include <fstream>
#include "common.h"

void readSize() // reading area size
{
    std::ifstream input("size.txt");
    input >> nX >> nY;
    N = nX * nY;
    input.close();
}

void readGrid() // reading gridX and gridY
{
    int i;
    std::ifstream inputX("gridX.txt");
    for (i = 0; i < nX; i++)
        inputX >> gridX[i];
    inputX.close();

    std::ifstream inputY("gridY.txt");
    for (i = 0; i < nY; i++)
        inputY >> gridY[i];
    inputY.close();
}