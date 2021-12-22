#include "common.h"

void boundaryCond()
{
    int edgeNum;

    readEdges(); // reading condition types on each edge

    // natural conditions
    for (edgeNum = 0; edgeNum < 4; edgeNum++) // iterate over the edges: 0 - lower, 1 - right, 2 - upper, 3 - left
    {
        if (edges[edgeNum] == 2) // if second type
            secondCond(edgeNum);
        else if (edges[edgeNum] == 3) // if third type
            thirdCond(edgeNum);
    }

    // first type condition
    for (edgeNum = 0; edgeNum < 4; edgeNum++)
        if (edges[edgeNum] == 1) 
            firstCond(edgeNum);
}

void firstCond(int edgeNum)
{
    if (edgeNum == 0)
        for (i = 0; i < nX; i++)
        {
            diOrig[i] = 1.0e+100;
            f[i] = 1.0e+100 * getAnalytical(gridX[i], gridY[0]);
        }

    if (edgeNum == 1)
        for (i = 0; i < nY; i++)
        {
            diOrig[nX * (i + 1) - 1] = 1.0e+100;
            f[nX * (i + 1) - 1] = 1.0e+100 * getAnalytical(gridX[nX - 1], gridY[i]);
        }

    if (edgeNum == 2)
        for (i = 0; i < nX; i++)
        {
            diOrig[nX * (nY - 1) + i] = 1.0e+100;
            f[nX * (nY - 1) + i] = 1.0e+100 * getAnalytical(gridX[i], gridY[nY - 1]);
        }

    if (edgeNum == 3)
        for (i = 0; i < nY; i++)
        {
            diOrig[nX * (i + 1) - 3] = 1.0e+100;
            f[nX * (i + 1) - 3] = 1.0e+100 * getAnalytical(gridX[0], gridY[i]);
        }
}

void secondCond(int edgeNum)
{
    double b[2];
    int edgeNode1 = 0, edgeNode2 = 0;
    double theta1 = 0, theta2 = 0;
    double h = 0;

    if (edgeNum == 0)
    {
        edgeNode1 = 0, edgeNode2 = nX - 1;
        h = gridX[nX - 1] - gridX[0];
        theta1 = getAnalytical(gridX[0], gridY[0]), theta2 = getAnalytical(gridX[nX - 1], gridY[0]);
    }

    if (edgeNum == 1)
    {
        edgeNode1 = nX - 1, edgeNode2 = nX * nY - 1;
        h = gridY[nY - 1] - gridY[0];
        theta1 = getAnalytical(gridX[nX - 1], gridY[0]), theta2 = getAnalytical(gridX[nX - 1], gridY[nY - 1]);
    }

    if (edgeNum == 2)
    {
        edgeNode1 = nX * (nY - 1), edgeNode2 = nX * (nY - 1) + nX - 1;
        h = gridX[nX - 1] - gridX[0];
        theta1 = getAnalytical(gridX[0], gridY[nY - 1]), theta2 = getAnalytical(gridX[nX - 1], gridY[nY - 1]);
    }

    if (edgeNum == 3)
    {
        edgeNode1 = 0, edgeNode2 = nX * (nY - 1);
        h = gridY[nY - 1] - gridY[0];
        theta1 = getAnalytical(gridX[0], gridY[0]), theta2 = getAnalytical(gridX[0], gridY[nY - 1]);
    }

    b[0] = (h / 6) * (2 * theta1 + theta2);
    b[1] = (h / 6) * (theta1 + 2 * theta2);

    f[edgeNode1] += b[0];
    f[edgeNode2] += b[1];
}

void thirdCond(int edgeNum)
{
    double b[2], A[2][2];
    int edgeNode1 = 0, edgeNode2 = 0;
    double uBeta1 = 0, uBeta2 = 0;
    double h = 0;

    if (edgeNum == 0)
    {
        edgeNode1 = 0, edgeNode2 = nX - 1;
        h = gridX[nX - 1] - gridX[0];
        uBeta1 = getAnalytical(gridX[0], gridY[0]), uBeta2 = getAnalytical(gridX[nX - 1], gridY[0]);
    }

    if (edgeNum == 1)
    {
        edgeNode1 = nX - 1, edgeNode2 = nX * nY - 1;
        h = gridY[nY - 1] - gridY[0];
        uBeta1 = getAnalytical(gridX[nX - 1], gridY[0]), uBeta2 = getAnalytical(gridX[nX - 1], gridY[nY - 1]);
    }

    if (edgeNum == 2)
    {
        edgeNode1 = nX * (nY - 1), edgeNode2 = nX * (nY - 1) + nX - 1;
        h = gridX[nX - 1] - gridX[0];
        uBeta1 = getAnalytical(gridX[0], gridY[nY - 1]), uBeta2 = getAnalytical(gridX[nX - 1], gridY[nY - 1]);
    }

    if (edgeNum == 3)
    {
        edgeNode1 = 0, edgeNode2 = nX * (nY - 1);
        h = gridY[nY - 1] - gridY[0];
        uBeta1 = getAnalytical(gridX[0], gridY[0]), uBeta2 = getAnalytical(gridX[0], gridY[nY - 1]);
    }

    b[0] = getBeta() * (h / 6) * (2 * uBeta1 + uBeta2);
    b[1] = getBeta() * (h / 6) * (uBeta1 + 2 * uBeta2);

    f[edgeNode1] += b[0];
    f[edgeNode2] += b[1];

    A[0][0] = getBeta() * h / 3, A[0][1] = getBeta() * h / 6;
    A[1][0] = getBeta() * h / 6, A[1][1] = getBeta() * h / 3;

    addGlobal(A[0][0], edgeNode1, edgeNode1);
    addGlobal(A[0][1], edgeNode1, edgeNode2);
    addGlobal(A[1][0], edgeNode2, edgeNode1);
    addGlobal(A[1][1], edgeNode2, edgeNode2);
}

