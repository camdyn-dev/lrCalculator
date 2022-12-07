#include <iostream>
#include <cmath>

using namespace std;

double calculateMean(int popSize, double *pop)
{
    double sum = 0;
    double mean = 0;
    for (int idx = 0; idx < popSize; idx++)
    {
        sum += pop[idx];
    }
    mean = sum / popSize;

    return mean;
}

double linearRegression()
{
    // holders for important values
    double xBar;
    double yBar;

    int popSize;
    double *xValues, *yValues;

    cout << "How many subjects are in your population?: ";
    cin >> popSize;

    // creating arrays for population values for calculations
    xValues = new double[popSize];
    yValues = new double[popSize];

    // receiving values for x and y population
    for (int xi = 0; xi < popSize; xi++)
    {
        cout << "What's xi" << xi + 1 << "?: ";
        cin >> xValues[xi];
    }

    cout << endl;
    for (int yi = 0; yi < popSize; yi++)
    {
        cout << "What's yi" << yi + 1 << "?: ";
        cin >> yValues[yi];
    }

    // calculation time
    xBar = calculateMean(popSize, xValues);
    yBar = calculateMean(popSize, yValues);

    cout << "xbar = " << xBar << '\t'
         << "ybar = " << yBar << endl;

    return 0;
}

int main()
{
    linearRegression();
    return 0;
}
