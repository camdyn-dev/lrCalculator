#include <iostream>
#include <cmath>

using namespace std;

double calculateMean(int popSize, double *pop)
{
    // used to calculate xbar and ybar
    double sum = 0;
    double mean = 0;
    for (int idx = 0; idx < popSize; idx++)
    {
        sum += pop[idx];
    }
    mean = sum / popSize;

    return mean;
}

void calculateIndexMinusMean(int popSize, double mean, double *pop, double *results)
{
    // used to calculate xi - xbar and yi - ybar
    // apparently I can't return a pointer, so I'll assign the values to the new function directly in here
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = (pop[idx] - mean);
    }
}

double linearRegression()
{
    // holders for important values
    double xBar;
    double yBar;

    int popSize;
    double *xValues, *yValues;
    double *xiMinusXBar, *yiMinusYBar;

    double *xiMinusXBarSquared, *yiMinusXBarSquared;

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

    // creating arrays for xi-xbar and yi-ybar
    xiMinusXBar = new double[popSize];
    yiMinusYBar = new double[popSize];

    // calculating values for xi-xbar and yi-ybar
    calculateIndexMinusMean(popSize, xBar, xValues, xiMinusXBar);
    calculateIndexMinusMean(popSize, yBar, yValues, yiMinusYBar);

    // creating arrays for (xi - xbar)^2 and (yi - ybar)^2
    xiMinusXBarSquared = new double[popSize];
    yiMinusXBarSquared = new double[popSize];

    // calculating values for (xi - xbar)^2 and (yi-ybar)^2

    return 0;
}

int main()
{
    linearRegression();
    return 0;
}
