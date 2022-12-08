#include <iostream>
#include <cmath>

using namespace std;

// these names boutta be verbose as FUCK

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

double calculateIndexMinusMeanSquaredAndSum(int popSize, double *idxMinusMean, double *results)
{
    // used to calculate (xi - xbar)^2 and (yi-ybar)^2
    // using double instead of void to return the sigma of (xi-xbar)^2 and (yi-ybar)^2
    double sum = 0;
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = pow(idxMinusMean[idx], 2);
        sum += results[idx];
        cout << "result: " << results[idx] << endl;
    }
    cout << "sum: " << sum << endl;
    return sum;
}

double calculateXiMinusXBarTimesYiMinusYBar(int popSize, double *xiMinusXBar, double *yiMinusYBar, double *results)
{
    double sum = 0;
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = (xiMinusXBar[idx] * yiMinusYBar[idx]);
        sum += results[idx];
        cout << "result: " << results[idx] << endl;
    }
    cout << "sum: " << sum << endl;
    return sum;
};

double linearRegression()
{
    // holders for important values
    // initializing all with zero just to be safe. They'll either be re-assigned or incremented/decremented
    double xBar = 0;
    double yBar = 0;

    int popSize = 0;
    double *xValues, *yValues;
    double *xiMinusXBar, *yiMinusYBar;

    double *xiMinusXBarSquared, *yiMinusYBarSquared;
    double sigmaXiMinusXBarSquared = 0, sigmaYiMinusYBarSquared = 0;

    // these variable names are getting fucking large
    double *xiMinusXBarTimesYiMinusYBar;
    double sigmaXiMinusXBarTimesYiMinusYBar = 0;

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
    yiMinusYBarSquared = new double[popSize];

    // calculating values for (xi - xbar)^2 and (yi-ybar)^2
    // described in function, but we're modifying arrays inside since pointers, and returning the sum
    // I guess I could use a pass-by reference to be more consistent? this seems safer, though.
    sigmaXiMinusXBarSquared = calculateIndexMinusMeanSquaredAndSum(popSize, xiMinusXBar, xiMinusXBarSquared);
    sigmaYiMinusYBarSquared = calculateIndexMinusMeanSquaredAndSum(popSize, yiMinusYBar, yiMinusYBarSquared);

    cout << "sigma of (xi - xbar)^2 : " << sigmaXiMinusXBarSquared << endl;
    cout << "sigma of (yi - ybar)^2 : " << sigmaYiMinusYBarSquared << endl;

    // creating array for (xi-xbar)(yi-ybar)
    xiMinusXBarTimesYiMinusYBar = new double[popSize];
    sigmaXiMinusXBarTimesYiMinusYBar = calculateXiMinusXBarTimesYiMinusYBar(popSize, xiMinusXBar, yiMinusYBar, xiMinusXBarTimesYiMinusYBar);

    return 0;
}

int main()
{
    linearRegression();
    return 0;
}
