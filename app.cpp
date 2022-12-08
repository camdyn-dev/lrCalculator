#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

/*
TODOS:
everything is a long double for precision i guess
fuck memory lol
*/

/*
TEST DATA
X values: 180, 150, 95, 70, 70, 35
Y values: 76, 71, 62, 57, 30, 34
*/

long double calculateMean(int popSize, long double *pop)
{
    // used to calculate xbar and ybar
    long double sum = 0;
    long double mean = 0;
    for (int idx = 0; idx < popSize; idx++)
    {
        sum += pop[idx];
    }
    mean = sum / popSize;

    return mean;
}

void calculateIndexMinusMean(int popSize, long double mean, long double *pop, long double *results)
{
    // used to calculate xi - xbar and yi - ybar
    // apparently I can't return a pointer, so I'll assign the values to the new function directly in here
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = (pop[idx] - mean);
    }
}

long double calculateIndexMinusMeanSquaredAndSum(int popSize, long double *idxMinusMean, long double *results)
{
    // used to calculate (xi - xbar)^2 and (yi-ybar)^2
    // using double instead of void to return the sigma of (xi-xbar)^2 and (yi-ybar)^2
    long double sum = 0;
    cout << "(idx - mean)^2" << endl;
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = pow(idxMinusMean[idx], 2);
        sum += results[idx];
        cout << results[idx] << endl;
    }
    cout << "sum: " << sum << endl;
    return sum;
}

long double calculateXiMinusXBarTimesYiMinusYBar(int popSize, long double *xiMinusXBar, long double *yiMinusYBar, long double *results)
{
    long double sum = 0;
    cout << "(xi - xbar)(yi - ybar)" << endl;
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = (xiMinusXBar[idx] * yiMinusYBar[idx]);
        sum += results[idx];
        cout << results[idx] << endl;
    }
    cout << "sum: " << sum << endl;
    return sum;
};

long double calculateVariance(int popSize, long double sigmaIdxMinusMeanSqaured, int decrement)
{
    // decrement is used for degrees of freedom
    //  for example, on standard variances of one set, it will be n -1
    //  but, when calculating slope variance (MSE), or for multiple sets, it'll be n - (number of sets)
    long double variance = 0;
    variance = sigmaIdxMinusMeanSqaured / (popSize - decrement);
    return variance;
}

long double calculateStdDev(long double variance)
{
    // using sqrtl and long double to ensure i get as much precision as possible
    // might wanna change most numbers to long doubles
    return sqrtl(variance);
}

long double calculateCorrelationCoefficient(long double covariance, long double xStdDev, long double yStdDev)
{
    long double bottom = xStdDev * yStdDev;
    long double correlationCoefficient = covariance / bottom;
    return correlationCoefficient;
}

long double calculateSlope(long double sigmaXiMinusXBarTimesYiMinusYBar, long double sigmaXiMinusXBarSquared)
{
    return (sigmaXiMinusXBarTimesYiMinusYBar / sigmaXiMinusXBarSquared);
}

long double calculateIntercept(long double yBar, long double xBar, long double slope)
{
    return (yBar - (slope * xBar));
}

void calculateYHatValues(int popSize, long double *xValues, long double intercept, long double slope, long double *results)
{
    cout << "Y Hat I values: " << endl;
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = intercept + (slope * xValues[idx]);
        cout << results[idx] << endl;
    }
}

void calculateYiMinusYHati(int popSize, long double *yValues, long double *yHatValues, long double *results)
{
    cout << "yi - y^i values: " << endl;
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = yValues[idx] - yHatValues[idx];
        cout << results[idx] << endl;
    }
}
// could just use idx-mean squared for
long double calculateYiMinusYHatiSquaredAndSum(int popSize, long double *yiMinusYHati, long double *results)
{
    long double sum = 0;
    cout << "(yi - y^i)^2 values: " << endl;
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = pow(yiMinusYHati[idx], 2);
        cout << results[idx] << endl;
        sum += results[idx];
    }
    return sum;
}

long double calculateTestStatistic(long double slope, long double slopeStdDev)
{
    return slope / slopeStdDev;
}

long double calculateSlopeStdDev(long double MSESquared, long double sigmaXiMinusXbarSquared)
{
    long double slopeStdDev = MSESquared / sqrt(sigmaXiMinusXbarSquared);
    return slopeStdDev;
}

void linearRegression()
{
    // holders for important values
    // initializing all with zero just to be safe. They'll either be re-assigned or incremented/decremented
    double xBar = 0;
    double yBar = 0;

    int popSize = 0;
    long double *xValues, *yValues;
    long double *xiMinusXBar, *yiMinusYBar;

    long double *xiMinusXBarSquared, *yiMinusYBarSquared;
    long double sigmaXiMinusXBarSquared = 0, sigmaYiMinusYBarSquared = 0;

    // these variable names are getting fucking large
    long double *xiMinusXBarTimesYiMinusYBar;
    long double sigmaXiMinusXBarTimesYiMinusYBar = 0;

    // variances, standard deviations, covariance
    // using long double for more accuracy
    long double xVariance = 0, yVariance = 0;
    long double xStdDev = 0, yStdDev = 0;
    long double covariance = 0;

    // correlation coefficienct
    long double correlationCoefficient = 0;

    // linear regression numbers
    long double intercept = 0, slope = 0, error = 0;

    // pointer arrays for yhat values
    long double *yHatValues;
    long double *yiMinusYHati, *yiMinusYHatiSquared;
    long double *yHatMinusYBar, *yHatMinusYBarSquared;

    // SSE, SSR and SST from yhats
    // SSE = sigma(yi - yHat)^2
    // SSR = sigma(yHati - ybar)^2
    // SST = SSE + SSR
    long double SSE = 0, SSR = 0, SST = 0;

    // coefficient of determination, MSE, MSESquared, slopeStdDev
    long double coefficientOfDetermination = 0;
    long double MSE = 0, MSESquared = 0;
    long double slopeStdDev = 0;

    // t test statistic POGGERRRRSSSS
    long double tTestStatistic = 0;

    //
    //
    //
    // PROGRAM START
    //
    //
    //

    cout << "How many subjects are in your population?: ";
    cin >> popSize;

    // creating arrays for population values for calculations
    xValues = new long double[popSize];
    yValues = new long double[popSize];

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
    xiMinusXBar = new long double[popSize];
    yiMinusYBar = new long double[popSize];

    // calculating values for xi-xbar and yi-ybar
    calculateIndexMinusMean(popSize, xBar, xValues, xiMinusXBar);
    calculateIndexMinusMean(popSize, yBar, yValues, yiMinusYBar);

    // creating arrays for (xi - xbar)^2 and (yi - ybar)^2
    xiMinusXBarSquared = new long double[popSize];
    yiMinusYBarSquared = new long double[popSize];

    // calculating values for (xi - xbar)^2 and (yi-ybar)^2
    // described in function, but we're modifying arrays inside since pointers, and returning the sum
    // I guess I could use a pass-by reference to be more consistent? this seems safer, though.
    sigmaXiMinusXBarSquared = calculateIndexMinusMeanSquaredAndSum(popSize, xiMinusXBar, xiMinusXBarSquared);
    sigmaYiMinusYBarSquared = calculateIndexMinusMeanSquaredAndSum(popSize, yiMinusYBar, yiMinusYBarSquared);

    cout << "sigma of (xi - xbar)^2 : " << sigmaXiMinusXBarSquared << endl;
    cout << "sigma of (yi - ybar)^2 : " << sigmaYiMinusYBarSquared << endl;

    // creating array for (xi-xbar)(yi-ybar)
    xiMinusXBarTimesYiMinusYBar = new long double[popSize];

    // calculating values for (xi-xbar)(yi-ybar)
    // same thing as above, return instead of pass-by reference
    sigmaXiMinusXBarTimesYiMinusYBar = calculateXiMinusXBarTimesYiMinusYBar(popSize, xiMinusXBar, yiMinusYBar, xiMinusXBarTimesYiMinusYBar);

    // calculating variances, standard deviation
    xVariance = calculateVariance(popSize, sigmaXiMinusXBarSquared, 1);
    yVariance = calculateVariance(popSize, sigmaYiMinusYBarSquared, 1);

    xStdDev = calculateStdDev(xVariance);
    yStdDev = calculateStdDev(yVariance);

    // outputting to check
    cout << "Sx2 = " << xVariance << endl;
    cout << "Sy2 = " << yVariance << endl;

    // using set precision to check how accurate
    cout << "Sx = " << setprecision(18) << xStdDev << endl;
    cout << "Sy = " << setprecision(18) << yStdDev << endl;

    // calculating covariance; sickly enough, we can just use the variance function with the (xi-xbar)(yi-ybar) val and it works
    covariance = calculateVariance(popSize, sigmaXiMinusXBarTimesYiMinusYBar, 1);

    cout << "Sxy = " << setprecision(18) << covariance << endl;

    // calculating correlation coefficient
    correlationCoefficient = calculateCorrelationCoefficient(covariance, xStdDev, yStdDev);

    cout << "Rxy = " << setprecision(18) << correlationCoefficient << endl;

    // o h  b o y  l i n e a r  r e g r e s s i o n  t i m e

    // calculating slope
    slope = calculateSlope(sigmaXiMinusXBarTimesYiMinusYBar, sigmaXiMinusXBarSquared);

    cout << "b1 = " << setprecision(18) << slope << endl;

    // calculating intercept

    intercept = calculateIntercept(yBar, xBar, slope);
    cout << "b0 = " << setprecision(18) << intercept << endl;

    // big stupid ass yhat table time now for SSE, SSR and SST :)))))))
    yHatValues = new long double[popSize];

    calculateYHatValues(popSize, xValues, intercept, slope, yHatValues);

    // calculating yi - y^i
    yiMinusYHati = new long double[popSize];
    calculateYiMinusYHati(popSize, yValues, yHatValues, yiMinusYHati);

    // calculating (yi - y^i)^2 and SSE
    yiMinusYHatiSquared = new long double[popSize];
    SSE = calculateYiMinusYHatiSquaredAndSum(popSize, yiMinusYHati, yiMinusYHatiSquared);

    cout << "sigma(yi - yhati)^2 OR SSE = " << SSE << endl;

    // onto SSR

    // calculating yHati - ybar
    // i can use the idx minis mean and idx minus mean squared functions for these
    yHatMinusYBar = new long double[popSize];
    calculateIndexMinusMean(popSize, yBar, yHatValues, yHatMinusYBar);

    yHatMinusYBarSquared = new long double[popSize];
    SSR = calculateIndexMinusMeanSquaredAndSum(popSize, yHatMinusYBar, yHatMinusYBarSquared);

    // okay, now that the ultra AAAAA calculations are done, time for some smaller ones
    coefficientOfDetermination = SSR / sigmaYiMinusYBarSquared;

    // MSE time
    MSE = calculateVariance(popSize, SSE, 2);
    MSESquared = calculateStdDev(MSE);

    // Slope Standard Deviation time
    slopeStdDev = calculateSlopeStdDev(MSESquared, sigmaXiMinusXBarSquared);

    // at last, test statistic
    tTestStatistic = calculateTestStatistic(slope, slopeStdDev);

    cout << "Test Statistic is: " << setprecision(3) << tTestStatistic << endl;
}

int main()
{
    linearRegression();
    return 0;
}
