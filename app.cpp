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

    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = pow(idxMinusMean[idx], 2);
        sum += results[idx];
    }

    return sum;
}

long double calculateXiMinusXBarTimesYiMinusYBar(int popSize, long double *xiMinusXBar, long double *yiMinusYBar, long double *results)
{
    long double sum = 0;

    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = (xiMinusXBar[idx] * yiMinusYBar[idx]);
        sum += results[idx];
    }

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
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = intercept + (slope * xValues[idx]);
    }
}

void calculateYiMinusYHati(int popSize, long double *yValues, long double *yHatValues, long double *results)
{

    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = yValues[idx] - yHatValues[idx];
    }
}
// could just use idx-mean squared for
long double calculateYiMinusYHatiSquaredAndSum(int popSize, long double *yiMinusYHati, long double *results)
{
    long double sum = 0;
    for (int idx = 0; idx < popSize; idx++)
    {
        results[idx] = pow(yiMinusYHati[idx], 2);

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
    long double MSE = 0, MSERoot = 0;
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

    // calculating covariance; sickly enough, we can just use the variance function with the (xi-xbar)(yi-ybar) val and it works
    covariance = calculateVariance(popSize, sigmaXiMinusXBarTimesYiMinusYBar, 1);

    // calculating correlation coefficient
    correlationCoefficient = calculateCorrelationCoefficient(covariance, xStdDev, yStdDev);

    // o h  b o y  l i n e a r  r e g r e s s i o n  t i m e

    // calculating slope
    slope = calculateSlope(sigmaXiMinusXBarTimesYiMinusYBar, sigmaXiMinusXBarSquared);

    // calculating intercept

    intercept = calculateIntercept(yBar, xBar, slope);

    // big stupid ass yhat table time now for SSE, SSR and SST :)))))))
    yHatValues = new long double[popSize];

    calculateYHatValues(popSize, xValues, intercept, slope, yHatValues);

    // calculating yi - y^i
    yiMinusYHati = new long double[popSize];
    calculateYiMinusYHati(popSize, yValues, yHatValues, yiMinusYHati);

    // calculating (yi - y^i)^2 and SSE
    yiMinusYHatiSquared = new long double[popSize];
    SSE = calculateYiMinusYHatiSquaredAndSum(popSize, yiMinusYHati, yiMinusYHatiSquared);

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
    MSERoot = calculateStdDev(MSE);

    // Slope Standard Deviation time
    slopeStdDev = calculateSlopeStdDev(MSERoot, sigmaXiMinusXBarSquared);

    // at last, test statistic
    tTestStatistic = calculateTestStatistic(slope, slopeStdDev);

    // header for dis shizzle
    cout << left << setprecision(10) << setw(10) << "xi" << setw(10) << "yi"
         << setw(16) << "xi - xbar" << setw(16) << "yi - ybar"
         << setw(20) << "(xi - xbar)^2" << setw(20) << "(yi - ybar)^2" << setw(24) << "(xi - xbar)(yi - ybar)"
         << setw(16) << "y^i"
         << setw(16) << "yi - y^i" << setw(20) << "(yi - y^i)^2"
         << setw(16) << "y^i - ybar" << setw(20) << "(y^i - ybar)^2" << endl;

    // printing the table
    for (int idx = 0; idx < popSize; idx++)
    {
        cout << left << setw(10) << xValues[idx] << setw(10) << yValues[idx]
             << setw(16) << xiMinusXBar[idx] << setw(16) << yiMinusYBar[idx]
             << setw(20) << xiMinusXBarSquared[idx] << setw(20) << yiMinusYBarSquared[idx] << setw(24) << xiMinusXBarTimesYiMinusYBar[idx]
             << setw(16) << yHatValues[idx]
             << setw(16) << yiMinusYHati[idx] << setw(20) << yiMinusYHatiSquared[idx]
             << setw(16) << yHatMinusYBar[idx] << setw(20) << yHatMinusYBarSquared[idx] << endl;
    }

    cout << left << setw(52) << "Sums: " // this is to clear the spaces, hopefully it works
         << setw(20) << sigmaXiMinusXBarSquared << setw(20) << sigmaYiMinusYBarSquared << setw(24) << sigmaXiMinusXBarTimesYiMinusYBar
         << setw(32) << " "
         << setw(20) << SSE
         << setw(16) << " "
         << setw(20) << SSR << endl;

    // print out other important values, including;
    //  xbar, ybar, sigma(xi-xbar)^2, sigma(yi-ybar)^2, sigma(xi-xbar)(yi-ybar),
    //  (Sx2) xVariance, (Sx) xStdDev, (Sy2) yVariance, (Sy) yStdDev, (Sxy) covariance
    //  correlation coefficient (Rxy)
    //  intercept (b1), slope (b0), yHat equation (b1 + (b0 * x))
    //  SSE, SSR, SST, coefficient of determination (r2)
    //  MSE (S2), MSERoot (S), slopeStdDev (Sb1)
    // test statistic (t)

    cout << endl
         << endl;

    cout << left << setprecision(8) << setw(8) << "x-bar: " << setw(12) << xBar << setw(8) << "y-bar: " << setw(12) << yBar << endl
         << setw(8) << "Sx2:" << setw(12) << xVariance << setw(8) << "Sy2:" << setw(12) << yVariance << endl
         << setw(8) << "Sx:" << setw(12) << xStdDev << setw(8) << "Sx:" << setw(12) << xStdDev << endl
         << setw(8) << "Sxy:" << setw(12) << covariance << setw(8) << "Rxy:" << setw(12) << correlationCoefficient << endl
         << setw(8) << "b1:" << setw(12) << intercept << setw(8) << "b0:" << setw(12) << slope << endl
         << setw(8) << "SSE:" << setw(12) << SSE << setw(8) << "SSR:" << setw(12) << SSR << endl
         << setw(8) << "SST:" << setw(12) << sigmaYiMinusYBarSquared << setw(8) << "R2:" << setw(12) << coefficientOfDetermination << endl
         << setw(8) << "MSE/S2:" << setw(12) << MSE << setw(8) << "S:" << setw(12) << MSERoot << endl
         << setw(8) << "Sb1:" << setw(12) << slopeStdDev << setw(8) << "T Stat:" << setw(12) << tTestStatistic << endl;

    cout << endl
         << endl;

    cout << "Y-Hat equation: " << setprecision(8) << intercept << " * (" << slope << " * x)" << endl;

    delete xValues;
    delete yValues;
    delete xiMinusXBar;
    delete yiMinusYBar;
    delete xiMinusXBarSquared;
    delete yiMinusYBarSquared;
    delete xiMinusXBarTimesYiMinusYBar;
    delete yHatValues;
    delete yiMinusYHati;
    delete yiMinusYHatiSquared;
    delete yHatMinusYBar;
    delete yHatMinusYBarSquared;
}

int main()
{
    linearRegression();
    return 0;
}
