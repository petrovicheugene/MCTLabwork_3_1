//=====================================================
#include <omp.h>
#include <iostream>
#include <math.h>
//=====================================================
using namespace std;
//#define Simpson

double f(double x);
double calcIntegral(int num_intervals);
double simpsonMethod(double a, double b);
double trapezeMethod(double a, double b);

const double XMin = -M_PI;
const double XMax = M_PI;

const int max_num_intervals = 4096;
//=====================================================
int main()
{
	int num_intervals = 8;
    double integral;
#ifdef Simpson
    cout << "Calculation of integral using Simpson method." << endl;

#else
    cout << "Calculation of integral using trapeze method." << endl;
#endif
    cout << "Function: pow(x, 2.0)*sin(x) / sqrt(pow(x, 2.0) + 1)" << endl;
	cout << endl << endl;

	for (; num_intervals <= max_num_intervals; num_intervals *= 2)
	{
		double startTime = omp_get_wtime();
        integral = calcIntegral(num_intervals);
		double endTime = omp_get_wtime();

        printf("Integral value: %14.16g .\n", integral);
        printf("Calculation took %14.16g sec.\n", endTime - startTime);
        cout << "Interval h: " << (XMax - XMin) / static_cast<double>(num_intervals)
			<< " Num of intervals " << num_intervals << endl << endl;
	}

	getchar();
	return 0;
}
//=====================================================
double f(double x)
{
	return pow(x, 2.0)*sin(x) / sqrt(pow(x, 2.0) + 1);
}
//=====================================================
double calcIntegral(int num_intervals)
{
	double sum = 0.0;
    double h = (XMax - XMin) / static_cast<double>(num_intervals);
#pragma omp parallel reduction (+: sum) 
	{

#pragma omp for
		for (int i = 0; i < num_intervals; i++)
		{
            double a = XMin + i * h;
			double b = a + h;
#ifdef Simpson

			sum += simpsonMethod(a, b);

#else
			sum += trapezeMethod(a, b);
#endif

		}
	}
	return sum;
}
//=====================================================
double simpsonMethod(double a, double b)
{
	return ((b - a) / 6) * (f(a) + 4.0*f((a + b) / 2.0) + f(b));
}
//=====================================================
double trapezeMethod(double a, double b)
{
	return (f(a) + f(b)) / 2.0 * (b - a);
}
//=====================================================
