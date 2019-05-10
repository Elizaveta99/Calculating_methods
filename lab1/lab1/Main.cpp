#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include<cstdio>
#include <random>

using namespace std;

const int n = 3;
const double eps = 0.00000000000000001;

int main()
{
	cin.tie(0);
	ios_base::sync_with_stdio(0);
	//freopen("input.txt", "r", stdin);
	//freopen("output.txt", "w", stdout);

	// 1
	srand(time(0));
	vector<vector<double> > a(n), a_start(n), e(n);
	vector<double> y, sum(n, 0), b(n, 0);
	for (int i = 0; i < n; i++)
	{
		double temp = 4.0 * ((double)rand() / (double)RAND_MAX) - 2.0;
		y.push_back(temp);
		a[i].resize(n, 0);
		a_start[i].resize(n, 0);
		e[i].resize(n, 0);
		e[i][i] = 1;
	}
	for (int i = 0; i < n; i++)
		for (int j = i + 1; j < n; j++)
		{
			double temp = 4.0 * ((double)rand() / (double)RAND_MAX) - 2.0;
			a[i][j] = temp;
			a[j][i] = a[i][j];
			sum[i] += a[i][j];
			sum[j] += a[j][i];
		}
	for (int i = 0; i < n; i++)
	{
		a[i][i] = sum[i];
	}

	a_start = a;

	for (int i = 0; i < n; i++) 
	{
		for (int j = 0; j < n; j++)
		{
			b[i] += a[i][j] * y[j];
		}
	}

	double norm1 = 0;
	for (int i = 0; i < n; i++)
	{
		double s1 = 0;
		for (int j = 0; j < n; j++)
			s1 += fabs(a[i][j]);
		if (s1 - norm1 > eps)
			norm1 = s1;
	}

	cout << "a = \n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << a[i][j] << ' ';
		cout << endl;
	}
	cout << "b = \n";
	for (int i = 0; i < n; i++)
		cout << b[i] << ' ';
	cout << endl;
	cout << "y = \n";
	for (int i = 0; i < n; i++)
		cout << y[i] << ' ';
	cout << endl;

	// 2
	for (int i = 0; i < n; i++)
	{
			for (int j = 0; j < n; j++) 
				if (j != i)
				{
					double temp = a[j][i] / a[i][i];
					for (int k = 0; k < n; k++)
					{
						a[j][k] -= a[i][k] * temp;
						e[j][k] -= e[i][k] * temp;
					}
				}
	}

	for (int i = 0; i < n; i++)
	{
		double temp = a[i][i];
		a[i][i] /= temp;
		for (int j = 0; j < n; j++)
		{
			e[i][j] /= temp;
			if (i != j)
				a[i][j] /= temp;
		}
	}

	/*cout << "a_reverse = \n";
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << e[i][j] << ' ';
		cout << endl;
	}
	cout << endl;*/

	double norm2 = 0;
	for (int i = 0; i < n; i++)
	{
		double s2 = 0;
		for (int j = 0; j < n; j++)
			s2 += fabs(e[i][j]);
		if (s2 - norm2 > eps)
			norm2 = s2;
	}
	double number_of_conditioning = 0;
	number_of_conditioning = norm1 * norm2;
	cout << "Number of conditioning : " << fixed << setprecision(13) << number_of_conditioning << "\n";

	//3
	a = a_start;

	vector<double> ans_gauss(n, 0);
	vector<int> old(n, -1);
	for (int i = 0; i < n; i++)
		old[i] = i;
	for (int i = 0; i < n; i++)
	{
		int mn_pos_i = i, mn_pos_j = i;
		double mn = a[i][i];
		for (int ii = i; ii < n; ii++)
			for (int j = i; j < n; j++)
				if (fabs(a[ii][j]) > fabs(mn))
				{
					mn = a[ii][j];
					mn_pos_i = ii;
					mn_pos_j = j;
				}
		if (fabs(mn) > eps)
		{
			for (int j = i; j < n; j++)
			{
				swap(a[i][j], a[mn_pos_i][j]);
				swap(b[i], b[mn_pos_i]);
			}
			for (int j = i; j < n; j++)
			{
				swap(a[j][i], a[j][mn_pos_j]);
			}

			swap(old[i], old[mn_pos_j]);

		for (int j = i + 1; j < n; j++) 
			{
				double temp = a[j][i] / a[i][i];
				for (int k = 0; k < n; k++)
				{
					a[j][k] -= a[i][k] * temp;					
				}
				b[j] -= b[i] * temp;
			}
		}
	}

	for (int i = 0; i < n; i++) //обратный ход
	{
		for (int j = 0; j < i; j++)
		{
			if (j != i)
			{
				double temp = a[j][i] / a[i][i];
				a[j][i] -= a[i][i] * temp;
				b[j] -= b[i] * temp;
			}
			
		}
	}

	for (int i = 0; i < n; i++) 
	{
		ans_gauss[i] = b[i] / a[i][i];
	}
	

	cout << "Answer Gauss : ";
	for (int i = 0; i < n; i++)
		cout << ans_gauss[old[i]] << ' ';
	cout << "\n";

	system("pause");
	return 0;
}