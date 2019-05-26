#define _USE_MATH_DEFINES
#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include<iomanip>
#include<cstdio>
#include<cmath>
#include <random>

using namespace std;

const int n = 10, N = 3;
const double eps = 0.7, eps_roots = 0.0001, eps_roots_n = 0.00000001, eps_roots_n_2 = 1e-25, h = 0.00001, pi = M_PI;

vector<double> f_multiplication(vector<vector<double> > a, vector<double> b)
{
	vector<double> ans(n);
	for (int i = 0; i < n; i++)
	{
		double val = 0;
		for (int j = 0; j < n; j++)
			val += a[i][j] * b[j];
		ans[i] = val;
	}
	return ans;
}

double func(double x)
{
	return (((pow(x, 9) + pi)*cos(log(x * x + 1.0))) / exp(x * x)) - (x / 2018.0);
}

double derivative_func(double x)
{
	return x*exp(-x*x)*cos(log(x*x+1))*(9.0*pow(x,7)-2.0*(pow(x,9)+pi))-2.0*(pow(x,9)+pi)*sin(log(x*x+1.0))*(-1.0/(x*x+1.0))+(-1.0/2018.0);
}

int main()
{
	setlocale(LC_ALL, "Russian");
	cin.tie(0);
	ios_base::sync_with_stdio(0);
	vector<pair<double, double> > norma_step(N);
	double time_average_eigenvalue = 0;

	// 1
	for (int glob_i = 0; glob_i < N; glob_i++)
	{
		srand(time(0));
		vector<vector<double> > a(n), a_start(n);
		vector<double> y, sum(n, 0), b(n, 0);
		for (int i = 0; i < n; i++)
		{
			a[i].resize(n, 0);
			a_start[i].resize(n, 0);
		}
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				a[i][j] = -2.0 + ((double)rand() / (double)RAND_MAX) * 4.0;
		a_start = a;

		/*cout << "a = \n";
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
				cout << a[i][j] << ' ';
			cout << endl;
		}
		cout << endl;*/

		//2
		vector<double> u_prev(n, 0), lambda1_vector(n);
		// lambda1
		u_prev[0] = 1.0;
		double prev_lambda = 0, lambda1 = 0;
		int k = 0;

		clock_t start, end;
		start = clock();

		while (1)
		{
			k++;
			vector<double> y_new(n), u_new(n);
			y_new = f_multiplication(a, u_prev);
			lambda1 = 0;
			for (int i = 0; i < n; i++)
			{
				if (fabs(y_new[i]) > fabs(lambda1))
					lambda1 = y_new[i];
			}
			for (int i = 0; i < n; i++)
				u_new[i] = y_new[i] / lambda1;
			if (fabs(prev_lambda - lambda1) < eps) break;
			if (k == 100000)
			{
				break;
			}
			u_prev = u_new;
		}
		lambda1_vector = u_prev;

		end = clock();
		time_average_eigenvalue += ((double)end - start) / ((double)CLOCKS_PER_SEC);

		vector<double> temp1;
		temp1 = f_multiplication(a, lambda1_vector);
		double norma;
		for (int i = 0; i < n; i++)
		{
			temp1[i] -= lambda1 * lambda1_vector[i];
			if (i == 0)
				norma = temp1[i];
			else
				norma = max(norma, temp1[i]);
		}
		norma_step[glob_i].first = norma;

		/*cout << fixed << setprecision(6) << lambda1 << "\n";
		for (int i = 0; i < n; i++)
			cout << lambda1_vector[i] << ' ';
		cout << "\n";*/

		// lambda2
		k = 0;
		vector<double> lambda2_vector(n);
		for (int i = 0; i < n; i++)
			a[i][i] -= lambda1,
			u_prev[i] = 0;
		u_prev[0] = 1;
		prev_lambda = 0;
		double lambda2 = 0;
		while (1)
		{
			k++;
			vector<double> y_new(n), u_new(n);
			y_new = f_multiplication(a, u_prev);
			lambda2 = 0;
			for (int i = 0; i < n; i++)
			{
				if (fabs(y_new[i]) > fabs(lambda2))
					lambda2 = y_new[i];
			}
			for (int i = 0; i < n; i++)
				u_new[i] = y_new[i] / lambda2;
			if (fabs(prev_lambda - lambda2) < eps) break;
			if (k == 100000)
			{
				break;
			}
			u_prev = u_new;
		}
		lambda2 += lambda1;
		lambda2_vector = u_prev;

		temp1.clear();
		temp1 = f_multiplication(a, lambda2_vector);
		for (int i = 0; i < n; i++)
		{
			temp1[i] -= lambda2 * lambda2_vector[i];
			if (i == 0)
				norma = temp1[i];
			else
				norma = max(norma, temp1[i]);
		}
		norma_step[glob_i].second = norma;

		/*cout << fixed << setprecision(6) << lambda2 << "\n";
		for (int i = 0; i < n; i++)
			cout << lambda2_vector[i] << ' ';
		cout << "\n";*/

		//3
		//a = a_start;
	}

		// 4, 5, 6
		vector<double> roots, ans_roots;
		vector<pair<double, double> > segments;
		segments.push_back({ -2.0, -1.5 });
		segments.push_back({ -1.4, -0.9 });
		segments.push_back({ 1.5, 2.0 });
		vector<int> cnt;

		for (int i = 0; i < segments.size(); i++)
		{
			double a = segments[i].first, b = segments[i].second, c;
			int count = 0;
			while (1)
			{
				c = (a + b) / 2.0;
				count++;
				if (fabs(a - b) <= eps_roots)
					break;
				if (func(a) * func(c) < 0)
					b = c;
				else
					if (func(b) * func(c) < 0)
						a = c;
			}
			segments[i].first = a, segments[i].second = b;
			ans_roots.push_back(c);
			cnt.push_back(count);
		}

		/*cout << "_____cnt__ segment____ans_roots_\n";
		for (int i = 0; i < cnt.size(); i++)
		{
			cout << fixed << setprecision(13) << cnt[i] << ' ' << segments[i].first << ' ' << segments[i].second << ' ' << ans_roots[i] << "\n";
		}*/

		// 4, 5, 7
		vector<double> x_start;
		x_start.push_back(-2.5);
		x_start.push_back(-1.0);
		x_start.push_back(2.0);
		vector<int> cnt_dn;
		vector<double> ans_roots_dn;
		for (int i = 0; i < segments.size(); i++)
		{
			double c;
			int count = 0;
			double x_prev = x_start[i];
			while (1)
			{
				c = x_prev - (h * func(x_prev) / (func(x_prev + h) - func(x_prev)));
				count++;
				if (fabs(c - x_prev) <= eps_roots_n)
					break;
				x_prev = c;
			}
			ans_roots_dn.push_back(c);
			cnt_dn.push_back(count);
		}

		/*cout << "_____cnt______ans_roots_dn\n";
		for (int i = 0; i < cnt_dn.size(); i++)
		{
			cout << fixed << setprecision(13) << cnt_dn[i] << ' ' << ans_roots_dn[i] << "\n";
		}*/

		//4, 5, 8
		x_start.clear();
		x_start.push_back(-2.5); 
		x_start.push_back(-1.0); 
		x_start.push_back(2.0);
		vector<int> cnt_n;
		vector<double> ans_roots_n;
		for (int i = 0; i < segments.size(); i++)
		{
			double c;
			int count = 0;
			double x_prev = x_start[i];
			while (1)
			{
				c = x_prev - func(x_prev) / derivative_func(x_prev);
				count++;
				if (fabs(c - x_prev) <= eps_roots)
					break;
				x_prev = c;
			}
			ans_roots_n.push_back(c);
			cnt_n.push_back(count);
		}

		/*cout << "_____cnt______ans_roots_n\n";
		for (int i = 0; i < cnt_n.size(); i++)
		{
			cout << fixed << setprecision(13) << cnt_n[i] << ' ' << ans_roots_n[i] << "\n";
		}*/
	

	cout << "Норма разности A*x1 - lambda*x1 :\n";
	for (int i = 0; i < N; i++)
		cout << fixed << setprecision(13) << norma_step[i].first << ' ' << norma_step[i].second << "\n";
	cout << "Среднее время нахождения собственного значения и вектора степенным методом :\n";
	cout << time_average_eigenvalue / N << "\n";
	cout << "Конечные отрезки отделённости корней методом бисекции :\n";
	for (int i = 0; i < segments.size(); i++)
		cout << fixed << setprecision(13) << segments[i].first << ' ' << segments[i].second << "\n";
	cout << "Количество шагов метода бисекции для каждого из корней :\n";
	for (int i = 0; i < segments.size(); i++)
		cout << cnt[i] << "\n";
	cout << "Корни уравнения дискретным методом Ньютона :\n";
	for (int i = 0; i < segments.size(); i++)
		cout << fixed << setprecision(20) << ans_roots_dn[i] << "\n";
	cout << "Количество шагов дискретного варианта метода Ньютона для каждого из корней :\n";
	for (int i = 0; i < segments.size(); i++)
		cout << cnt_dn[i] << "\n";
	cout << "Корни уравнения методом Ньютона :\n";
	for (int i = 0; i < segments.size(); i++)
		cout << fixed << setprecision(20) << ans_roots_n[i] << "\n";
	cout << "Количество шагов метода Ньютона для каждого из корней :\n";
	for (int i = 0; i < segments.size(); i++)
		cout << cnt_n[i] << "\n";
	system("pause");
	return 0;
}