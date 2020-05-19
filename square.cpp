#include "pch.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

double X = 10, Y = 4, Z = 10;
int Nx = 60, Ny = 30, Nz = 60;

double _x = X / Nx, _y = Y / Ny, _z = Z / Nz;
double _x_sq = _x * _x, _y_sq = _y * _y, _z_sq = _z * _z;
double _a = 0.5 / (1 / _x + 1 / _y + 1 / _z);
double a_sq = _a*_a;

double c = 1;
double w = 15;

bool in(int i, int j, int s)
{
	return !((j > 6) && (j < 9) && !((i >= 27) && (i <= 33)) && !((s >= 27) && (s <= 33)));
//	return 1;
}

int main()
{
	cout << _a << endl;
	cout << _a / c << endl;

	vector<double> l_1(Nx*Ny*Nz), l_2(Nx*Ny*Nz);

	double E0 = 1;

	for (int i = 0; i < Nx; i++)
	{
		for (int s = 0; s < Nz; s++)
		{
			l_1[i + Nx * Ny*s] = E0;
		}
	}
	for (int j = 1; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			for (int s = 0; s < Nz; s++)
			{
				int n = i + Nx * j + Nx * Ny*s;
				l_1[n] = 0;
			}
		}
	}

	l_2 = l_1;

	int T = 500;

	cout << T * _a / c << endl;

	vector<double> l(Nx*Ny*Nz);

	for (int k = 0; k < T; k++)
	{
		//Заполнение граничного слоя
		double boundary = E0 * cos(w*k*sqrt(a_sq) / c);
		cout << k << " " << boundary << endl;

		l[0] = boundary;
		l[Nx - 1] = boundary;
		l[Nx*Ny*(Nz - 1)] = boundary;
		l[Nx - 1 + Nx * Ny*(Nz - 1)] = boundary;

		l[Nx*(Ny - 1)] = l_2[1 + Nx * (Ny - 2) + Nx * Ny];
		l[Nx - 1 + Nx * (Ny - 1)] = l_2[Nx - 2 + Nx * (Ny - 2) + Nx * Ny];
		l[Nx*(Ny - 1) + Nx * Ny*(Nz - 1)] = l_2[1 + Nx * (Ny - 2) + Nx * Ny*(Nz - 2)];
		l[Nx - 1 + Nx * (Ny - 1) + Nx * Ny*(Nz - 1)] = l_2[Nx - 2 + Nx * (Ny - 2) + Nx * Ny*(Nz - 2)];

		for (int i = 1; i < Nx - 1; i++)
		{
			l[i] = boundary;
			l[i + Nx * Ny*(Nz - 1)] = boundary;
			l[i + Nx * (Ny - 1)] = l_2[i + Nx * (Ny - 2) + Nx * Ny];
			l[i + Nx * (Ny - 1) + Nx * Ny*(Nz - 1)] = l_2[i + Nx * (Ny - 2) + Nx * Ny*(Nz - 2)];
		}
		for (int s = 1; s < Nz - 1; s++)
		{
			l[Nx*Ny*s] = boundary;
			l[Nx - 1 + Nx * Ny*s] = boundary;;
			l[Nx * (Ny - 1) + Nx * Ny*s] = l_2[1 + Nx * (Ny - 2) + Nx * Ny*s];
			l[Nx - 1 + Nx * (Ny - 1) + Nx * Ny*s] = l_2[Nx - 2 + Nx * (Ny - 2) + Nx * Ny*s];
		}
		for (int j = 1; j < Ny - 1; j++)
		{
			l[Nx*j] = l_2[1 + Nx * j + Nx * Ny];
			l[Nx - 1 + Nx * j] = l_2[Nx - 2 + Nx * j + Nx * Ny];
			l[Nx*j + Nx * Ny*(Nz - 1)] = l_2[1 + Nx * j + Nx * Ny*(Nz - 2)];
			l[Nx - 1 + Nx * j + Nx * Ny*(Nz - 1)] = l_2[Nx - 2 + Nx * j + Nx * Ny*(Nz - 2)];
		}
		
		for (int i = 1; i < Nx - 1; i++)
		{
			for (int s = 1; s < Nz - 1; s++)
			{
				l[i + Nx * Ny*s] = boundary;
				l[i + Nx * (Ny - 1) + Nx * Ny*s] = l_2[i + Nx * (Ny - 2) + Nx * Ny*s];
			}
		}
		for (int i = 1; i < Nx - 1; i++)
		{
			for (int j = 1; j < Ny - 1; j++)
			{
				l[i + Nx * j] = l_2[i + Nx * j + Nx * Ny];
				l[i + Nx * j + Nx * Ny*(Nz - 1)] = l_2[i + Nx * j + Nx * Ny*(Nz - 2)];
			}
		}
		for (int s = 1; s < Nz - 1; s++)
		{
			for (int j = 1; j < Ny - 1; j++)
			{
				l[Nx*j + Nx * Ny*s] = l_2[1 + Nx * j + Nx * Ny*s];
				l[Nx - 1 + Nx * j + Nx * Ny*s] = l_2[Nx - 2 + Nx * j + Nx * Ny*s];
			}
		}

		// Заполнение уровня

		for (int j = 1; j < Ny - 1; j++)
		{
			for (int i = 1; i < Nx - 1; i++)
			{
				for (int s = 1; s < Nz - 1; s++)
				{
					double E_ijsk = l_2[(i)+Nx * (j)+Nx * Ny*(s)];
					if (in(i,j,s) == 1)
					{
						l[(i)+Nx * (j)+Nx * Ny*(s)] = (a_sq / _x_sq) * (l_2[(i + 1) + Nx * (j)+Nx * Ny*(s)] + l_2[(i - 1) + Nx * (j)+Nx * Ny*(s)] - 2 * E_ijsk) + \
							(a_sq / _y_sq) * (l_2[(i)+Nx * (j + 1) + Nx * Ny*(s)] + l_2[(i)+Nx * (j - 1) + Nx * Ny*(s)] - 2 * E_ijsk) + \
							(a_sq / _z_sq) * (l_2[(i)+Nx * (j)+Nx * Ny*(s + 1)] + l_2[(i)+Nx * (j)+Nx * Ny*(s - 1)] - 2 * E_ijsk) - \
							l_1[(i)+Nx * (j)+Nx * Ny*(s)] + 2 * E_ijsk;
					}
					else
					{
						l[(i)+Nx * (j)+Nx * Ny*(s)] = 0;
					}
					
				}
			}
		}

		l_1 = l_2;
		l_2 = l;
	}

	std::ofstream out;
	out.open("Output_E_sample.txt");

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			out << i << " " << j << " " << l[i + Nx * j + Nx * Ny* Nz / 2] << endl;
		}
	}
	
	std::ofstream out2;
	out2.open("on_screen.txt");

	for (int s = 0; s < Nz; s++)
	{
		for (int i = 0; i < Nx; i++)
		{
			out2 << i << " " << s << " " << l[i + Nx * (Ny-10) + Nx * Ny*(s)] << endl;
		}
	}
	/*
	std::ofstream out3;
	out3.open("geom.txt");

	for (int s = 0; s < Nz; s++)
	{
		for (int i = 0; i < Nx; i++)
		{
			out3 << i << " " << s << " " << l[i + Nx * (17) + Nx * Ny*(s)] << endl;
		}
	}
	*/
}
