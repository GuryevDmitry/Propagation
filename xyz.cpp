#include "pch.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

double X = 2, Y = 2, Z = 2;
int Nx = 20, Ny = 20, Nz = 20;

double _x = X / Nx, _y = Y / Ny, _z = Z / Nz;
double _x_sq = _x * _x, _y_sq = _y * _y, _z_sq = _z*_z;
double a_sq = 1 / (1 / _x_sq + 1 / _y_sq + 1 / _z_sq);

double c = 1;
double w = 10;

int main()
{
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

	int T = 100;

	vector<double> l(Nx*Ny*Nz);

	for (int k = 0; k < T; k++)
	{
		//Заполнение граничного слоя
		double boundary = E0 * cos(w*k*sqrt(a_sq) / c);

		l[0] = boundary;
		l[Nx-1] = boundary;
		l[Nx*Ny*(Nz-1)] = boundary;
		l[Nx-1 + Nx * Ny*(Nz - 1)] = boundary;

		l[Nx*(Ny - 1)] = l_2[1 + Nx * (Ny - 2) + Nx * Ny];
		l[Nx - 1 + Nx * (Ny - 1)] = l_2[Nx - 2 + Nx * (Ny - 2) + Nx * Ny];
		l[Nx*(Ny - 1) + Nx * Ny*(Nz - 1)] = l_2[1 + Nx * (Ny - 2) + Nx * Ny*(Nz - 2)];
		l[Nx - 1 + Nx * (Ny - 1) + Nx * Ny*(Nz - 1)] = l_2[Nx - 2 + Nx * (Ny - 2) + Nx * Ny*(Nz - 2)];

		for (int i = 1; i < Nx-1; i++)
		{
			for (int s = 1; s < Nz-1; s++)
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
					l[(i)+Nx * (j)+Nx * Ny*(s)] = (a_sq / _x_sq) * (l_2[(i + 1) + Nx * (j)+Nx * Ny*(s)] + l_2[(i - 1) + Nx * (j)+Nx * Ny*(s)]) + \
						(a_sq / _y_sq) * (l_2[(i)+Nx * (j + 1) + Nx * Ny*(s)] + l_2[(i)+Nx * (j - 1) + Nx * Ny*(s)]) + \
						(a_sq / _z_sq) * (l_2[(i)+Nx * (j)+Nx * Ny*(s + 1)] + l_2[(i)+Nx * (j)+Nx * Ny*(s - 1)]) - \
						l_1[(i)+Nx * (j)+Nx * Ny*(s)];

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
			out << i << " " << j << " " << l[i + Nx * j+Nx*Ny*(int)(Nz/2)] << endl;
		}
	}
}
