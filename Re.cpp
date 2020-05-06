//Real part

#include "pch.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

double koeff = 0.6;

double X = 1, Y = 1;
int Nx = (int)(koeff * 100), Ny = (int)(koeff*100);

double _x = X / Nx, _y = Y / Ny;
double _x_sq = _x * _x, _y_sq = _y * _y;
double a_sq = _x_sq * _y_sq / (_x_sq + _y_sq);

double c = 1;
double w = 50;

double diel_const(int n)
{
	int h = 40;
	int j = (int)(n / Nx);
	int i = n - j * Nx;
	if ((j > koeff*30 && j < koeff*80) && (i > koeff*(50 - h) && i < koeff*(50 + h)))
		return 1.15;
	return 1;
}

int main()
{
	vector<double> layer_1(Nx*Ny), layer_2(Nx*Ny);

	double E0 = 1;

	for (int i = 0; i < Nx; i++)
	{
		int n = i;
		layer_1[n] = E0;
	}
	for (int j = 1; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			int n = i + Nx * j;
			layer_1[n] = 0;
		}
	}

	layer_2 = layer_1; //Можно ли так сделать?

	int T = 100;

	vector<double> layer(Nx*Ny);

	for (int k = 0; k < T; k++)
	{
		//Заполнение граничного слоя
		layer[0] = cos(w*k*sqrt(a_sq) / c);
		layer[Nx*(Ny - 1)] = layer_2[1 + Nx * (Ny - 2)];
		layer[Nx - 1 + Nx * (Ny - 1)] = layer_2[Nx - 2 + Nx * (Ny - 2)];
		layer[Nx - 1] = cos(w*k*sqrt(a_sq) / c);

		for (int i = 1; i < Nx - 1; i++)
		{
			layer[i] = cos(w*k*sqrt(a_sq) / c);
			layer[i + Nx * (Ny - 1)] = layer_2[i + Nx * (Ny - 2)];
		}
		for (int j = 1; j < Ny - 1; j++)
		{
			layer[Nx*j] = layer_2[1 + Nx * j];
			layer[Nx - 1 + Nx * j] = layer_2[Nx - 2 + Nx * j];
		}

		// Заполнение уровня

		for (int j = 1; j < Ny - 1; j++)
		{
			for (int i = 1; i < Nx - 1; i++)
			{
				layer[i + Nx * j] = a_sq / (_x_sq*diel_const(i + Nx * j)) * (layer_2[i + 1 + Nx * j] + layer_2[i - 1 + Nx * j]) + a_sq / (_y_sq*diel_const(i + Nx * j)) * (layer_2[i + Nx * (j + 1)] + layer_2[i + Nx * (j - 1)]) - layer_1[i + Nx * j];
			}
		}

		layer_1 = layer_2;
		layer_2 = layer;
	}

	std::ofstream out;
	out.open("Output_E_sample.txt");

	for (int j = 0; j < Ny; j++)
	{
		for (int i = 0; i < Nx; i++)
		{
			out << i << " " << j << " " << layer[i + Nx * j] << endl;
		}
	}
}
