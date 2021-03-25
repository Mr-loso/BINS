#include <math.h>
#include <omp.h>
#include <iostream>
using namespace std;

//Число Пи
#define M_PI	3.14159265358979323846

//перевод из градусов в радианы
double gr = M_PI / 180;

//перевод из радиан в градусы
double rg = 180 / M_PI;

//угловая скорость вращения Земли
double wzem = 0.0042 * gr;

//вектор ускорения силы тяжести
double g = 9.81;

//радиус Земли
double R = 6500 * 1e3;

//начальные значения для gamma(крен)
double gamma0 = 0 * gr, A_gamma = 3 * gr, w_gamma = 0.1 * gr;
//начальные значения для psi(курс)
double psi0 = 45 * gr, A_psi = 1 * gr, w_psi = 0.2 * gr;
//начальные значения для teta(дифферент)
double teta0 = 0 * gr, A_teta = 2 * gr, w_teta = 0.4 * gr;

//начальные значения для формирования вектора кажущегося ускорения
double A_ax = 0.5, w_ax = 0.6;
double A_ay = 1, w_ay = 0.2;
double A_az = 1.5, w_az = 0.4;

//Указатель на функцию
void(*fff)(double*, double*, double);
//функция вычисления правых частей системы ОДУ
void fprav(double*, double*, double);
//Метод Эйлера с контролем точности
void ejler(double*, int, double, double, double, double*, void(*func)(double*, double*, double));
//Функция вычисления wxiz wyiz wziz akxiz akyiz akziz для формирования заданного движения
void matrix_calc(double* wx, double* wy, double* wz, double* akx, double* aky, double* akz, double t);

//Указатель на вектор
double* y, * x, * z, * z1, * y11;

//количество потоков
#define threads 4

int main()
{
	omp_set_num_threads(threads);
	int N;
	double t, t0, tK, eps, step;
	double* h, h0;

	//Исходные данные
	t0 = 0;		//начальный момент
	tK = 0.01;		//конечный момент 
	step = h0 = 1e-6;	//начальный шаг интегрирования
	eps = 1e-15;		//требуемая точность

	N = 17; //размерность вектора состояния

//выделение памяти 
	y = (double*)malloc(sizeof(double) * N);
	x = (double*)malloc(sizeof(double) * N);
	z = (double*)malloc(sizeof(double) * N);
	z1 = (double*)malloc(sizeof(double) * N);
	y11 = (double*)malloc(sizeof(double) * N);

	//начальные значения вектора состояния
	y[0] = cos(psi0) * cos(teta0);
	y[1] = -sin(psi0) * sin(gamma0) - cos(psi0) * sin(teta0) * cos(gamma0);
	y[2] = -sin(psi0) * cos(gamma0) + cos(psi0) * sin(teta0) * sin(gamma0);
	y[3] = sin(teta0);
	y[4] = cos(teta0) * cos(gamma0);
	y[5] = -cos(teta0) * sin(gamma0);
	y[6] = sin(psi0) * cos(teta0);
	y[7] = cos(psi0) * sin(gamma0) - sin(psi0) * cos(gamma0) * sin(teta0);
	y[8] = cos(psi0) * cos(gamma0) + sin(psi0) * sin(gamma0) * sin(teta0);

	//начальные условия для скоростей vksi veta vdzeta
	double vx_0 = 100, vy_0 = 0, vz_0 = 0;
	y[9] = y[0] * vx_0 + y[1] * vy_0 + y[2] * vz_0;
	y[10] = y[3] * vx_0 + y[4] * vy_0 + y[5] * vz_0;
	y[11] = y[6] * vx_0 + y[7] * vy_0 + y[8] * vz_0;

	//начальные условия для координат
	y[12] = 0; y[13] = 0; y[14] = 0;

	//начальные условия для широты и долготы
	double fi0 = 50 * gr;
	double lyambda0 = 30 * gr;
	y[15] = fi0;
	y[16] = lyambda0;

	//Формирование указателей
	h = &h0;
	fff = fprav;

	//Численное решение системы ОДУ
	for (t = t0; t < tK; t += *h)
	{
		ejler(y, N, eps, t, step, h, fff);
	}

	double lyambda1 = atan((-1) * ((cos(z[15] * sin(z[16])) / (sin(z[15])))));
	double fi1 = asin((cos(z[15])) * (cos(z[16])));
	double per = atan((-1) * ((sin(z[16])) / (sin(z[15]) * cos(z[16]))));
	double psi1 = (atan(x[6] / x[0])) + per;

	y11[0] = cos(psi1) * cos(asin(x[3]));
	y11[1] = -sin(psi1) * sin(-atan(x[5] / x[4])) - cos(psi1) * sin(asin(x[3])) * cos(-atan(x[5] / x[4]));
	y11[2] = -sin(psi1) * cos(-atan(x[5] / x[4])) + cos(psi1) * sin(asin(x[3])) * sin(-atan(x[5] / x[4]));
	y11[3] = sin(asin(x[3]));
	y11[4] = cos(asin(x[3])) * cos(-atan(x[5] / x[4]));
	y11[5] = -cos(asin(x[3])) * sin(-atan(x[5] / x[4]));
	y11[6] = sin(psi1) * cos(asin(x[3]));
	y11[7] = cos(psi1) * sin(-atan(x[5] / x[4])) - sin(psi1) * cos(-atan(x[5] / x[4])) * sin(asin(x[3]));
	y11[8] = cos(psi1) * cos(-atan(x[5] / x[4])) + sin(psi1) * sin(-atan(x[5] / x[4])) * sin(asin(x[3]));

	free(x); free(z); free(z1); free(y); 	//освобождение памяти
	return 0;


}//метод Эйлера с контролем точности
void ejler(double* y, int n, double eps, double t, double step, double* h, void(*func)(double*, double*, double))
{
	double  delt, q;
	*h = step; //начальный шаг

	for (int i = 0; i < n; i++) x[i] = y[i];
	do {
		for (int i = 0; i < n; i++) y[i] = x[i];
		(*func)(z, y, t); //вычисление z от y(n)
		for (int i = 0; i < n; i++) y[i] += *h * z[i]; //вычисление y(n+1)
		(*func)(z1, y, t); //вычисление z от y(n+1)
		//вычисление нормы
		delt = 0;
		for (int i = 0; i < n; i++) delt += (z[i] - z1[i]) * (z[i] - z1[i]);
		delt = 0.5 * (*h) * sqrt(delt);
		//вычисление q
		q = sqrt(eps / delt);
		//конроль q и изменение шага
		if (q < 1) *h *= q / 1.1;
	} while (q < 1);
}




//функция вычисления правых частей системы ОДУ
void fprav(double* z, double* x, double t)
{
	double wx = 0, wy = 0, wz = 0, akx = 0, aky = 0, akz = 0;
	double* wx0, * wy0, * wz0, * akx0, * aky0, * akz0;
	wx0 = &wx; wy0 = &wy; wz0 = &wz; akx0 = &akx; aky0 = &aky; akz0 = &akz;
	//формирования заданного движения
	matrix_calc(wx0, wy0, wz0, akx0, aky0, akz0, t);
	double gamma = 0, Kq = 0, teta = 0;
	gamma = -atan(x[5] / x[4]);
	Kq = atan(x[6] / x[0]);
	teta = asin(x[3]);


	//вычисление матрицы B
#pragma omp parallel sections
	{
#pragma omp section

		{
			z[0] = (cos(teta) * sin(Kq)) * (wy + (x[10] / R) * (-cos(teta) * sin(gamma)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (sin(teta)) - cos(gamma) * cos(teta) * sin(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem))
				+ (sin(teta)) * (wz + (x[10] / R) * (cos(gamma) * cos(Kq) + sin(gamma) * sin(Kq) * sin(teta)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (cos(teta) * sin(Kq)) + sin(x[16]) * (cos(Kq) * sin(gamma) - cos(gamma) * sin(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem));
		}
#pragma omp section
		{
			z[1] = -(cos(teta) * sin(Kq)) * (wx - sin(x[16]) * (sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem) + (x[10] / R) * cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta) - cos(x[16]) * (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem))
				- (cos(Kq) * cos(teta)) * (wz + (x[10] / R) * (cos(gamma) * cos(Kq) + sin(gamma) * sin(Kq) * sin(teta)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (cos(teta) * sin(Kq)) + sin(x[16]) * (cos(Kq) * sin(gamma) - cos(gamma) * sin(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem));
		}
#pragma omp section
		{
			z[2] = -(sin(teta)) * (wx - sin(x[16]) * (sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem) + (x[10] / R) * cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta) - cos(x[16]) * (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem))
				+ (cos(Kq) * cos(teta)) * (wy + (x[10] / R) * (-cos(teta) * sin(gamma)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (sin(teta)) - cos(gamma) * cos(teta) * sin(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem));
		}
#pragma omp section
		{
			z[3] = -(cos(Kq) * sin(gamma) - cos(gamma) * sin(Kq) * sin(teta)) * (wy + (x[10] / R) * (-cos(teta) * sin(gamma)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (sin(teta)) - cos(gamma) * cos(teta) * sin(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem))
				+ (cos(teta) * cos(gamma)) * (wz + (x[10] / R) * (cos(gamma) * cos(Kq) + sin(gamma) * sin(Kq) * sin(teta)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (cos(teta) * sin(Kq)) + sin(x[16]) * (cos(Kq) * sin(gamma) - cos(gamma) * sin(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem));
		}
#pragma omp section
		{
			z[4] = -(sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * (wz + (x[10] / R) * (cos(gamma) * cos(Kq) + sin(gamma) * sin(Kq) * sin(teta)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (cos(teta) * sin(Kq)) + sin(x[16]) * (cos(Kq) * sin(gamma) - cos(gamma) * sin(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem))
				+ (cos(Kq) * sin(gamma) - cos(gamma) * sin(Kq) * sin(teta)) * (wx - sin(x[16]) * (sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem) + (x[10] / R) * cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta) - cos(x[16]) * (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem));
		}
#pragma omp section
		{
			z[5] = (sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * (wy + (x[10] / R) * (-cos(teta) * sin(gamma)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (sin(teta)) - cos(gamma) * cos(teta) * sin(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem))
				- (cos(teta) * cos(gamma)) * (wx - sin(x[16]) * (sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem) + (x[10] / R) * cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta) - cos(x[16]) * (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem));
		}
#pragma omp section
		{
			z[6] = (-cos(teta) * sin(gamma)) * (wz + (x[10] / R) * (cos(gamma) * cos(Kq) + sin(gamma) * sin(Kq) * sin(teta)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (cos(teta) * sin(Kq)) + sin(x[16]) * (cos(Kq) * sin(gamma) - cos(gamma) * sin(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem)) - (cos(gamma) * cos(Kq) + sin(gamma) * sin(Kq) * sin(teta)) * (wy + (x[10] / R) * (-cos(teta) * sin(gamma)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (sin(teta)) - cos(gamma) * cos(teta) * sin(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem));
		}
#pragma omp section
		{
			z[7] = (cos(gamma) * cos(Kq) + sin(gamma) * sin(Kq) * sin(teta)) * (wx - sin(x[16]) * (sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem) + (x[10] / R) * cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta) - cos(x[16]) * (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem))
				- (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * (wz + (x[10] / R) * (cos(gamma) * cos(Kq) + sin(gamma) * sin(Kq) * sin(teta)) + cos(x[16]) * ((x[12] / (R * cos(x[16]))) + wzem) * (cos(teta) * sin(Kq)) + sin(x[16]) * (cos(Kq) * sin(gamma) - cos(gamma) * sin(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem));
		}
#pragma omp section
		{
			z[8] = (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * (wx - sin(x[16]) * (sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem) + (x[10] / R) * cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta) - cos(x[16]) * (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem))
				- (-cos(teta) * sin(gamma)) * (wx - sin(x[16]) * (sin(gamma) * sin(Kq) + cos(gamma) * cos(Kq) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem) + (x[10] / R) * cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta) - cos(x[16]) * (cos(gamma) * sin(Kq) - cos(Kq) * sin(gamma) * sin(teta)) * ((x[12] / (R * cos(x[16]))) + wzem));
		}

		//скорости vksi veta vdzeta
#pragma omp section
		{
			z[9] = x[10] * (wx * x[6] + wy * x[7] + wz * x[8])
				+ akx * x[0] + aky * x[1] + akz * x[2]
				- x[11] * (wx * x[3] + wy * x[4] + wz * x[5]
					+ wzem * sin(x[15]));
		}
#pragma omp section
		{
			z[10] = g - x[9] * (wx * x[6] + wy * x[7] + wz * x[8])
				+ akx * x[3] + aky * x[4] + akz * x[5]
				+ x[11] * (wx * x[0] + wy * x[1] + wz * x[2]
					+ wzem * cos(x[15]));
		}
#pragma omp section
		{
			z[11] = akx * x[6] + aky * x[7] + akz * x[8]
				- x[10] * (wx * x[0] + wy * x[1] + wz * x[2]
					+ wzem * cos(x[15])) + x[9] * (wx * x[3]
						+ wy * x[4] + wz * x[5] + wzem * sin(x[15]));
		}
		//координаты ksi eta dzeta
#pragma omp section
		{
			z[12] = x[9];
		}
#pragma omp section
		{
			z[13] = x[10];
		}
#pragma omp section
		{
			z[14] = x[11];
		}
		//широта fi долгота lyambda
#pragma omp section
		{
			z[15] = x[9] / R;
		}
#pragma omp section
		{
			z[16] = x[11] / (R * cos(x[15]));
		}
	}
}

//Функция вычисления wxiz wyiz wziz akxiz akyiz akziz для формирования заданного движения
void matrix_calc(double* wx, double* wy, double* wz, double* akx, double* aky, double* akz, double t)
{
	//wxyziz
	double wxiz, wyiz, wziz;
	wxiz = A_gamma * w_gamma * cos(w_gamma * t) + A_psi * w_psi * cos(w_psi * t) * sin(teta0 + A_teta * sin(w_teta * t));

	wyiz = A_teta * w_teta * cos(w_teta * t) * sin(gamma0 + A_gamma * sin(w_gamma * t)) + A_psi * w_psi * cos(w_psi * t)
		* cos(gamma0 + A_gamma * sin(w_gamma * t)) * cos(teta0 + A_teta * sin(w_teta * t));

	wziz = A_teta * w_teta * cos(w_teta * t) * cos(gamma0 + A_gamma * sin(w_gamma * t))
		- A_psi * w_psi * cos(w_psi * t)
		* sin(gamma0 + A_gamma * sin(w_gamma * t)) * cos(teta0 + A_teta * sin(w_teta * t));

	//akxyziz
	double akxiz = A_ax * sin(w_ax * t);
	double akyiz = A_ay * sin(w_ay * t);
	double akziz = A_az * sin(w_az * t);

	//пересчет wxyziz в wxyz
	double B[3][3] = { {x[0], x[1], x[2] }, {x[3], x[4], x[5]}, {x[6], x[7], x[8]} };
	double Btr[3][3] = { {x[0], x[3], x[6] }, {x[1], x[4], x[7]}, {x[2], x[5], x[8]} };
	double wzemM[] = { wzem * cos(x[15]), wzem * sin(x[15]), 0 };
	double wdev[] = { cos(x[15]) * (x[11] / (R * cos(x[15]))), sin(x[15]) * (x[11] / (R * cos(x[15]))), -(x[9] / R) };

	double buffer_1[3] = { B[0][0] * wxiz + B[0][1] * wyiz + B[0][2] * wziz,
	B[1][0] * wxiz + B[1][1] * wyiz + B[1][2] * wziz,
	B[2][0] * wxiz + B[2][1] * wyiz + B[2][2] * wziz };

	double buffer_2[3] = {
	buffer_1[0] + wzemM[0] + wdev[0],
	buffer_1[1] + wzemM[1] + wdev[1],
	buffer_1[2] + wzemM[2] + wdev[2] };

	double wxyz[3] = { Btr[0][0] * buffer_2[0] + Btr[0][1] * buffer_2[1] + Btr[0][2] * buffer_2[2],
		Btr[1][0] * buffer_2[0] + Btr[1][1] * buffer_2[1] + Btr[1][2] * buffer_2[2],
		Btr[2][0] * buffer_2[0] + Btr[2][1] * buffer_2[1] + Btr[2][2] * buffer_2[2] };

	*wx = wxyz[0];
	*wy = wxyz[1];
	*wz = wxyz[2];

	//пересчет akxyziz в akxyz
	double G[3] = { 0, -g, 0 };
	double A[3] =
	{ buffer_2[0] + wzemM[0] ,
		buffer_2[1] + wzemM[1] ,
		buffer_2[2] + wzemM[2] };

	double Ax = A[0];
	double Ay = A[1];
	double Az = A[2];
	double P[3][3] = { {0, -Az , Ay}, {Az, 0, -Ax}, {-Ay, Ax, 0} };
	double X[3] =
	{ P[0][0] * x[9] + P[0][1] * x[10] + P[0][2] * x[11] ,
		P[1][0] * x[9] + P[1][1] * x[10] + P[1][2] * x[11] ,
		P[2][0] * x[9] + P[2][1] * x[10] + P[2][2] * x[11] };

	double B_mul_akxyziz[3] = { B[0][0] * akxiz + B[0][1] * akyiz + B[0][2] * akziz,
	B[1][0] * akxiz + B[1][1] * akyiz + B[1][2] * akziz,
	B[2][0] * akxiz + B[2][1] * akyiz + B[2][2] * akziz };
	double newBakxyziz[3] = {
	G[0] + B_mul_akxyziz[0] + X[0],
	G[1] + B_mul_akxyziz[1] + X[1],
	G[2] + B_mul_akxyziz[2] + X[2] };

	double akxyz[3] = { Btr[0][0] * newBakxyziz[0] + Btr[0][1] * newBakxyziz[1] + Btr[0][2] * newBakxyziz[2],
	Btr[1][0] * newBakxyziz[0] + Btr[1][1] * newBakxyziz[1] + Btr[1][2] * newBakxyziz[2],
	Btr[2][0] * newBakxyziz[0] + Btr[2][1] * newBakxyziz[1] + Btr[2][2] * newBakxyziz[2] };

	*akx = akxyz[0];
	*aky = akxyz[1];
	*akz = akxyz[2];
}