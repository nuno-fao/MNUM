#include "funcs.h"


double m(double x) {
	return pow(x, 3) -10*sin(x)+2.8;
}

void bissec() {
	double a = 1.5, b = 4.2, x;
	int counter = 0;

	for (int i = 0; i < 3; i++) {
		x = (a + b) / 2;
		cout << "a: " << a << endl;
		cout << "b: " << b << endl;
		cout << "x: " << x << endl << endl;
		if (m(a) * m(x) < 0) {
			b = x;
		}
		else {
			a = x;
		}
		
	}

}

double j(double x) {
	return pow(x, 7) + 0.5 * x - 0.5;
}

void corda() {
	double a = 0, b = 0.8, x = 0.656044;
	int counter = 1;

	for (int i = 0; i < 3; i++) {
		if (j(a) * j(x) < 0) {
			b = x;
		}
		else {
			a = x;
		}
		x = (j(a) * b - j(b) * a) / (j(a) - j(b));
		cout << "a: " << a << endl;
		cout << "b: " << b << endl;
		cout << "x: " << x << endl << endl;
	}

}

double g(double x) {
	return log(5+x);
}
void picard_peano() {
	double x = 3;
	int counter = 0;
	double temp = g(x); 
	cout << "Iter: " << counter << endl;
	cout << "Valor: " << temp << endl;
	while (abs(temp - x) > pow(10,-5)) {
		x = temp;
		temp = g(x);
		counter++;
		cout << "Iter: " << counter << endl;
		cout << "Valor: " << temp << endl;
	}
}

double g1(double x, double y) {
	return sqrt((x * y + 5 * x - 1) / 2);
}

double g2(double x, double y) {
	return sqrt(x + 3 * log(x));
}

void picard_peano_sys() {
	double x=1, y=1;
	int counter = 0;
	double temp = g1(x, y);
	double tempy = g2(x, y);
	while (abs(temp - x) > pow(10,-3)&& abs(tempy - y) > pow(10, -3)) {
		x = temp;
		y = tempy;
		temp = g1(x, y);
		tempy = g2(x, y);
		counter++;
		cout << "Iter: " << counter << endl;
		cout << "Valor x: " << temp << "  \t Valor y: " << tempy << endl;
	}
}

double o(double x) {
	return exp(x)-x-5;
}

double o_d(double x) {
	return exp(x)-1;
}

void newton() {
	double a = 3;
	double x = -o(a) / o_d(a) + a;
	double temp = 0;
	int counter = 0;
	while (true) {
		temp = -o(x) / o_d(x) + x;
		cout << "iter n:" << counter << " x:" << x << endl;
		if (abs(x - temp) <= pow(10,-5)) {
			counter++;
			x = temp;
			cout << "iter n:" << counter << " x:" << x << endl;
			break;
		}
		x = temp;
		counter++;
	}
}

double f1(double x, double y) {
	return x*x-y-1.2;
}
double f2(double x, double y) {
	return -x+y*y-1;
}
double f1_x(double x, double y) {
	return 2*x;
}
double f2_x(double x, double y) {
	return -1;
}
double f1_y(double x, double y) {
	return -1;
}
double f2_y(double x, double y) {
	return +2 * y;
}
double J(double x, double y) {
	return f1_x(x, y) * f2_y(x, y) - f1_y(x, y) * f2_x(x, y);
}
void newtonsys() {
	double x=1, y=1;
	int counter = 0;
	double temp = x - (f1(x, y) * f2_y(x, y) - f2(x, y) * f1_y(x, y)) / J(x, y);
	double tempy = y - (f2(x, y) * f1_x(x, y) - f1(x, y) * f2_x(x, y)) / J(x, y);
	cout << "Iter: " << 0 << endl;
	cout << "Valor x: " << temp << "  \t Valor y: " << tempy << endl;
	while (abs(temp - x) > pow(10,-3)&& abs(tempy - y) > pow(10, -3)) {
		x = temp;
		y = tempy;
		temp = x - (f1(x, y) * f2_y(x, y) - f2(x, y) * f1_y(x, y)) / J(x, y);
		tempy = y - (f2(x, y) * f1_x(x, y) - f1(x, y) * f2_x(x, y)) / J(x, y);
		counter++;
		cout << "Iter: " << counter << endl;
		cout << "Valor x: " << temp << "  \t Valor y: " << tempy << endl;
	}
}

void externalError() {
	vector<vector<double>> matrix = { {0.1,0.5,3,0.25},{1.2,0.2,0.25,0.2},{-1,0.25,0.3,2},{2,0.00001,1,0.4} };
	double coefficient_error = 0.3;
	vector<double> err_list = { 0,0,0,0 };
	vector<double> sol_list = { 0.97263,-3.06443,0.32662,1.82038 };
	double div = 0, mul = 0;
	for (size_t i = 0; i < matrix.size(); i++)
	{
		err_list[i] = coefficient_error * (1 - sol_list[0] - sol_list[1] - sol_list[2] - sol_list[3]);	//com 4 acrescentar sol_list[3]
	}
	for (size_t i = 0; i < matrix.size(); i++)
	{
		div = matrix[i][i];
		for (size_t j = i; j < matrix.size(); j++)
		{
			matrix[i][j] /= div;
		}
		err_list[i] /= div;
		for (size_t j = i + 1; j < matrix.size(); j++)
		{
			mul = matrix[j][i];
			for (size_t k = i; k < matrix.size(); k++)
			{
				matrix[j][k] -= matrix[i][k] * mul;
			}
			err_list[j] -= err_list[i] * mul;
		}
	}
	for (size_t i = 0; i < matrix.size(); i++)
	{
		for (size_t j = 0; j < matrix.size(); j++)
		{
			cout << matrix[i][j] << "\t";
		}
		cout << "|" << err_list[i] << endl;
	}
	for (int i = matrix.size() + 1; i >= 0; i--)
	{
		for (size_t j = i + 1; j < matrix.size(); j++)
		{
			err_list[i] -= matrix[i][j] * err_list[j];
		}
	}
	cout << "Delta 1: " << err_list[0] << endl;
	cout << "Delta 2: " << err_list[1] << endl;
	cout << "Delta 3: " << err_list[2] << endl;
	cout << "Delta 4: " << err_list[3] << endl;
}

float r1(float x, float y, float z) {
	return (7 - y - z) / 3;
}
float r2(float x, float y, float z) {
	return (4 - x - 2 * z) / 4;
}
float r3(float x, float y, float z) {
	return (5 - 2 * y) / 5;
}

void Gauss_Jacobi() {
	float x=1, y=1, z=1;
	int counter = 0;
	float tempx, tempy, tempz;
	do {
		tempx = x; tempy = y; tempz = z;
		x = r1(tempx, tempy, tempz);
		y = r2(tempx, tempy, tempz);
		z = r3(tempx, tempy, tempz);
		counter++;
		cout << "Iter: " << counter << endl;
		cout << "Valor x: " << x << "  \t Valor y: " << y << "  \t Valor z: " << z << endl;
	} while (abs(tempx - x) > pow(10, -5) || abs(tempy - y) > pow(10, -5) || abs(tempz - z) > pow(10, -5));
}

void Gauss_Seidel() {
	float x=1, y=1, z=1;
	int counter = 0;
	float tempx, tempy, tempz;
	do {
		tempx = x; tempy = y; tempz = z;
		x = r1(tempx, tempy, tempz);
		y = r2(x, tempy, tempz);
		z = r3(x, y, tempz);
		counter++;
		cout << "Iter: " << counter << endl;
		cout << "Valor x: " << x << "  \t Valor y: " << y << "  \t Valor z: " << z << endl;
	} while (abs(tempx - x) > pow(10, -5) || abs(tempy - y) > pow(10, -5) || abs(tempz - z) > pow(10, -5));
}

double sty(double x) {
	return exp(1.5*x);
}

void TrapNorm() {
	double result = 0, h = 0.03125, x = 0, xf = 1; //mudar h, x, xf e FUNÇÃO ACIMA
	for (double i = x + h; i < xf; i += h) {
		result += sty(i);
	}
	result *= 2;
	result += sty(x) + sty(xf);
	result *= h * 0.5;
	cout << fixed<<setprecision(5)<<result;
}
void SimpsonNorm() {
	double result = 0, h = 0.03125, x = 1, xf = 1.5;
	int counter = 1;
	for (double i = x + h; i < xf; i += h) {
		if (counter % 2 == 0) {
			result += 2*sty(i);
		}
		else {
			result += 4*sty(i);
		}
		counter++;
		
	}
	result += sty(x) + sty(xf);
	result *= h;
	result /= 3;
	cout << fixed << setprecision(5) << result;
}

void Simpson() {
	vector<double> results = { 0,0,0 };
	double result, h = 1, err, dx = 0.25;
	int counter;
	vector<double> list = { 0.36,1.19,1.32,0.21,1.15,1.39,0.12,1.22,0.60 };
	for (int j = 0; j < 3; j++)
	{
		result = list[0];
		counter = 1;
		for (int i = (h / dx); i < list.size() - h / dx; i += (h / dx))
		{
			if (counter % 2 == 0)
				result += list[i] * 2;
			else
				result += list[i] * 4;
			counter++;
		}
		result += list[list.size() - 1];
		result *= h / 3;
		cout << "Step: " << h << "; Result: " << result << endl;
		results[j] = result;
		h /= 2;
	}
	err = (results[2] - results[1]) / 15;
	cout << "Error: " << err;

}

void Trap3d() {
	vector<vector<double>> xy = { {1.1,2.1,7.3 },  { 1.4,3.1,1.5 }, { 7.7,2.2,1.2 } };
	double result=0, ha = 1, hb = 1;
	result += xy[0][0] + xy[2][0] + xy[0][2] + xy[2][2];
	result += 2*(xy[1][0] + xy[0][1] + xy[1][2] + xy[2][1]);
	result += 4 * xy[1][1];
	result *= ha * hb;
	result /= 4;
	cout << "Result: " << result;
}

void Simpson3d() {
	vector<vector<double>> xy = { {1.1,2.1,7.3 },  { 1.4,3.1,1.5 }, { 7.7,2.2,1.2 } };
	double result = 0, ha = 1, hb = 1;
	result += xy[0][0] + xy[2][0] + xy[0][2] + xy[2][2];
	result += 4 * (xy[1][0] + xy[0][1] + xy[1][2] + xy[2][1]);
	result += 16 * xy[1][1];
	result *= ha * hb;
	result /= 9;
	cout << "Result: " << result;
}

float wx(float x, float y) {
	return -0.25*(y-37);
}

void Euler2d() {
	float x, y, h = 0.4;
	x = 5;
	y = 3;
	//float n = (x2 - x1) / h;
	for (int i = 1; i <= 2; i++) {
		y = y + h * wx(x, y);
		x = x + h;
		cout << "Iteracao " << i << endl;
		cout << "y: " << y << endl;
	}
}

double cy(double x, double y, double z) {
	return -exp(-0.5/(z+273))*y;
}
double cz(double x, double y, double z) {
	return 30* exp(-0.5 / (z + 273)) * y-0.5*(z-20);
}
void EulerOrd3() {
	double x, y, z, ztemp, h = 0.0625;
	x = 0;
	y = 2.5;
	z = 25;
	//float n = (x2 - x1) / h;
	cout << fixed << setprecision(5);
	for (int i = 0; i <9; i++) {
		
		cout << "Iteracao " << i << endl;
		cout << "x: " << x << endl;
		cout << "y: " << y << endl;
		cout << "z: " << z << endl;
		ztemp = z;
		z = z + h * cz(x, y, z);
		y = y + h * cy(x, y, ztemp);
		x = x + h;
	}
	
}

double rx(double x, double y) {
	return sin(2 * x) + sin(2 * y);
}

void RK4() {
	vector<double> sol_list = { 0,0,0 };
	double t = 1, C = 1, h = 0.125;
	double dC1, dC2, dC3, dC4;
	for (int i = 0; i <= 6; i++)
	{
		cout << "Iteracao " << i << endl;
		cout << "t: " << t << endl;
		cout << "C: " << C << endl;
		dC1 = h * rx(t, C);
		dC2 = h * rx(t + h / 2, C + dC1 / 2);
		dC3 = h * rx(t + h / 2, C + dC2 / 2);
		dC4 = h * rx(t + h, C + dC3);
		C += dC1 / 6 + dC2 / 3 + dC3 / 3 + dC4 / 6;
		t += h;
	}
}

double ey(double x, double y, double z)
{
	return -exp(-0.5 / (z + 273)) * y;
}

double ez(double x, double y, double z)
{
	return 30 * exp(-0.5 / (z + 273)) * y - 0.5 * (z - 20);
}

void RK43d() {
	double h = 0.25;
	double x = 0, y = 2.5, z = 25;
	cout << fixed << setprecision(5);
	double dY1, dY2, dY3, dY4, dZ1, dZ2, dZ3, dZ4;
	for (int i = 0; i < 3; i++) {
		cout << "iteracao " << i << endl;
		cout << "x: " << x << endl;
		cout << "y: " << y << endl;
		cout << "z: " << z << endl;
		dY1 = h * ey(x, y, z);
		dZ1 = h * ez(x, y, z);
		dY2 = h * ey(x + h / 2, y + dY1 / 2, z + dZ1 / 2);
		dZ2 = h * ez(x + h / 2, y + dY1 / 2, z + dZ1 / 2);
		dY3 = h * ey(x + h / 2, y + dY2 / 2, z + dZ2 / 2);
		dZ3 = h * ez(x + h / 2, y + dY2 / 2, z + dZ2 / 2);
		dY4 = h * ey(x + h, y + dY3, z + dZ3);
		dZ4 = h * ez(x + h, y + dY3, z + dZ3);
		z += dZ1 / 6 + dZ2 / 3 + dZ3 / 3 + dZ4 / 6;
		y += dY1 / 6 + dY2 / 3 + dY3 / 3 + dY4 / 6;
		x += h;
	}
}

double rx(double x) {
	return pow(x-5,2)+pow(x,4)	;
}

void aurea() {
	double x1 = 1, x2 = 3; //alterar as bordas e a FUNÇÃO ACIMA
	double x3, x4;
	double B, A;

	while (abs(x2 - x1) > pow(10, -3)) {
		B = (sqrt(5) - 1) / 2;
		A = pow(B, 2);
		x3 = x1 + A * (x2 - x1);
		x4 = x1 + B * (x2 - x1);
		cout << "x1: " << x1 << endl;
		cout << "x2: " << x2 << endl;
		cout << "x3: " << x3 << endl;
		cout << "x4: " << x4 << endl << endl;
		if (rx(x3) < rx(x4)) {
			x2 = x4;
		}
		else {
			x1 = x3;
		}
	}
	cout << "Min: ";
	if (rx(x1) < rx(x2)) {
		cout << rx(x1) << endl;
	}
	else {
		cout << rx(x2) << endl;
	}
}

float qx(float x, float y) {
	return -1.1*x*y+12*y+7*x*x-8*x;
}

void gradiente() {
	float x = 3, y = 1, tempx, tempy;
	float h = 0.1;
	while (true) {
		tempx = x - h * (-1.1*y+14*x-8);
		tempy = y - h * (-1.1*x+12);
		cout << "x: " << x << endl;
		cout << "y: " << y << endl;
		cout << "grad: " << qx(tempx,tempy) << endl << endl;
		if (abs(x - tempx) < pow(10, -3) && abs(y - tempy) < pow(10, -3)) {
			break;
		}
		if (qx(x, y) < qx(tempx, tempy)) {
			h /= 2;
		}
		else {
			h *= 2;
			x = tempx;
			y = tempy;
		}
		break;
	}
	cout << qx(x,y);	//3,6
}

void quadratic() {
	float x = 1, y = 1, tempx, tempy;
	while (true) {
		tempx = x - (4 * x - 2 * y) / 4;
		tempy = y - (2 * y - 2 * x - 6) / 4;
		if (abs(x - tempx) < pow(10, -3) && abs(y - tempy) < pow(10, -3)) {
			break;
		}
		x = tempx;
		y = tempy;
	}
	cout << x << " " << y;	//3,6
}

