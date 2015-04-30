#define _CRT_SECURE_NO_DEPRECATE
#include "stdio.h"
#include "math.h"
#include <iostream>
using namespace std;
 
 
double t = 0.0001;
double h = 0.03;
 
double time_current = 0.0;	//Текущее время
 
 
double a = 0.0;		//Левая граница
double b = 10.0;	//Правая граница
 
double Ua;	//Граничное условие слева
double Ub = 0.0;										//Граничное условие справа
 
int size = 1+(b-a)/h;
 
double function(double x)	//T(t = 0) - начальное распределение температуры
{
	return 0.0;
}
 
void set_conditions(double * u)		//задание начальные условия
{
	for (int i = 0; i*h <= (b-a); i++)
		*(u+i) = function(i*h+a);
}
 
void massive_set(double * a, double * b, double * c, double * f,double * u)		//Схема Кранка-Николсона
{
	Ua = 1.0 + 0.5*sin(2*3.1415926*time_current);	//По заданию, условие на правой границе меняется со временем
	for(int i = 0; i < size-1; i++)
		{
			*(a+i) = -1.0/(2*h*h);
			*(b+i) = -1.0/(h*h)-1.0/t;
			*(c+i) = -1.0/(2*h*h);
			*(f+i) = -*(u+i)*(1/t-1/(h*h))-*(u+i+1)/(2*h*h)-*(u+i-1)/(2*h*h);
		}
	*a = 0.0;
	*(c+size-1) = 0.0;	
	
	//Зададим условие на значение  слева (T=Ua, x=0) и значение функции справа (T=Ub, x=1)
 
	
	//Set left border conditions: - для значения на левой стороне
	*b = 1.0;
	*c = 0.0;
	*f = Ua;
	
	
	//Set right border condition
	*(a+size-1) = 0;
	*(b+size-1) = 1;
	*(f+size-1) = Ub;
	
 
	/*
	//Set right border conditions: (derivat.) - откомментить в случае надобности задания производной (двойная точность O(h*h)
	*(b+size-1)= 2.0/3.0;
	*(c+size-1) = 2.0/3.0 - h*h/(3.0*t);
	*(f+size-1) = *(u+1)*h*h/(3.0*t) - 2.0/3.0 * h * Ub;
	*/	
	return;
}
 
void massive_get(double * a, double * b, double * c, double * f, double * beta, double * z)		//Заполнение массивов бета, альфа и Z
{
	*beta =  *c / *b;
	*z = *f / *b;
 
	for(int i = 1; i < size-1; i++)
	{
		*(beta+i) = *(c+i) / (*(b+i) - (*(a+i))*(*(beta+i-1)));
		*(z+i) = (*(f+i) +(*(a+i))*(*(z+i-1)))/(*(b+i) - (*(a+i))*(*(beta+i-1)));
	}
	return;
}
 
void get_solution(double * beta, double * z,double * u)		//конечный пересчет
{
	*(u+size-1) = Ub;
 
	for (int i = size-2; i >= 0; i--)
		*(u+i) = *(beta+i) * (*(u+i+1)) + *(z+i);
 
	return;
}
 
void main()
{
	FILE * file = fopen("file.txt", "w");
	
	double * u = new double[size * sizeof(double)];
	double * u_next = new double[size * sizeof(double)];
	
	set_conditions(u); //Set start temperature distribution
	
	double * a = new double[size * sizeof(double)];	//Подсчет коэффициентов a,b,c из схемы и beta,z - промежуточных коэффициентов в неявной схеме
	double * b = new double[size * sizeof(double)];
	double * c = new double[size * sizeof(double)];
	double * f = new double[size * sizeof(double)];
 
	double * beta = new double[size * sizeof(double)];
	double * z = new double[size * sizeof(double)];
 
 
	double time = 10.0;		//До какого момента времени считать решение
	
	while (time_current < time)
	{
		time_current = time_current + t;
		massive_set(a,b,c,f,u);	//Задание массивов a,b,c,f
		massive_get(a,b,c,f,beta,z);	//Заполнение промежуточных коэффициентов
		get_solution(beta,z,u_next);	//Подсчет решения
		for (int i = 0; i < size; i++)	//Сохранение решения в u_next
			*(u+i) = *(u_next + i);
 
		fprintf(file, "%e %e %e\n",time_current,*u, *(u+50));	//Запишем в файл значение функции в момент времени time_current в точках x=0 (*u)
																//и точке x = 0.5 = h*50 (*(u+50) - пятидесятый элемент массива)
	}
 
	/*for (int i = 0; i < size; i++)				- если надо вывести распределение температуры в конечный момент времени
		fprintf(file, "%e %e\n", i*h, *(u+i));*/
 
	return;
}
