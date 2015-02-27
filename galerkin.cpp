//#define PROG_WIN PROG_WIN //закомментировать данную строку,
						  //если работаем в среде Linux
#ifndef PROG_WIN
	#include "mpi.h"
#else
	#include "conio.h"
#endif
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;

stringstream strFileName;

const double mu22 = 1.01;
const double pi = 3.14159265358979323846;

float a = 0.0, //Левая нижняя точка на границе тела 
	  b = 1.0; //Правая верхняя точка на границе тела (короче, единичный куб R2 => единичный квадрат)

//Правая часть интегрального уравнения
double f0(float *X)
{
	//X[0] -- первая координата
	//X[1] -- вторая координата

	return sin(8*X[0]);
}

//Ядро. Функция от двух(из R2) переменных
double K (float *X, float *Y)
{
	//X[0] -- первая координата первой переменной
	//X[1] -- вторая координата первой переменной
	//Y[0] -- первая координата второй переменной
	//Y[1] -- вторая координата второй переменной

	return (
		1.0
		/
		(	4.0*pi
			*
			sqrt(((X[0]-Y[0])*(X[0]-Y[0])+(X[1]-Y[1])*(X[1]-Y[1])))
		)
			);
}


//Правая часть СЛАУ. *X -- левая нижняя точка кубика, h -- длина грани кубика
double intf0(float *X, int n, float h)
{
	//X[0] -- первая координата первой переменной
	//X[1] -- вторая координата первой переменной
	//Y[0] -- первая координата второй переменной
	//Y[1] -- вторая координата второй переменной

	float hx, hy;
	double summ = 0.0;

	hx = h/n;
	hy = hx; //Шаг интегрирования у нас одинаковый (только для правой части)

	float Xh[2];//Для текущего значения

	//Далее -- Метод прямоугольников с выбором средней точки

	Xh[0] = X[0];
	Xh[1] = X[1];
	Xh[0] += (hx/2.0);//Попадаем в среднюю точку интервала разбиения численного интегрирования
	Xh[1] += (hy/2.0);//Попадаем в среднюю точку интервала разбиения численного интегрирования

	int i, j;
	
	for(i = 0; i < n; i++)
	{// 1  for(i=0; i< n ;i++)
		for(j = 0; j < n; j++)
		{// 11  for(j=0; j<n ;j++)
			summ += f0(Xh);
			Xh[1] += hy;
		}// 11  for(j=0; j< n ;j++)

		//Подготовка к следующему шагу по первой координате

		//Сбиваем в первоначальное значение вторую координату
		Xh[1] = X[1];
		Xh[1] += (hy/2.0);

		Xh[0] += hx;
	
	}// 1  for(i=0;i< n ;i++)
	summ *= (hx*hy);
	return (summ);
}

//Четырехкратный интеграл по двум (из R2) переменным
double intK(float *X, float *Y, int n, float h)
{
	float hx1, hy1, hx2, hy2;
	double summ = 0.0;

	hx1 = h/n;
	hy1 = hx1; //Шаг у нас одинаковый

	hx2 = h/(n+1);//Для разноса точек от особенности
	hy2 = hx2; //Шаг у нас одинаковый
	
	float Xh[2];//Для текущего значения
	Xh[0] = X[0];
	Xh[1] = X[1];
	Xh[0] += (hx1/2.0);//Попадаем в среднюю точку интервала разбиения численного интегрирования
	Xh[1] += (hy1/2.0);//Попадаем в среднюю точку интервала разбиения численного интегрирования

	float Yh[2];//Для текущего значения по второй переменной
	Yh[0] = Y[0];
	Yh[1] = Y[1];
	Yh[0] += (hx2/2.0);//Попадаем в среднюю точку интервала разбиения численного интегрирования
	Yh[1] += (hy2/2.0);//Попадаем в среднюю точку интервала разбиения численного интегрирования


	int i, j, ii, jj;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j< n; j++)
		{
			for(ii = 0; ii < (n + 1); ii++) // Мы сделали разнос (шаг интергирования поменьше, а интервлов по второй переменной на единицу больше)
			{
				for(jj = 0; jj < (n + 1); jj++)
				{
					summ += K(Xh,Yh);
					Yh[1] += hy2;

				}
				
				//Подготовка к следующему шагу по первой координате второй переменной
				
				//Сбиваем в первоначальное значение вторую координату второй переменной
				Yh[1] = Y[1];
				Yh[1] += (hy2/2.0);
								
				Yh[0] += hx2;
			}
			
			//Подготовка к следующему шагу по второй координате первой переменной
			
			//Сбиваем в первоначальное значение первую координату второй переменной
			Yh[0] = Y[0];
			Yh[0] += (hx2/2.0);
			
			Xh[1] += hy1;
		}
		
		//Подготовка к следующему шагу по первой координате первой переменной
		
		//Сбиваем в первоначальное значение вторую координату первой переменной
		Xh[1] = X[1];
		Xh[1] += (hy1/2.0);
		
		Xh[0] += hx1;
	}
	summ*=(hx1*hx2*hy1*hy2);
	summ*=(mu22-1);
	return (summ);
}


int main (int argc, char* argv[])
{
	int numProc = 0; //номер процесса
	int groupsize = 1; //количество процессов
	int* displs;//смещение
	int* recvcounts;//число элементов получаемых от каждого процесса
	double startTime, endTime;

#ifndef PROG_WIN
	MPI_Init(&argc, &argv);		
	MPI_Comm_size(MPI_COMM_WORLD, &groupsize);
	MPI_Comm_rank(MPI_COMM_WORLD, &numProc);
		
	startTime = MPI_Wtime();
#endif
	// Количество интервалов разбиения по одной координате
	int N = 10;
	//Количество интервалов разбиения для численного интегрирования
	int n = 6;

	//Передаем количество интервалов разбиения в качестве аргументов
	//командной строки
	if (1 < argc)
	{
		 N = atoi(argv[1]);
	}
	if(2 < argc)
	{
		n = atoi(argv[2]);
	}

	//генерируем имя файла, в который будет записано решение: "N_n_groupsize.txt"
	strFileName << argv[1] << "_" << argv[2] << "_" << groupsize << ".txt";
	
	int NN = N*N; //кол-во базисных функций

	int ost;
	
	int NNNN_Local; //NNNN_Local - для массива А на каждом процессе
	int indexBeginNNNN_Local;
	NNNN_Local = NN*NN / groupsize;
	ost = NN*NN % groupsize;
	indexBeginNNNN_Local = numProc * NNNN_Local;
	if (numProc < ost)
	{
		NNNN_Local++;
		indexBeginNNNN_Local += numProc;
	}
	else
	{
		indexBeginNNNN_Local += ost;
	}

	int NN_Local; //для массива B на каждом процессе
	int indexBeginNN_Local;
	NN_Local = NN / groupsize;
	ost = NN % groupsize;
	indexBeginNN_Local = numProc * NN_Local;
	if (numProc < ost)
	{
		NN_Local++;
		indexBeginNN_Local += numProc;
	}
	else
	{
		indexBeginNN_Local += ost;
	}
	
	//Основная матрица системы	
	double* Matrix;
	//Вектор правой части (а после прямого хода и вектор решения)
	double* vectorRight;

	int* recvcountsNNNN;
	int* recvcountsNN;
	int* displsNNNN;
	int* displsNN;
	
	if (numProc == 0)
	{
		Matrix = new double [NN*NN];
		vectorRight = new double [NN];
		recvcountsNNNN = new int [groupsize];
		recvcountsNN = new int [groupsize];

		displsNNNN = new int [groupsize];
		displsNN = new int [groupsize];			

		for (int i = 0; i < groupsize; i++)
		{
			ost = NN*NN % groupsize;

			recvcountsNNNN[i] = NN*NN / groupsize;
			displsNNNN[i] = i * recvcountsNNNN[i];
			if (i < ost)
			{
				recvcountsNNNN[i]++;
				displsNNNN[i] += i;
			}
			else
			{
				displsNNNN[i] += ost;
			}

			recvcountsNN[i] = NN / groupsize;
			ost = NN % groupsize;
			displsNN[i] = i * recvcountsNN[i];
			if (i < ost)
			{
				recvcountsNN[i]++;
				displsNN[i] += i;
			}
			else
			{
				displsNN[i] += ost;
			}
		}
	}
	else
	{
		Matrix = new double [NNNN_Local];
		vectorRight = new double [NN_Local];
	}	

	int i1, j1; //Индексы для первой переменной (в методе Гаусса не используются)
	int i2, j2; //Индексы для второй переменной (в методе Гаусса не используются)

	//Для метода Гаусса
	//(для хранения индекса строки, где найден коэффициент c максимальным по модулю значением)
	int Imax; 	
	
	//Для метода Гаусса
	//(для поиска коэффициента c максимальным по модулю значением)
	double valueMax;
	double valueMMax;

	float h = (b-a) / N;
	
	float X[2]; //хранение нижней левой точки интересующего кубика

	int I = 0; //Общий индекс по двум координатам (одной переменной) для правой части СЛАУ
			//и для двух координат первой переменной в матрице

	//Сначала разберемся с вектором правой части для СЛАУ, потом заполним матрицу

	for (i1 = 0; i1 < N; i1++)
	{// 1  for (i1=0; i1<N; i1++)
		
		X[0] = h*i1;
				
		for (j1 = 0; j1 < N; j1++)
		{// 11  for (j1=0; j1<N; j1++)
				if(I >= indexBeginNN_Local && I < indexBeginNN_Local + NN_Local)
				{
					X[1] = h*j1;
					vectorRight[I-indexBeginNN_Local] = intf0(X, n, h);
				}
				I++;	
		}// 11  for (j1=0; j1<N; j1++)
	}// 1  for (i1=0; i1<N; i1++)


	float Y[2]; //хранение нижней левой точки интересующего кубика
				//по второй переменной
	
	I = 0;
	int J;//Общий индекс по двум координатам второй переменной

	int II = 0;

	//Теперь (после заполнения вектора правой части) заполняю матрицу (смещаемся по квадратикам
	//независимо по двум переменным, но пройти надо все пары квадратиков)
	for (i1 = 0; i1 < N; i1++)
	{// 1  for (i1=0; i1<N; i1++) Первая переменная первая координата
		
		X[0] = h*i1;
		
		for (j1 = 0; j1 < N; j1++)
		{// 11  for (j1=0; j1<N; j1++) Первая переменная вторая координата
			
			X[1] = h*j1;
			J = 0;
			for (i2 = 0; i2 < N; i2++)
			{// 111  for (i2=0; i2<N; i2++) Вторая переменная первая координата
				
				Y[0] = h*i2;
		
				for (j2 = 0; j2 < N; j2++)
				{// 111 1  for (j2=0; j2<N; j2++) Вторая переменная вторая координата
						
					Y[1] = h*j2;
					II = I*N*N+J;

					if(II >= indexBeginNNNN_Local && II < indexBeginNNNN_Local + NNNN_Local)
					{
						//Внимание! Здесь знак "минус"
						Matrix[II - indexBeginNNNN_Local]=-intK(X,Y,n,h);
						if(I == J)
						{// 111 11  if(I==J)
							//Добавляем "площадь", если индексы одинаковые
							Matrix[II - indexBeginNNNN_Local] += (h*h);
						}// 111 11  if(I==J)
			
					}
					J++;
				}// 111 1  for (j2=0; j2<N; j2++) Вторая переменная вторая координата
			}// 111  for (i2=0; i2<N; i2++) Вторая переменная первая координата
			I++;
		}// 11  for (j1=0; j1<N; j1++) Первая переменная вторая координата
	}// 1  for (i1=0; i1<N; i1++) Первая переменная первая координата
	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef PROG_WIN
			//сбор данных на нулевом процессе
			MPI_Gatherv(Matrix, NNNN_Local, MPI_DOUBLE, Matrix, recvcountsNNNN, displsNNNN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Gatherv(vectorRight, NN_Local, MPI_DOUBLE, vectorRight, recvcountsNN, displsNN, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//Метод Гаусса с выделением максимального значения по строке
	if (numProc == 0)
	{
		int i, j;
		int k;
		for(k = 0; k < NN; k++)
		{// 1  for(k=0; k<NN; k++)
	
			// Сначала ведущий элемент. Запоминаем его индекс строки
	
			Imax = k;
			valueMMax = ::abs(Matrix[k*NN+k]);

			//Поиск максимального по модулю элемента
			for (i = k+1; i < NN; i++)
			{// 11  for (i=k; i<NN; i++) бежим по каждой строчке
				if (::abs(Matrix[i*N*N+k]) > valueMMax)
				{
					valueMMax = ::abs(Matrix[i*N*N+k]);
					Imax = i;
				}
			}
			if (k < Imax)
			{// 11  if (k < Imax)
			
				//Значение valueMax нам уже не сильно нужно,
				//так как значение Imax уже нашли.
				//Используем valueMax уже для перестановки местами
				//коэффициентов матрицы и вектора правой части			
			
				for(j = k; j < NN; j++)
				{// 111  for(j=0; j< NN; j++)
					valueMMax = Matrix[Imax*NN+j];
					Matrix[Imax*NN+j] = Matrix[k*NN+j];
					Matrix[k*NN+j] = valueMMax;
				}// 111  for(j=0; j< NN; j++)

				valueMMax = vectorRight[Imax];
				vectorRight[Imax] = vectorRight[k];
				vectorRight[k] = valueMMax;
			}// 11  if (k < Imax)

			//Ну вот. Мaксимальный элемент, найденный в строке "k" или в строчках ниже строки "k"
			//теперь стоит на главной диагонали матрицы в строке "k"

			//Делаем унитреугольную матрицу. Новое значение диагонального элемента не интересует
			//(оно всегда будет равно единице), а вот старым значением пользуемся
			vectorRight[k] /= Matrix[k*NN+k];
			for (j = NN - 1; j > k ; j--)
			{
				Matrix[k*NN+j] /= Matrix[k*NN+k];
			}


			//Элементы под диагональю не трогаем, так как знаем, что там нули будут.
			//Поэтому начинаем с индекса "(k+1)", а не "k"
			for (i = (k+1); i < NN; i++)
			{
				vectorRight[i] -= (vectorRight[k] * Matrix[i*NN+k]);
		
				for (j= (NN -1); j > k; j--)
				{
					Matrix[i*NN+j] = Matrix[i*NN+j]-Matrix[k*NN+j]*Matrix[i*NN+k];
				}
				//j==k не интересует (если сделать эту операцию, то будет нуль)
			}

		}// 1  for(k=0; k<NN; k++)

			//Обратный ход (записываю результат)
	
			for (i= (NN-2); i > -1; i--)
			{
				for (j = (NN-1); j > k; j--)
				{
					vectorRight[i] -= (Matrix[i*NN+j] * vectorRight[j]);
				}
			}
		}	
	//Конец метода Гаусса с выделением максимального значения по строке

#ifndef PROG_WIN
	//время завершения расчетов
	endTime = MPI_Wtime();
#endif
	
	if (numProc == 0)
	{
		//Создаем файл куда будем выводить полученное решение
		ofstream outFile;
		outFile.open(strFileName.str().c_str(), ofstream::out|ofstream::trunc);
		if (outFile == NULL)
			return 0;
	
		I=0;
		for(int i = 0; i < N; i++)
		{// 1  for(i=0;i<N;i++)
			for(int j=0;j<N;j++)
			{// 11  for(j=0;j<N;j++)
				outFile << vectorRight[I] << " ";
				I++;
			}// 11  for(j=0;j<N;j++)
			outFile << endl;
		}// 1  for(i=0;i<N;i++)
		outFile.close();

		delete[] recvcountsNNNN;
		delete[] recvcountsNN ;

		delete[] displsNNNN;
		delete[] displsNN;
	
		cout << endTime - startTime;
	}
			
	delete []Matrix;
	delete []vectorRight;

#ifdef PROG_WIN
	getch();
#else
	MPI_Finalize();
#endif

	return 0;
}
