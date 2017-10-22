/*  PROGRAM: szum skorelowany
 *  WERSJA: 0.1
 *  Autor: Artur GRACKI mailto: arteq(at)arteq(dot)org
 *  OSTATNIA MODYFIKACJA: 2007/04/16 (pon) 12:25:21
 *  KEYWORDS: szum skorelowany, gauss, szum biały, kolorowy, proces stochastyczny, ornstein-uhlenbeck, FFT
 */

#include	<iostream>
#include	<fstream>
#include	<cmath>
#include	<cstdlib>
#include	<cstdlib>
#include	<cstdio>
#include	<ctime>
#include	"fft.h"

// ilość punktów do generowania, potęga dwójki
#define N 2048

// żeby nie pisać co chwilę deklarujemy globalnie
#define dwa_PI 6.2831853

using namespace std;

/* ====================================================================== */

// tablice na dane
float dane_X[N], dane_Y[N], dane_A[N];

// parametry początkowe, h ~ 0.09 << 1
double tau=10.0, h=200.0/N, hh=h;

// na potrzeby randoma
time_t  t; 

/* ====================================================================== */

// zwraca liczbę z rozkładu Gaussa, wykorzystując systemowy random
float gauss()
{
	double losowanko;
	losowanko = sqrt(2*fabs( log((rand()+1.0)/RAND_MAX) ) ) * cos(dwa_PI*(rand()+1.0)/RAND_MAX);

	return losowanko;
}

/* ====================================================================== */

// zapis tablicy z danymi do pliku
void zapis(char plik[], float dat_x[], float dat_y[], int ile)
{
	// Otwieramy plik do zapisu
	ofstream OUT(plik); 
	
	for(int i=0; i<ile; i++)
	{
	 OUT << dat_x[i] << "\t" << dat_y[i] << "\n";
	}

	// Zamykamy pliczek
	OUT.close();
}

/* ====================================================================== */

int main()
{

	// Odpalamy generator liczb losowych, żeby random był bardziej randomowy
	srand( time(&t) );
	
	// Warunki początkowe
	dane_X[0] = 0.0;
	dane_Y[0] = 0.0;
	dane_A[0] = 0.0;

	for(int i=2; i<N-2; i+=2)
	{
		// Losujemy liczbę z rozkładu Gaussa
		dane_Y[i] = gauss();

		// Do nieparzystych elementów tablicy zapisujemy zera
		dane_X[i+1]=0;
		dane_Y[i+1]=0;

		dane_A[i/2] = hh;
		dane_X[i] = dane_X[i-2] * exp(-h/tau) + sqrt(2.0/tau)/2.0 * (1.0-exp(-2.0*h/tau))*gauss();
		hh += h;
	}

	// Liczymy FFT dla wygenerowanych danych
	fft(dane_X-1, N/2, 1);
	fft(dane_Y-1, N/2, 1);
	
	// Liczymy widmo mocy (kwadraty wsp)
	for(int i=0; i<N; i+=2)
	{
		dane_X[i] *= dane_X[i];
		dane_Y[i] *= dane_Y[i];
		dane_X[i/2] = dane_X[i];
		dane_Y[i/2] = dane_Y[i];
	}   

	// Zapisujemy dane do plików
	zapis("szuum.dat", dane_A, dane_X, N/2);
	zapis("gauss.dat", dane_A, dane_Y, N/2);
	
	return 0;
}
