#include<stdio.h>
#include<conio.h>
#include<math.h>
#pragma warning(disable : 4996)
// 水の物性値,重力加速度
#define DENSIT 1000.0	// 密度 kg/m3
#define VISCOS 1.0E-3	// 粘性係数 Pa s
#define GRAVITY 9.81

	double calc_flow(double z0,double z1,double L,double d,double U0);

int main(int argc,char *argv[])
 {
	double z0=0.0;
 	double z1=0.5;	//	高低差 z1-z0=0.5 m
 	double L=1.0;	// パイプ長 1 m
 	double d=0.01;	// パイプ内径 10 mm
 	double U0=1;	// 流速の初期予測値 1 m/s
 	printf("u=%8.4lf [m/s] \n",calc_flow(z0,z1,L,d,U0));
	printf("Strike any key to terminate ");
 	getch();
 }


//	関数 calc_flow()
// 水面高低差，パイプ長，パイプ内径より流速を求める
	double calc_flow(double z0,double z1,double L,double d,double U0)
{
	int k;
	double Re,u,lambda;	// レイノルズ数，流速，λ

	u=U0;	// 初期予測値をuに代入
	for(k=0;k<10;k++){		// 10回繰り返す
		Re = (DENSIT * u * d) / VISCOS;
		if (Re < 2300)
		{
			lambda = 64 / Re;
			u = sqrt((2 * GRAVITY * (z1 - z0)) / (1 + ((lambda * L) / d)));
		}
		else
		{
			lambda = pow((-1.8) * log10(6.9 / Re), (-2));
			u = sqrt((2 * GRAVITY * (z1 - z0)) / (1 + ((lambda * L) / d)));
		}
	}
	return u;
}	// End of calc_flow()


