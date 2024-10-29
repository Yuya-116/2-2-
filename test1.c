#include<stdio.h>
#include<conio.h>
#include<math.h>
#pragma warning(disable : 4996)
// ���̕����l,�d�͉����x
#define DENSIT 1000.0	// ���x kg/m3
#define VISCOS 1.0E-3	// �S���W�� Pa s
#define GRAVITY 9.81

	double calc_flow(double z0,double z1,double L,double d,double U0);

int main(int argc,char *argv[])
 {
	double z0=0.0;
 	double z1=0.5;	//	���፷ z1-z0=0.5 m
 	double L=1.0;	// �p�C�v�� 1 m
 	double d=0.01;	// �p�C�v���a 10 mm
 	double U0=1;	// �����̏����\���l 1 m/s
 	printf("u=%8.4lf [m/s] \n",calc_flow(z0,z1,L,d,U0));
	printf("Strike any key to terminate ");
 	getch();
 }


//	�֐� calc_flow()
// ���ʍ��፷�C�p�C�v���C�p�C�v���a��藬�������߂�
	double calc_flow(double z0,double z1,double L,double d,double U0)
{
	int k;
	double Re,u,lambda;	// ���C�m���Y���C�����C��

	u=U0;	// �����\���l��u�ɑ��
	for(k=0;k<10;k++){		// 10��J��Ԃ�
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


