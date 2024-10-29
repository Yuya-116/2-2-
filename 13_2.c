/*�M�`���v�Z�v���O�����i���K13.2�@����ݒ�ۑ�j
by Noguchi Yuya 2024�N1��24��*/

#include<stdio.h>
#include<math.h>
#define NCVX 5         //������N
#define NX (NCVX +2)    //�m�ۂ���z��̌�
#define LENGTHX 0.50    //���_�̒��� 0.1m
#define RHO 8920.0      //���̖��x
#define LAMBDA 398.0    //���_�̔M�`����
#define C 385           //���̔�M
#define TEMP_INI 10.0   //�_�̏������x(��������)
#define TAIR 25.0      //���[�O�����̉��xTa 10��
#define HAIR 120.0     //���[���M�`�B�� h 120W
#define DT 10.0     //���ԍ��ݕ�
#define RATE_T 1.0E-5   //���x�̎��ԕω���臒l
#define MAX_TIME 60*60*24

void grid(double x[NX], int* ni)
{
    int I;
    double dx;

    *ni = NCVX + 1;          //��ԉE�[�̊i�q�_�ԍ�
    dx = LENGTHX / NCVX;        //���W�̒���L��CV�̌��Ŋ����CV�̉���
    x[0] = 0;           //���W���_
    x[1] = x[0] + dx / 2;        //���[�̍s�g�Ԋu�͒ʏ�̔���
    for (I = 2; I < *ni; I++) x[I] = x[I - 1] + dx;  //�ʏ�̊i�q�Ԋu
    x[*ni] = x[*ni - 1] + dx / 2;
    printf("���W�v�Z�I��\n");
}
void init(double qin[NX], double temp_new[NX], int ni)
{
    int i;

    qin[0] = 0;
    qin[1] = 30000;
    qin[2] = 30000;
    qin[3] = 60000;
    qin[4] = 90000;
    qin[5] = 100000;
    qin[6] = 200000;
    qin[7] = 400000;
    qin[8] = 400000;
    qin[9] = 400000;
    qin[10] = 400000;
    qin[11] = 0;

    for (i = 0; i <= ni; i++)temp_new[i] = TEMP_INI;     //���[�܂߂ď������x��10��
}
//���x�X�V�֐�copynew2old//
void copynew2old(double temp_new[NX], double temp_old[NX], int ni)
{
    int i;
    for (i = 0; i <= ni; i++)
    {
        temp_old[i] = temp_new[i];
    }
}
//���x�v�Z�֐�calculate//
void calculate(double x[NX], double qin[NX], double temp_new[NX], double temp_old[NX], int ni)
{
    int i;
    double dxE, dxW, dx, delta_temp;
    double dx01;

    dx01 = x[1] - x[0];
    dx = LENGTHX / NCVX;
    for (i = 1; i < ni; i++)
    {
        dxW = x[i] - x[i - 1];
        dxE = x[i + 1] - x[i];
        delta_temp = LAMBDA * DT / (RHO * C * dx) * ((temp_old[i + 1] - temp_old[i]) / dxE - (temp_old[i] - temp_old[i - 1]) / dxW) + DT * qin[i] / (RHO * C);
        temp_new[i] = temp_old[i] + delta_temp;
    }
    temp_new[0] = (LAMBDA / dx01 * temp_new[1] + HAIR * TAIR) / (LAMBDA / dx01 + HAIR);
    temp_new[ni] = temp_new[ni - 1];
}

//���ʏo�͊֐�output2file//
void output2file(double x[NX], double temp_new[NX], int ni)
{
    int i;
    FILE* pfile;
    errno_t err;
    err = fopen_s(&pfile, "data1.csv", "w");

    for (i = 0; i <= ni; i++)
    {
        fprintf(pfile, "%9.5e,%10.5f\n", x[i] / LENGTHX, temp_new[i]);
    }
    fclose(pfile);
}

int main(void)
{
    double t;         // ���ԃX�e�b�v
    double qin[NX], x[NX];      // ���M�ʂ�CV�̍��W�l�̔z��
    double temp_new[NX], temp_old[NX];   // �V�����x�̔z��
    int ni;          // ��ԉE�[�̊i�q�_�̔ԍ�
    FILE* pfile;
    errno_t err;
    err = fopen_s(&pfile, "data2.csv", "w");
    grid(x, &ni);        // (1)�i�q�_�̍��W�lx[]�̌v�Z
    init(qin, temp_new, ni);     // (2)�������x,���M�ʂ̐ݒ�
    t = 0.0;         // �v�Z�J�n�����ݒ� t=0
    printf("�v�Z�J�n\n");
    do {          // ����ԂɒB����܂Ŏ��Ԑi�s���J��Ԃ�
        t += DT;        // (3)���������֐i�߂�

        copynew2old(temp_new, temp_old, ni); // (4)�Â����x�̍X�V (new��old�ɃR�s�[)
        calculate(x, qin, temp_new, temp_old, ni);
        printf("%10.2lf\t%12.4lf\n", t, temp_new[NCVX / 2]);
        fprintf(pfile, "%9.5e,%10.5f\n", t, temp_new[5]);//�ۑ�Q�͂�����temp_new[ni]��ni��0����11�ɕύX����B��Ɠ��l�ɘA���Ńf�o�b�N���čs����OK�Ō��CSV�m�F����Α��v//

    } while (t < MAX_TIME);      // ���x�ω�����RATE_T�ȏ�̂������͌v�Z���J��Ԃ�
    printf("�v�Z�I��\n");
    output2file(x, temp_new, ni);    // (6)�v�Z���ʂ��t�@�C���ɏ����o��
    printf("�t�@�C���o�͏I��\n");
    fclose(pfile);
    return 0;
}