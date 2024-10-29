/*熱伝導計算プログラム（演習13.2　自主設定課題）
by Noguchi Yuya 2024年1月24日*/

#include<stdio.h>
#include<math.h>
#define NCVX 5         //分割数N
#define NX (NCVX +2)    //確保する配列の個数
#define LENGTHX 0.50    //銅棒の長さ 0.1m
#define RHO 8920.0      //銅の密度
#define LAMBDA 398.0    //銅棒の熱伝導率
#define C 385           //銅の比熱
#define TEMP_INI 10.0   //棒の初期温度(初期条件)
#define TAIR 25.0      //左端外部流体温度Ta 10℃
#define HAIR 120.0     //左端部熱伝達率 h 120W
#define DT 10.0     //時間刻み幅
#define RATE_T 1.0E-5   //温度の時間変化の閾値
#define MAX_TIME 60*60*24

void grid(double x[NX], int* ni)
{
    int I;
    double dx;

    *ni = NCVX + 1;          //一番右端の格子点番号
    dx = LENGTHX / NCVX;        //座標の長さLをCVの個数で割るとCVの横幅
    x[0] = 0;           //座標原点
    x[1] = x[0] + dx / 2;        //左端の行使間隔は通常の半分
    for (I = 2; I < *ni; I++) x[I] = x[I - 1] + dx;  //通常の格子間隔
    x[*ni] = x[*ni - 1] + dx / 2;
    printf("座標計算終了\n");
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

    for (i = 0; i <= ni; i++)temp_new[i] = TEMP_INI;     //両端含めて初期温度は10℃
}
//温度更新関数copynew2old//
void copynew2old(double temp_new[NX], double temp_old[NX], int ni)
{
    int i;
    for (i = 0; i <= ni; i++)
    {
        temp_old[i] = temp_new[i];
    }
}
//温度計算関数calculate//
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

//結果出力関数output2file//
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
    double t;         // 時間ステップ
    double qin[NX], x[NX];      // 加熱量とCVの座標値の配列
    double temp_new[NX], temp_old[NX];   // 新旧温度の配列
    int ni;          // 一番右端の格子点の番号
    FILE* pfile;
    errno_t err;
    err = fopen_s(&pfile, "data2.csv", "w");
    grid(x, &ni);        // (1)格子点の座標値x[]の計算
    init(qin, temp_new, ni);     // (2)初期温度,加熱量の設定
    t = 0.0;         // 計算開始時刻設定 t=0
    printf("計算開始\n");
    do {          // 定常状態に達するまで時間進行を繰り返し
        t += DT;        // (3)時刻を一つ先へ進める

        copynew2old(temp_new, temp_old, ni); // (4)古い温度の更新 (newをoldにコピー)
        calculate(x, qin, temp_new, temp_old, ni);
        printf("%10.2lf\t%12.4lf\n", t, temp_new[NCVX / 2]);
        fprintf(pfile, "%9.5e,%10.5f\n", t, temp_new[5]);//課題２はここのtemp_new[ni]のniを0から11に変更する。上と同様に連続でデバックして行けばOK最後にCSV確認すれば大丈夫//

    } while (t < MAX_TIME);      // 温度変化率がRATE_T以上のあいだは計算を繰り返す
    printf("計算終了\n");
    output2file(x, temp_new, ni);    // (6)計算結果をファイルに書き出す
    printf("ファイル出力終了\n");
    fclose(pfile);
    return 0;
}