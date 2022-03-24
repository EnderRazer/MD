//Частицы
//Координаты частиц
double *Cx;
double *Cy;
double *Cz;
//Скорости частиц
double *Vx;
double *Vy;
double *Vz;
//Сила действующая на частицы
double *Fx;
double *Fy;
double *Fz;

double *Fxnew;
double *Fynew;
double *Fznew;

//Глобальные переменные системы
double Epot1=0,Ekin1=0,Eterm1=0,Eint1=0,E1=0;//Энергии системы на 1 частицу
double E_av,Ekin_av,Epot_av;//Энергии системы усредненные по шагам
double sumVirials[3][3];//Сумма вириалов для 3-мерной задачи
double sumMV[3][3];
double VCM[3];//Скорость центра масс системы
double T,P;//Температура и Давление системы
double T_av,P_av;//Температура и Давление системы усредненные по шагам
double P_tensors[3][3];//Тензоры давления
double PV=0;//Импульс системы

double *V0x;
double *V0y;
double *V0z;
double Cdiff;//Коэффициент диффузии
double DiffIntegral=0;
double Diff_av=0;

double sumMV0[3][3];
double sumVirials0[3][3];

double Zkk,ZkkIntegral,Zkk_av=0;
double Zkp,ZkpIntegral,Zkp_av=0;
double Zpk,ZpkIntegral,Zpk_av=0;
double Zpp,ZppIntegral,Zpp_av=0;

double Zkkxy_av=0,Zkkyz_av=0,Zkkzx_av=0;
double Zkpxy_av=0,Zkpyz_av=0,Zkpzx_av=0;
double Zpkxy_av=0,Zpkyz_av=0,Zpkzx_av=0;
double Zppxy_av=0,Zppyz_av=0,Zppzx_av=0;

double p0[3];
double CVisc;
double ViscIntegral=0;
double Visc_av=0;


FILE *total_data;//Главный файл вывода программы
FILE *diff_data;
FILE *el_data;
FILE *p_data;
FILE *z_data;
const int BACKUP_FREQ = 100;//Частота бекапа
time_t start; 
time_t endt;
int threads=0;
int nstep=0;