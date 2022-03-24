
const double MASSA = 66.317;// Масса элемента в наших единицах (кг^-27)

double REBROKR=3.338339;//0.360276//2.287075;//3.338339;//4.968253

//Размеры области LX * LY * LZ
double LX=NUMKRIST_X * REBROKR;
double LY=NUMKRIST_Y * REBROKR;
double LZ=NUMKRIST_Z * REBROKR;
double HALFLX = LX/2;
double HALFLXSQR=HALFLX*HALFLX;
double VOLUME=LX*LY*LZ;

const double D = 3;//Число степеней свободы
const double KBOLTZMAN = 1.380648528; //Постоянная Больцмана
const double START_TEMP = 2.7315;//Стартовая температура системы
const double PREF_TEMP = 2.7315;//Предпочтительная температура системы для термостата
const double PREF_PRESSURE = 10.085;//Предпочтительное давление системы для баростата
const double SIGMA_Maxwell = sqrt(KBOLTZMAN*START_TEMP/MASSA);
const double Tau_Ber=1.0;
//Параметры потенциала Леннарда Джонса
const double EPSILON_LJ = 1.65401750; //Глубина потенциальной ямы
const double SIGMA_LJ = 0.3405;
const double SIGMA_LJ2=SIGMA_LJ*SIGMA_LJ;
const double RCUT = 2.5 * SIGMA_LJ; //Радиус обрезания потенциала

const double RCUT2 = RCUT*RCUT;

const double DELTA_T = 0.002;//0.01;//Шаг по времени

const double MT2= DELTA_T*DELTA_T/(2*MASSA);//Константа для расчета координат(t^2/2m)
const double MT=DELTA_T/(2*MASSA);//Константа для расчета скорости(t/2m)
//Константы для расчета потенциала Леннарда Джонса
const double EPSILON_LJ4= 4*EPSILON_LJ;
const double EPSILON_LJ24 = 24*EPSILON_LJ;

const double sigmarcut = pow(SIGMA_LJ2/RCUT2,3);
const double RCUT_POT = EPSILON_LJ4*((sigmarcut*sigmarcut)-sigmarcut);

const double T_CONST=(2/(D*KBOLTZMAN));//Константа для расчета температуры