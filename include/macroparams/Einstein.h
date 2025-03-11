#if defined(MACROPAR_ON) && defined(MACROPAR_EIN)
double CxCopy[PARTICLENUMBER] = {}, CyCopy[PARTICLENUMBER] = {}, CzCopy[PARTICLENUMBER] = {};

double Cx0[PARTICLENUMBER] = {}, Cy0[PARTICLENUMBER] = {}, Cz0[PARTICLENUMBER] = {};
double Vx0[PARTICLENUMBER] = {}, Vy0[PARTICLENUMBER] = {}, Vz0[PARTICLENUMBER] = {};

double EinD = 0;
double EinDNoTime = 0;
double EinVisc = 0;
double EinViscTconst = 0;
double MQVisc = 0;
double MQViscTconst = 0;

double EinViscNoTime = 0;
double EinViscTconstNoTime = 0;
double MQViscNoTime = 0;
double MQViscTconstNoTime = 0;
void RangeInit(int thread, int *i1, int *i2);

void initEin() { // Инициализация для формулы Ейнштейна
    for (int i = 0; i < PARTICLENUMBER; i++) {
        // Копируем координаты
        CxCopy[i] = Cx[i];
        CyCopy[i] = Cy[i];
        CzCopy[i] = Cz[i];
        // Определяем начальные координаты
        Cx0[i] = Cx[i];
        Cy0[i] = Cy[i];
        Cz0[i] = Cz[i];

        Vx0[i] = Vx[i];
        Vy0[i] = Vy[i];
        Vz0[i] = Vz[i];
    }
}

void CoordVerleEin() { // Расчет координат по схеме Верле(Без ПГУ)
    for (int i = 0; i < PARTICLENUMBER; i++) {
        CxCopy[i] += Vx[i] * DELTA_T + Fx[i] * MT2;
        CyCopy[i] += Vy[i] * DELTA_T + Fy[i] * MT2;
        CzCopy[i] += Vz[i] * DELTA_T + Fz[i] * MT2;
    }
}

void getEinD(int step) { // Расчет диффузии по формуле Ейнштейна
    double r;
    double sumR = 0;
    for (int i = 0; i < PARTICLENUMBER; i++) {
        r = (CxCopy[i] - Cx0[i]) * (CxCopy[i] - Cx0[i]) + (CyCopy[i] - Cy0[i]) * (CyCopy[i] - Cy0[i]) +
            (CzCopy[i] - Cz0[i]) * (CzCopy[i] - Cz0[i]);
        sumR += r;
    }
    EinDNoTime = sumR / (2 * D * PARTICLENUMBER);
    EinD = EinDNoTime / (step * DELTA_T);
}
void getEinVisc(int step) { // Расчеи вязкости по формуле Ейнштейна
    double rvxy = 0, rvyz = 0, rvxz = 0;
    double sumRV;
    for (int i = 0; i < PARTICLENUMBER; i++) {
        rvxy += (CxCopy[i] * Vy[i] - Cx0[i] * Vy0[i]);
        rvyz += (CyCopy[i] * Vz[i] - Cy0[i] * Vz0[i]);
        rvxz += (CxCopy[i] * Vz[i] - Cx0[i] * Vz0[i]);
    }
    sumRV = rvxy * rvxy + rvyz * rvyz + rvxz * rvxz;
    EinViscNoTime = (MASSA * MASSA * sumRV) / (2 * 3 * KBOLTZMAN * T * VOLUME);
    EinViscTconstNoTime = (MASSA * MASSA * sumRV) / (2 * 3 * KBOLTZMAN * PREF_TEMP * VOLUME);

    EinVisc = EinViscNoTime / (step * DELTA_T);
    EinViscTconst = EinViscTconstNoTime / (step * DELTA_T);
}

void getMQVisc(int step) { // Расчеи вязкости по формуле Мак-Кьюри
    double rvxy = 0, rvyz = 0, rvxz = 0;
    double sumRV;
    for (int i = 0; i < PARTICLENUMBER; i++) {
        rvxy += (CxCopy[i] * Vy[i] - Cx0[i] * Vy0[i]) * (CxCopy[i] * Vy[i] - Cx0[i] * Vy0[i]);
        rvyz += (CyCopy[i] * Vz[i] - Cy0[i] * Vz0[i]) * (CyCopy[i] * Vz[i] - Cy0[i] * Vz0[i]);
        rvxz += (CxCopy[i] * Vz[i] - Cx0[i] * Vz0[i]) * (CxCopy[i] * Vz[i] - Cx0[i] * Vz0[i]);
    }
    sumRV = rvxy + rvyz + rvxz;
    MQViscNoTime = (MASSA * MASSA * sumRV) / (2 * 3 * KBOLTZMAN * T * VOLUME);
    MQViscTconstNoTime = (MASSA * MASSA * sumRV) / (2 * 3 * KBOLTZMAN * PREF_TEMP * VOLUME);

    MQVisc = MQViscNoTime / (step * DELTA_T);
    MQViscTconst = MQViscTconstNoTime / (step * DELTA_T);
}

void doStepEin(int step) { // Отдельный расчет шага
    CoordVerleEin();
    getEinD(step);
    getEinVisc(step);
    getMQVisc(step);
}

#endif