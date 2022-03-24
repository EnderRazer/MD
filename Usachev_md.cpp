#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <assert.h>
#include <omp.h>
#include "malloc.h"


#include "params.h"
#include "gas_constants.h"
#include "global_var.h"
#include "Speed_generation.h"
#include "start_cond.h"
#include "Backup_func.h"
#include "Memory.h"

using namespace std;
//Функции
double PotLJ(double r)//Вычисление потенциала Леннарда-Джонса(используем r^2)
{
    double sigmar = (SIGMA_LJ2/r)*(SIGMA_LJ2/r)*(SIGMA_LJ2/r);
    return EPSILON_LJ4*((sigmar*sigmar)-sigmar);
}

double FPotLJ(double r)//Вычисление силы потенциала(используем r^2)
{
    double sigmar = (SIGMA_LJ2/r)*(SIGMA_LJ2/r)*(SIGMA_LJ2/r);
    return EPSILON_LJ24*(2*(sigmar*sigmar) - sigmar);
}
//
void CoordVerle()//Расчет координат по схеме Верле
{
        for(int i=0;i<PARTICLENUMBER;i++){
                Cx[i] += Vx[i]*DELTA_T + Fx[i]*MT2;
                Cy[i] += Vy[i]*DELTA_T + Fy[i]*MT2;
                Cz[i] += Vz[i]*DELTA_T + Fz[i]*MT2;
                //ПГУ
                if(PGU){
                    if(Cx[i]>=LX) Cx[i]-=LX;
                    if(Cx[i]<0) Cx[i]+=LX;

                    if(Cy[i]>=LY) Cy[i]-=LY;
                    if(Cy[i]<0) Cy[i]+=LY;

                    if(Cz[i]>=LZ) Cz[i]-=LZ;
                    if(Cz[i]<0) Cz[i]+=LZ;
                }
        }
}

void outputInFile()//Запись данных в файл
{
    fprintf(el_data,"%10.6f;%10.6f;%10.6f;%10.6f;\n",LX,LY,LZ,VOLUME);
    fprintf(total_data,"%10.24f;%10.24f;%10.24f;%10.24f;%10.24f;%10.24f;%10.24f;;%10.24f;%10.24f;%10.24f;%10.24f;%10.24f;%10.24f;\n",T,P,Ekin1,Epot1,Eint1,Eterm1,E1,T_av/nstep,P_av/nstep,Ekin_av/nstep,Epot_av/nstep,E_av/nstep,PV);
    fprintf(diff_data,"%10.8f;%10.8f;%10.8f;%10.8f;\n",Cdiff,Diff_av,CVisc,Visc_av);
    fprintf(p_data,"%10.8f;%10.8f;%10.8f\n",P_tensors[0][1],P_tensors[1][2],P_tensors[2][0]);



}
void RangeInit(int thread, int *i1, int *i2){//Определение диапазона потока
    	int range = PARTICLENUMBER/threads;
	*i1 = thread*range;
	if(thread<threads-1) *i2=*i1+range; else *i2=PARTICLENUMBER;
}

void ForceJob(int nt){
	int i1, i2;//Диапазон индексов потока
	RangeInit(nt,&i1,&i2);
	//Локальные переменные каждого потока
        double rVecx=0,rVecy=0,rVecz=0,r2=0;//Векторное расстояние 
        double U=0,FU=0;//Потенциал и сила потенциала
        double localFx=0,localFy=0,localFz=0;//Вектор силы
		double localVirial[3][3]={0};
		double localEpot=0;
        double VCx=0,VCy=0,VCz=0;//Радиус-вектор виртуальной частицы
        for(int i=i1;i<i2;i++){
            for(int j=0;j<PARTICLENUMBER;j++){
                if(j!=i){
                    rVecx = Cx[i]-Cx[j];
                    rVecy = Cy[i]-Cy[j];
                    rVecz = Cz[i]-Cz[j];
                    r2 = (rVecx*rVecx) + (rVecy*rVecy) + (rVecz*rVecz);//Расстояние между частицами в квадрате
                    if(r2>=HALFLXSQR)//Если расстояние больше половины бокса - отражаем как виртуальную
                    {
                        //Отражение по X
                        if(abs(rVecx)>=HALFLX)//если изменение по X больше половины бокса отражаем по Х
                        {
                            if(Cx[j]<=HALFLX)
                            {
                                VCx=Cx[j]+LX;
                            }else{
                                VCx=Cx[j]-LX;
                            }
                        }else{
                            VCx=Cx[j];
                        }
                        //Отражение по Y
                        if(abs(rVecy)>=HALFLX)//если изменение по Y больше половины бокса отражаем по Y
                        {
                            if(Cy[j]<=HALFLX)
                            {
                                VCy=Cy[j]+LY;
                            }else{
                                VCy=Cy[j]-LY;
                            }
                        }else{
                            VCy=Cy[j];
                        }
                        //Отражение по Z
                        if(abs(rVecz)>=HALFLX)//если изменение по Z больше половины бокса отражаем по Z
                        {
                            if(Cz[j]<=HALFLX)
                            {
                                VCz=Cz[j]+LZ;
                            }else{
                                VCz=Cz[j]-LZ;
                            }
                        }else{
                            VCz=Cz[j];
                        }
                        rVecx = Cx[i]-VCx;
                        rVecy = Cy[i]-VCy;
                        rVecz = Cz[i]-VCz;
                        r2 = (rVecx*rVecx) + (rVecy*rVecy) + (rVecz*rVecz);
                    }
                    if(r2<=RCUT2){//Учет обрезания потенциала
                        //Вычисление потенциала Леннарда-Джонса(U(r))(Со сдвигом при обрезании потенциала) 
                        U = PotLJ(r2)-RCUT_POT;
                        localEpot += U;
                        //Вычисление силы потенциала(U`(r))(Используем r^2)
                        FU = FPotLJ(r2);
                        //Вычисление вектора силы
                        localFx = FU * rVecx /r2;
                        localFy = FU * rVecy /r2;
                        localFz = FU * rVecz /r2;
						//Вычисление вириалов
                        localVirial[0][0]+=rVecx*localFx;
                        localVirial[0][1]+=rVecx*localFy;
                        localVirial[0][2]+=rVecx*localFz;

                        localVirial[1][0]+=rVecy*localFx;
                        localVirial[1][1]+=rVecy*localFy;
                        localVirial[1][2]+=rVecy*localFz;

                        localVirial[2][0]+=rVecz*localFx;
                        localVirial[2][1]+=rVecz*localFy;
                        localVirial[2][2]+=rVecz*localFz;

                        Fxnew[i] += localFx;
                        Fynew[i] += localFy;
                        Fznew[i] += localFz;
			
                    }

                }
            }//end j
        }//end i
		#pragma omp critical
		{
            Epot1 += localEpot;
                        
            sumVirials[0][0]+=localVirial[0][0];
            sumVirials[0][1]+=localVirial[0][1];
            sumVirials[0][2]+=localVirial[0][2];

            sumVirials[1][0]+=localVirial[1][0];
            sumVirials[1][1]+=localVirial[1][1];
            sumVirials[1][2]+=localVirial[1][2];

            sumVirials[2][0]+=localVirial[2][0];
            sumVirials[2][1]+=localVirial[2][1];
            sumVirials[2][2]+=localVirial[2][2];
		}//end critical
}
void ForceCalc()//Вычисление силы, вириалов и потенциальной энергии
{
    #pragma omp parallel num_threads(threads)
    {
		int thread = omp_get_thread_num();//Номер потока
		ForceJob(thread);
    }//end parallel
}

void VelocityCalc()//Расчет скорости молекулы
{
    for(int i=0;i<PARTICLENUMBER;i++){
        Vx[i] = Vx[i] + (Fx[i]+Fxnew[i])*MT;
        Vy[i] = Vy[i] + (Fy[i]+Fynew[i])*MT;
        Vz[i] = Vz[i] + (Fz[i]+Fznew[i])*MT;
    }
}

void getVCM()//Расчет скорости центра масс системы
{
    VCM[0]=0,VCM[1]=0,VCM[2]=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        VCM[0]+=Vx[i];
        VCM[1]+=Vy[i];
        VCM[2]+=Vz[i];
    }
    VCM[0]/=PARTICLENUMBER;
    VCM[1]/=PARTICLENUMBER;
    VCM[2]/=PARTICLENUMBER;
}

double getAvgEkin()//Расчет кинетичетичекой энергии системы
{
    double Ekin=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        Ekin+=MASSA*(Vx[i]*Vx[i]+Vy[i]*Vy[i]+Vz[i]*Vz[i])/2;
    };
    return Ekin/PARTICLENUMBER;
}

double getAvgEterm()//Расчет тепловой энергии системы
{
    getVCM();
    double Eterm=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        Eterm+=MASSA*((Vx[i]-VCM[0])*(Vx[i]-VCM[0])+(Vy[i]-VCM[1])*(Vy[i]-VCM[1])+(Vz[i]-VCM[2])*(Vz[i]-VCM[2]))/2;
    };
    return Eterm/PARTICLENUMBER;
}

double getTemp()//Расчет температуры 1 молекулы системы
{
    return Eterm1*T_CONST;
}

double PressureCalc()//Расчет тензоров давления и давления системы 
{
    getVCM();
    for(int i=0;i<3;i++)
	for(int j=0;j<3;j++)
		sumMV[i][j]=0;
    
    for(int i=0;i<PARTICLENUMBER;i++){
        sumMV[0][0]+=(Vx[i]-VCM[0])*(Vx[i]-VCM[0]);//xx
        sumMV[0][1]+=(Vx[i]-VCM[0])*(Vy[i]-VCM[1]);//xy
        sumMV[0][2]+=(Vx[i]-VCM[0])*(Vz[i]-VCM[2]);//xz

        sumMV[1][0]+=(Vy[i]-VCM[1])*(Vx[i]-VCM[0]);//yx
        sumMV[1][1]+=(Vy[i]-VCM[1])*(Vy[i]-VCM[1]);//yy
        sumMV[1][2]+=(Vy[i]-VCM[1])*(Vz[i]-VCM[2]);//yz

        sumMV[2][0]+=(Vz[i]-VCM[2])*(Vx[i]-VCM[0]);//zx
        sumMV[2][1]+=(Vz[i]-VCM[2])*(Vy[i]-VCM[1]);//zy
        sumMV[2][2]+=(Vz[i]-VCM[2])*(Vz[i]-VCM[2]);//zz
    }
    sumMV[0][0]*=MASSA;
    sumMV[0][1]*=MASSA;
    sumMV[0][2]*=MASSA;

    sumMV[1][0]*=MASSA;
    sumMV[1][1]*=MASSA;
    sumMV[1][2]*=MASSA;

    sumMV[2][0]*=MASSA;
    sumMV[2][1]*=MASSA;
    sumMV[2][2]*=MASSA;

    P_tensors[0][0] = (sumMV[0][0] + 0.5 * sumVirials[0][0]) / VOLUME;
    P_tensors[0][1] = (sumMV[0][1] + 0.5 * sumVirials[0][1]) / VOLUME;
    P_tensors[0][2] = (sumMV[0][2] + 0.5 * sumVirials[0][2]) / VOLUME;

    P_tensors[1][0] = (sumMV[1][0] + 0.5 * sumVirials[1][0]) / VOLUME;
    P_tensors[1][1] = (sumMV[1][1] + 0.5 * sumVirials[1][1]) / VOLUME;
    P_tensors[1][2] = (sumMV[1][2] + 0.5 * sumVirials[1][2]) / VOLUME;

    P_tensors[2][0] = (sumMV[2][0] + 0.5 * sumVirials[2][0]) / VOLUME;
    P_tensors[2][1] = (sumMV[2][1] + 0.5 * sumVirials[2][1]) / VOLUME;
    P_tensors[2][2] = (sumMV[2][2] + 0.5 * sumVirials[2][2]) / VOLUME;
    //Расчет давления по XX,YY,ZZ компонентам
    return (P_tensors[0][0]+P_tensors[1][1]+P_tensors[2][2])/3;
}

void Zcalc(){
	double Zkkxy,Zkkyz,Zkkzx,ZkkxyIntegral,ZkkyzIntegral,ZkkzxIntegral;
	double Zkpxy,Zkpyz,Zkpzx,ZkpxyIntegral,ZkpyzIntegral,ZkpzxIntegral;
	double Zpkxy,Zpkyz,Zpkzx,ZpkxyIntegral,ZpkyzIntegral,ZpkzxIntegral;
	double Zppxy,Zppyz,Zppzx,ZppxyIntegral,ZppyzIntegral,ZppzxIntegral;

	Zkkxy=sumMV0[0][1]*sumMV[0][1];
	Zkkyz=sumMV0[1][2]*sumMV[1][2];
	Zkkzx=sumMV0[2][0]*sumMV[2][0];

	Zkpxy=sumMV0[0][1] * 0.5*sumVirials[0][1];
	Zkpyz=sumMV0[1][2] * 0.5*sumVirials[1][2];
	Zkpzx=sumMV0[2][0] * 0.5*sumVirials[2][0];

	Zpkxy=0.5*sumVirials0[0][1] * sumMV[0][1];
	Zpkyz=0.5*sumVirials0[1][2] * sumMV[1][2];
	Zpkzx=0.5*sumVirials0[2][0] * sumMV[2][0];

	Zppxy=0.5*sumVirials0[0][1] * 0.5*sumVirials[0][1];
	Zppyz=0.5*sumVirials0[1][2] * 0.5*sumVirials[1][2];
	Zppzx=0.5*sumVirials0[2][0] * 0.5*sumVirials[2][0];

	ZkkxyIntegral=Zkkxy*DELTA_T;
	Zkkxy_av+=ZkkxyIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	ZkkyzIntegral=Zkkyz*DELTA_T;
	Zkkyz_av+=ZkkyzIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	ZkkzxIntegral=Zkkzx*DELTA_T;
	Zkkzx_av+=ZkkzxIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	
	ZkpxyIntegral=Zkpxy*DELTA_T;
	Zkpxy_av+=ZkpxyIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	ZkpyzIntegral=Zkpyz*DELTA_T;
	Zkpyz_av+=ZkpyzIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	ZkpzxIntegral=Zkpzx*DELTA_T;
	Zkpzx_av+=ZkpzxIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	
	ZpkxyIntegral=Zpkxy*DELTA_T;
	Zpkxy_av+=ZpkxyIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	ZpkyzIntegral=Zpkyz*DELTA_T;
	Zpkyz_av+=ZpkyzIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	ZpkzxIntegral=Zpkzx*DELTA_T;
	Zpkzx_av+=ZpkzxIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	
	ZppxyIntegral=Zppxy*DELTA_T;
	Zppxy_av+=ZppxyIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	ZppyzIntegral=Zppyz*DELTA_T;
	Zppyz_av+=ZppyzIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	ZppzxIntegral=Zppzx*DELTA_T;
	Zppzx_av+=ZppzxIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);
	
	fprintf(z_data,"%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f;%10.8f\n",Zkkxy,Zkkyz,Zkkzx,Zkpxy,Zkpyz,Zkpzx,Zpkxy,Zpkyz,Zpkzx,Zppxy,Zppyz,Zppzx,ZkkxyIntegral,ZkkyzIntegral,ZkkzxIntegral,ZkpxyIntegral,ZkpyzIntegral,ZkpzxIntegral,ZpkxyIntegral,ZpkyzIntegral,ZpkzxIntegral,ZppxyIntegral,ZppyzIntegral,ZppzxIntegral);
	
	//Zkk=Zkkxy+Zkkyz+Zkkzx;
	//Zkp=Zkpxy+Zkpyz+Zkpzx;
	//Zpk=Zpkxy+Zpkyz+Zpkzx;
	//Zpp=Zppxy+Zppyz+Zppzx;

}
void Thermostat()//Реализация термостата Берендсена
{
    double lambda = sqrt(1+(DELTA_T/Tau_Ber)*((PREF_TEMP/T)-1));
    for(int i=0;i<PARTICLENUMBER;i++){
        Vx[i]*=lambda;
        Vy[i]*=lambda;
        Vz[i]*=lambda;
    }
}
void Barostat(){
    double hi = 1 - DELTA_T/Tau_Ber*(PREF_PRESSURE - P);
    double mu = pow(hi, 0.33333333);
    for (int i = 0; i < PARTICLENUMBER; i++){
        Cx[i] *= mu;
        Cy[i] *= mu;
        Cz[i] *= mu;
    }
    LX *= mu;
    LY *= mu;
    LZ *= mu;
    HALFLX=LX/2;
    HALFLXSQR=HALFLX*HALFLX; 
    VOLUME = LX * LY * LZ;
}

double Integral(double *f, double step,int size)//Решение интеграла методом трапеций
{
    double value=0.5*(f[0]+f[size-1]);
    for(int i=1;i<size-1;i++)
	value+=f[i];
    return value*step;
}

void DiffusionCDiff()//Расчет автокорелляционной функции скорости
{
    double CDiffx=0,CDiffy=0,CDiffz=0;
    for(int i=0;i<PARTICLENUMBER;i++){
        CDiffx+=Vx[i]*V0x[i];
        CDiffy+=Vy[i]*V0y[i];
        CDiffz+=Vz[i]*V0z[i];
    }

    Cdiff=(CDiffx+CDiffy+CDiffz);
}

void ViscosityCVisc()
{
    double CViscxy=P_tensors[0][1]*p0[0];
    double CViscyz=P_tensors[1][2]*p0[1];
    double CVisczx=P_tensors[2][0]*p0[2];
    
    CVisc=(CViscxy+CViscyz+CVisczx);

}
void SetOrigs(){
    for(int i=0;i<PARTICLENUMBER;i++){
        V0x[i]=Vx[i];
        V0y[i]=Vy[i];
        V0z[i]=Vz[i];
    }

    for(int i=0;i<3;i++)
	for(int j=0;j<3;j++){
	    sumMV0[i][j]=sumMV[i][j];
	    sumVirials0[i][j]=sumVirials[i][j];
	}

    p0[0]=P_tensors[0][1];
    p0[1]=P_tensors[1][2];
    p0[2]=P_tensors[2][0];
}

double PulseCalc(){
	double pulse=0;
	for(int i=0;i<PARTICLENUMBER;i++){
		pulse+=Vx[i]+Vy[i]+Vz[i];
	}
	return pulse*MASSA;
}

void ZeroStep(){
    Epot1=0,Ekin1=0,Eterm1=0,Eint1=0,E1=0;// Обнуление энергии на каждом шаге
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            sumVirials[i][j]=0;
        }
    }
    for(int i=0;i<PARTICLENUMBER;i++){
            Fxnew[i]=0;
            Fynew[i]=0;
            Fznew[i]=0;
    }
    ForceCalc();
    for(int i=0;i<PARTICLENUMBER;i++){
        Fx[i]=Fxnew[i];
        Fy[i]=Fynew[i];
        Fz[i]=Fznew[i];
    }
    //Расчет нужных параметров
    Epot1/=2*PARTICLENUMBER;//Расчет потенциальной энергии на 1 частицу
    Ekin1 = getAvgEkin();//Расчет кинетической энергии на 1 частицу
    Eterm1 = getAvgEterm();//Расчет тепловой энергии на 1 частицу
    E1=Ekin1+Epot1;//Расчет полной энергии на 1 частицу
    T = getTemp();//Расчет температуры системы
    P = PressureCalc();//Расчет давления системы
    PV=PulseCalc();
    cout<<"T= "<<T<<" P= "<<P<<" Pulse= "<<PV<<endl;
    nstep=1;
    SetOrigs();
}

void MD()//Основная функция расчетов МД
{
	Epot_av = 0.0;
	Ekin_av = 0.0;
	E_av = 0.0;
	T_av = 0.0;
	P_av = 0.0;
    Diff_av=0,Visc_av=0;
    ZeroStep();
    for(int n=1;n<NSTEPS;n++){
		Epot1=0,Ekin1=0,Eterm1=0,Eint1=0,E1=0;// Обнуление энергии на каждом шаге
        for(int i=0;i<3;i++){
            for(int j=0;j<3;j++){
                sumVirials[i][j]=0;
            }
        }
        for(int i=0;i<PARTICLENUMBER;i++){
            Fxnew[i]=0;
            Fynew[i]=0;
            Fznew[i]=0;
        }
        CoordVerle();//Расчет координат по схеме Верле
        ForceCalc();//Расчет сил и потенциальной энергии
        VelocityCalc();//Расчет скоростей
        
        for(int i=0;i<PARTICLENUMBER;i++){//Замена вектора силы предыдущего шага на силу текущего
            Fx[i] = Fxnew[i];
            Fy[i] = Fynew[i];
            Fz[i] = Fznew[i];
        }

	if(THERMOSTAT){//Перерасчет скоростей с учетом термостата
            Thermostat();
	    //Barostat();
        }
        //Бекап после расчета координат и скоростей
        if(n%BACKUP_FREQ==0){
            do_backup(n);
        }
        //Расчет нужных параметров
        Epot1/=2*PARTICLENUMBER;//Расчет потенциальной энергии на 1 частицу
        Epot_av+=Epot1;

        Ekin1 = getAvgEkin();//Расчет кинетической энергии на 1 частицу
        Ekin_av+=Ekin1;

        Eterm1 = getAvgEterm();//Расчет тепловой энергии на 1 частицу
        
        Eint1 = Eterm1+Epot1;//Расчет внутренней энергии на 1 частицу
        
        E1=Ekin1+Epot1;//Расчет полной энергии на 1 частицу
        E_av+=E1;
        
        T = getTemp();//Расчет температуры системы
        T_av+=T;
        
        P = PressureCalc();//Расчет давления системы
        P_av+=P;

        if(!THERMOSTAT){
            DiffusionCDiff();
            ViscosityCVisc();
			Zcalc();
            DiffIntegral=Cdiff*DELTA_T;
            Diff_av+=DiffIntegral/(3*PARTICLENUMBER);
            ViscIntegral=CVisc*DELTA_T;
            Visc_av+=ViscIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);

	    //ZkkIntegral=Zkk*DELTA_T;
	    //Zkk_av+=ZkkIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);

	    //ZkpIntegral=Zkp*DELTA_T;
	    //Zkp_av+=ZkpIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);

	    //ZppIntegral=Zpp*DELTA_T;
	    //Zpp_av+=ZppIntegral/(3*KBOLTZMAN*PREF_TEMP*VOLUME);

        }
        PV=PulseCalc();
        outputInFile();//Вывод в файл
        time(&endt);
          cout<<"\r"<<n<<"/"<<NSTEPS<<" Passed: "<<difftime(endt,start)<<" ET: "<<round(NSTEPS*difftime(endt,start)/n)<<flush;
	nstep++;
    }
}

int main(int argc, char* argv[]){
    if (argc > 1){
        sscanf(argv[1],"%d",&threads);
    }else{
        threads=1;
    }
    char filename[50];
    char filenameD[50];
    char filenameV[50];
    char filenameP[50];
    char filenameZ[50];

    sprintf(filename, "TPE_%06d_%10.4fK.txt",PARTICLENUMBER,PREF_TEMP);
    sprintf(filenameD, "Diff_%06d_%2.4fK.txt",PARTICLENUMBER,PREF_TEMP);
    sprintf(filenameV, "XYZ_%06d_%2.4fK.txt",PARTICLENUMBER,PREF_TEMP);
    sprintf(filenameP, "Tensors_%06d_%2.4fK.txt",PARTICLENUMBER,PREF_TEMP);
    sprintf(filenameZ, "Z_%06d_%2.4fK.txt",PARTICLENUMBER,PREF_TEMP);


    total_data = fopen(filename, "w");
    diff_data= fopen(filenameD,"w");
    el_data = fopen(filenameV, "w");
    p_data = fopen(filenameP, "w");
    z_data = fopen(filenameZ,"w");

    AllocArrays();
    if(BACKUP){
        restore_backup(BACKUPSTEP);
    }else{
        start_cond(TYPE,PGU);
    }
    cout<<threads<<endl;
    time(&start);
    MD();
    fclose(total_data);
    cout<<"Done in "<<difftime(endt, start)<<endl;
    return 0;
}
