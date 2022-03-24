//Количество шагов
const int NSTEPS = 2000000;

//Количество кристалических решеток по осям
const int NUMKRIST_X = 10; 
const int NUMKRIST_Y = NUMKRIST_X;
const int NUMKRIST_Z = NUMKRIST_X;

//Число частиц
const int PARTICLENUMBER=NUMKRIST_X*NUMKRIST_Y*NUMKRIST_Z;
const int TYPE=1;//(0-для 2х частиц, 1-примитивная кубическая решетка, 2-ГЦК)
const bool PGU=true;//(0 - без ПГУ, 1- с ПГУ)
const bool THERMOSTAT=false;//Использование термостата
const bool BACKUP=true;
const int BACKUPSTEP=2999900; //Начальный шаг с бекапа
