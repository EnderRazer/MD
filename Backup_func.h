//Функции бекапа
#define OSIZE sizeof(int)
//Выводим любые данные в виде набора unsigned int в десятичной записи
static void myfwrite(void *ptr, int size, int n, FILE *file){
    int i;
    unsigned int *p = (unsigned int *)ptr;
    assert(size%OSIZE==0);
    size /= OSIZE;
    for(i=0;i<size*n;i++)
        fprintf(file, "%u ", p[i]);
    return;
}
//Читаем данные из набора unsigned int в десятичной записи
static void myfread(void *ptr, int size, int n, FILE *file){
    int i;
    unsigned int *p = (unsigned int *)ptr;
    assert(size%OSIZE==0);
    size /=OSIZE;
    for(i=0; i<size*n;i++){
        assert(1==fscanf(file, "%u", &p[i]));
    }
    return;
}
#undef OSIZE
//Резервное сохранение в файл
static void do_backup(int step/*, double sumtime*/){
    char filename[50];
    //Имя текущего файла
    sprintf(filename, "backup%06d_%06d_%10.4fK.txt",step,PARTICLENUMBER,PREF_TEMP);
    FILE *file_backup = fopen(filename, "w");
    //Пишем данные в файл
    myfwrite(&step, sizeof(step),1,file_backup);
    //myfwrite(&sumtime, sizeof(sumtime),1,file_backup);
    myfwrite(Cx,sizeof(Cx[0]),PARTICLENUMBER,file_backup);
    myfwrite(Cy,sizeof(Cy[0]),PARTICLENUMBER,file_backup);
    myfwrite(Cz,sizeof(Cz[0]),PARTICLENUMBER,file_backup);

    myfwrite(Vx,sizeof(Vx[0]),PARTICLENUMBER,file_backup);
    myfwrite(Vy,sizeof(Vy[0]),PARTICLENUMBER,file_backup);
    myfwrite(Vz,sizeof(Vz[0]),PARTICLENUMBER,file_backup);
    fclose(file_backup);
    //Удаление старых файлов при необходимости
    if(step > 2*BACKUP_FREQ){
        sprintf(filename, "backup%06d_%06d_%10.4fK.txt", step-2*BACKUP_FREQ,PARTICLENUMBER,PREF_TEMP);
        remove(filename);
    }
    return;
}
//Восстановление из резервной копии
static int restore_backup(int step/*, double *sumtime*/){
    //Сохраняем только на итерациях кратных BACKUP_FREQ
    assert(step % BACKUP_FREQ==0);
    char filename[50];
    int step2;
    //Имя файла
    sprintf(filename, "backup%06d_%06d_%10.4fK.txt", step, PARTICLENUMBER,PREF_TEMP);
    FILE *file_backup = fopen(filename, "r");
    if(!file_backup){
        fprintf(stderr, "Error: no restore file (%s)\n", filename);
        return 1;
    }
    //Считываем данные из файла
    myfread(&step2, sizeof(step2),1,file_backup);
    assert(step2==step);
    //myfread(sumtime, sizeof(sumtime),1,file_backup);
    myfread(Cx, sizeof(Cx[0]),PARTICLENUMBER,file_backup);
    myfread(Cy, sizeof(Cy[0]),PARTICLENUMBER,file_backup);
    myfread(Cz, sizeof(Cz[0]),PARTICLENUMBER,file_backup);
    
    myfread(Vx, sizeof(Vx[0]),PARTICLENUMBER,file_backup);
    myfread(Vy, sizeof(Vy[0]),PARTICLENUMBER,file_backup);
    myfread(Vz, sizeof(Vz[0]),PARTICLENUMBER,file_backup);

    fclose(file_backup);
    return 0;
}
