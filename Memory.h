#define AllocMem(arr,n,type) arr=(type*) malloc((n)*sizeof(type))

#define checkAlloc(arr) arr==NULL ? std::cout<<"Error memalloc!\n":std::cout<<"Memalloc successful!\n"

void AllocArrays(){
    //Выделение памяти для координат всех частиц
    AllocMem(Cx,PARTICLENUMBER,double);
    checkAlloc(Cx);
    AllocMem(Cy,PARTICLENUMBER,double);
    checkAlloc(Cy);
    AllocMem(Cz,PARTICLENUMBER,double);
    checkAlloc(Cz);
    //Выделение памяти для скоростей всех частиц
    AllocMem(Vx,PARTICLENUMBER,double);
    checkAlloc(Vx);
    AllocMem(Vy,PARTICLENUMBER,double);
    checkAlloc(Vy);
    AllocMem(Vz,PARTICLENUMBER,double);
    checkAlloc(Vz);
    //Выделение памяти для сил действующих на частицы
    AllocMem(Fx,PARTICLENUMBER,double);
    checkAlloc(Fx);
    AllocMem(Fy,PARTICLENUMBER,double);
    checkAlloc(Fy);
    AllocMem(Fz,PARTICLENUMBER,double);
    checkAlloc(Fz);
    //Выделение памяти для новых сил действующих на частицы
    AllocMem(Fxnew,PARTICLENUMBER,double);
    checkAlloc(Fxnew);
    AllocMem(Fynew,PARTICLENUMBER,double);
    checkAlloc(Fynew);
    AllocMem(Fznew,PARTICLENUMBER,double);
    checkAlloc(Fznew);
    
    AllocMem(V0x,PARTICLENUMBER,double);
    checkAlloc(V0x);
    AllocMem(V0y,PARTICLENUMBER,double);
    checkAlloc(V0y);
    AllocMem(V0z,PARTICLENUMBER,double);
    checkAlloc(V0z);
}