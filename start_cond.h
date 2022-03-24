void startingCoordsPrimCube(){
    double x=0,y=0,z=0;
    if(!PGU){
        for(int i=0;i<PARTICLENUMBER;i++){
            if(z>NUMKRIST_X){
                y++;
                z=0;
            }
            if(y>NUMKRIST_X){
                x++;
                z=0;
            }
            Cx[i]=x*REBROKR;
            Cy[i]=y*REBROKR;
            Cz[i]=z*REBROKR;
            z++;
        }
    }else{
        for(int i=0;i<PARTICLENUMBER;i++){
            if(z>NUMKRIST_X-1){
                y++;
                z=0;
            }
            if(y>NUMKRIST_X-1){
                x++;
                y=0;
            }
            Cx[i]=x*REBROKR;
            Cy[i]=y*REBROKR;
            Cz[i]=z*REBROKR;
            z++;
        }
    }
}

void startPrimCube(){
    startingCoordsPrimCube();
    speedGeneration();
}

void start_cond(int type,int pgu){
    switch(type){
        case 0: //2 молекулы
            break;
        case 1: //Примитивная кубическая решетка
            startPrimCube();
            break;
        case 2: //Гранецентрированная кубическая решетка    
            break;
    }
}
