#define _USE_MATH_DEFINES

#include "kernel11.h"

__constant__ double SENSOR2;
__constant__ int NAHO;

__device__ int selectedCounts[MACRO_NMAX];
__device__ double tmpPhero_d[MACRO_MAX][MACRO_MAX];
__device__ curandState rnd_state[MACRO_NMAX];

//Misc
__device__ bool isGotFood(Food& food);
__device__ double atomicAddDouble(double* address, double val);
__device__ enum Direction genDirRand(int id);
__device__ double genProbRand(int id);
__device__ int genAntNumRand(int id);
__device__ double degToRad(double a);
__device__ double dist(Cell a,Cell b);
__device__ double distCandP(Cell a,double x,double y);
__device__ bool isOppositeDir(enum Direction nestDir,enum Direction dir);
__device__ bool isOppositeDir(Cell& cell, enum Direction dir);
__device__ enum Direction selectNextDir(Cell& cell, enum Direction dir);
__device__ double hilFunc(double x,double alpha);

//Initializer
__host__ void getDevicePtrs();
__global__ void randInit();
__global__ void antsInit();
__global__ void cellsInit();
__global__ void setNest();
__global__ void setDistFromNest();
__global__ void setNestDirs();
__global__ void setNearestDirFromNest();
__global__ void setFoodsDir();

//Calculation functions
__global__ void selectAnts();
__global__ void naturalFoodDecrease();
__global__ void evapolation();
__global__ void chemotaxis();
__global__ void diffusion();
__global__ void pheroUpdate();
__global__ void gauss_path();


__host__ void calculation(){
    naturalFoodDecrease<<<1,MACRO_NUM_FOODS>>>();
    evapolation<<<MACRO_MAX,MACRO_MAX>>>();

    //sortKeyInit<<<1,MACRO_NMAX>>>();
    //thrust::sort_by_key(sort_key_d_ptr, sort_key_d_ptr + MACRO_NMAX, ants_d_ptr);

    selectAnts<<<1,MACRO_NMAX>>>();
    chemotaxis<<<1,MACRO_NMAX>>>();
    //cudaMemcpyFromSymbol(cells,cells_d,MACRO_MAX*MACRO_MAX*sizeof(Cell),0);
    //chemotaxis();
    //cudaMemcpyToSymbol(cells_d,cells,MACRO_MAX*MACRO_MAX*sizeof(Cell),0);
    diffusion<<<MACRO_MAX,MACRO_MAX>>>();
    pheroUpdate<<<MACRO_MAX,MACRO_MAX>>>();
}

//Initialize

__host__ void getDevicePtrs(){
    cudaGetSymbolAddress((void**)&sort_key_d_ptr_raw, sort_key_d);
    sort_key_d_ptr = thrust::device_ptr<unsigned int>(sort_key_d_ptr_raw);

    cudaGetSymbolAddress((void**)&seeds_d_ptr_raw, seeds_d);
    seeds_d_ptr = thrust::device_ptr<unsigned long long int>(seeds_d_ptr_raw);

    cudaGetSymbolAddress((void**)&ants_d_ptr_raw, ants_d);
    ants_d_ptr = thrust::device_ptr<Ant>(ants_d_ptr_raw);

    cudaGetSymbolAddress((void**)&cells_d_ptr_raw, cells_d);
    cells_d_ptr = thrust::device_ptr<Cell>(cells_d_ptr_raw);

    cudaGetSymbolAddress((void**)&gauss_path1_d_ptr_raw, gauss_path1_d);
    gauss_path1_d_ptr = thrust::device_ptr<Cell>(gauss_path1_d_ptr_raw);

    cudaGetSymbolAddress((void**)&gauss_path2_d_ptr_raw, gauss_path2_d);
    gauss_path2_d_ptr = thrust::device_ptr<Cell>(gauss_path2_d_ptr_raw);

    cudaGetSymbolAddress((void**)&foods_d_ptr_raw, foods_d);
    foods_d_ptr = thrust::device_ptr<Food>(foods_d_ptr_raw);
}

__global__ void gauss_path(){
    const int i = threadIdx.x;
    const int j = blockIdx.x;

    double vx = cells_d[i][j].cart.x - cells_d[MACRO_NEST_Y][MACRO_NEST_X].cart.x;
    double vy = cells_d[i][j].cart.y - cells_d[MACRO_NEST_Y][MACRO_NEST_X].cart.y;

    double inner[2] = {vx * x_Nest2Food_d[0] + vy * y_Nest2Food_d[0], vx * x_Nest2Food_d[1] + vy * y_Nest2Food_d[1]};
    double px[2] = {vx - inner[0] * x_Nest2Food_d[0], vx - inner[1] * x_Nest2Food_d[1]};
    double py[2] = {vy - inner[0] * y_Nest2Food_d[0], vy - inner[1] * y_Nest2Food_d[1]};

    double dij[2] = {px[0]*px[0]+py[0]*py[0], px[1]*px[1]+py[1]*py[1]};

    gauss_path1_d[i][j].phero = exp(-dij[0]/(4.0*MACRO_DIFFE))/(2.0*sqrt(M_PI*MACRO_DIFFE));
    gauss_path2_d[i][j].phero = exp(-dij[1]/(4.0*MACRO_DIFFE))/(2.0*sqrt(M_PI*MACRO_DIFFE));
}

__global__ void randInit(){
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    curand_init(seeds_d[id],0,0,&rnd_state[id]);
}

__global__ void antsReset(){
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    ants_d[id].status = FORAGE;
    ants_d[id].i = MACRO_NEST_Y;
    ants_d[id].j = MACRO_NEST_X;
    ants_d[id].searchTime = 0;
    ants_d[id].dir = genDirRand(id);
    for (int i=0; i<MACRO_NUM_FOODS; i++){
        ants_d[id].homing[i] = 0;
    }
    if(id<NAHO){
        ants_d[id].ch = FOOL_CH;
    }
    else {
        ants_d[id].ch = NORMAL_CH;
    }
}

__global__ void cellsReset(){
    const int i = threadIdx.x;
    const int j = blockIdx.x;
    cells_d[i][j].phero = 0.0;
}

__global__ void cellsInit(){
    const int i = threadIdx.x;
    const int j = blockIdx.x;

    cells_d[i][j].foodNo = -1;
    cells_d[i][j].status = NORMAL_CELL;

    //Cell number initialize
    cells_d[i][j].i = i;
    cells_d[i][j].j = j;

    //Cartesian initialize
    cells_d[i][j].cart.x = (j-MACRO_CART_X_ZERO)*(sqrt(3.0)/2.0);
    cells_d[i][j].cart.y = (abs(j-MACRO_CART_X_ZERO)%2)/2.0+(i-MACRO_CART_Y_ZERO);
    //Edge initialize
    cells_d[i][j].edge = NONE;

    //Nest Dir initialize
    cells_d[i][j].nestDir = NONE;

    cells_d[i][j].distFromNest = 0.0;

    // Gauss path for food source 1
    gauss_path1_d[i][j].foodNo = -1;
    gauss_path1_d[i][j].status = NORMAL_CELL;

    //Cell number initialize
    gauss_path1_d[i][j].i = i;
    gauss_path1_d[i][j].j = j;

    //Cartesian initialize
    gauss_path1_d[i][j].cart.x = cells_d[i][j].cart.x;
    gauss_path1_d[i][j].cart.y = cells_d[i][j].cart.y;
    //Edge initialize
    gauss_path1_d[i][j].edge = NONE;

    //Nest Dir initialize
    gauss_path1_d[i][j].nestDir = NONE;

    gauss_path1_d[i][j].distFromNest = 0.0;

    // Gauss path for food source 1
    gauss_path2_d[i][j].foodNo = -1;
    gauss_path2_d[i][j].status = NORMAL_CELL;

    //Cell number initialize
    gauss_path2_d[i][j].i = i;
    gauss_path2_d[i][j].j = j;

    //Cartesian initialize
    gauss_path2_d[i][j].cart.x = cells_d[i][j].cart.x;
    gauss_path2_d[i][j].cart.y = cells_d[i][j].cart.y;
    //Edge initialize
    gauss_path2_d[i][j].edge = NONE;

    //Nest Dir initialize
    gauss_path2_d[i][j].nestDir = NONE;

    gauss_path2_d[i][j].distFromNest = 0.0;
}


__global__ void setEdges(){
    const int i = threadIdx.x;
    const int j = blockIdx.x;
    if(i==MACRO_MAX-1){ //For upper edge
        cells_d[i][j].edge |= UP;
        if(abs((j-MACRO_CART_X_ZERO)%2)==1){
            cells_d[i][j].edge |= (UPLEFT | UPRIGHT);
        }
    }
    else if(i==0){//For lower edge
        cells_d[i][j].edge |= LOW;
        if(abs((j-MACRO_CART_X_ZERO)%2)==0){
            cells_d[i][j].edge |= LOWLEFT | LOWRIGHT;
        }
    }

    if(j==0){//For left edge
        cells_d[i][j].edge |= LEFT;
    }
    else if(j==MACRO_MAX-1){//For right edge
        cells_d[i][j].edge |= RIGHT;
    }
    gauss_path1_d[i][j].edge = cells_d[i][j].edge;
    gauss_path2_d[i][j].edge = cells_d[i][j].edge;
}

__global__ void setNest(){
    const int i = threadIdx.x;
    const int j = blockIdx.x;

    Cell* c;
    if(i==MACRO_NEST_Y && j==MACRO_NEST_X){
        cells_d[MACRO_NEST_Y][MACRO_NEST_X].status |= NEST_CELL;

        gauss_path1_d[MACRO_NEST_Y][MACRO_NEST_X].edge = cells_d[MACRO_NEST_Y][MACRO_NEST_X].edge;
        gauss_path2_d[MACRO_NEST_Y][MACRO_NEST_X].edge = cells_d[MACRO_NEST_Y][MACRO_NEST_X].edge;

        for(enum Direction d = UP; d<=UPLEFT; (d<<=1) ){
            c = getCell(cells_d,MACRO_NEST_Y,MACRO_NEST_X,d);
            c->status |= NEST_NEIGHBOUR_CELL;

            c = getCell(gauss_path1_d,MACRO_NEST_Y,MACRO_NEST_X,d);
            c->status |= NEST_NEIGHBOUR_CELL;

            c = getCell(gauss_path2_d,MACRO_NEST_Y,MACRO_NEST_X,d);
            c->status |= NEST_NEIGHBOUR_CELL;
        }
    }
}

__global__ void setDistFromNest(){
    const int i = threadIdx.x;
    const int j = blockIdx.x;

    Cell *nest_c;
    nest_c = &cells_d[MACRO_NEST_Y][MACRO_NEST_X];
    double d = dist(cells_d[i][j],*nest_c);
    cells_d[i][j].distFromNest = d;
    gauss_path1_d[i][j].distFromNest = d;
    gauss_path2_d[i][j].distFromNest = d;
}

__device__ double dot(Cartesian a, Cartesian b) {
    return (a.x * b.x + a.y * b.y);
}

__device__ double cross(Cartesian a, Cartesian b) {
    return (a.x * b.y - a.y * b.x);
}

__global__ void setCriticalAngle() {
    const int i = threadIdx.x;
    const int j = blockIdx.x;


    cells_d[i][j].criticalAngle = NONE;
    gauss_path1_d[i][j].criticalAngle = NONE;
    gauss_path2_d[i][j].criticalAngle = NONE;

    if( (cells_d[i][j].status&NEAR_NEST)!=NORMAL_CELL ){
        return;
    }

    Cartesian c = cells_d[i][j].cart;
    c.x = -c.x/cells_d[i][j].distFromNest;
    c.y = -c.y/cells_d[i][j].distFromNest;

    for(enum Direction dir = UP; dir<=UPLEFT; (dir<<=1) ) {
        Cartesian d;

        switch (dir) {
            case UP:
                d.x = 0;
                d.y = 1;
                break;
            case UPRIGHT:
                d.x = 1;
                d.y = tan(M_PI/4.0);
                break;
            case LOWRIGHT:
                d.x = 1;
                d.y = -tan(M_PI/4.0);
                break;
            case LOW:
                d.x = 0;
                d.y = -1;
                break;
            case LOWLEFT:
                d.x = -1;
                d.y = -tan(M_PI/4.0);
                break;
            case UPLEFT:
                d.x = -1;
                d.y = tan(M_PI/4.0);
                break;
            default:
                break;
        }

        d.x = d.x/sqrt(dot(d,d));
        d.y = d.y/sqrt(dot(d,d));

        double dotVal = dot(c,d);
        if (dotVal<=0.3){
            cells_d[i][j].criticalAngle |= dir;
            gauss_path1_d[i][j].criticalAngle |= dir;
            gauss_path2_d[i][j].criticalAngle |= dir;
        }
    }
}


__global__ void setNearestDirFromNest(){
    const int i = threadIdx.x;
    const int j = blockIdx.x;
    Cell& c = cells_d[i][j];
    Cell& gp1 = gauss_path1_d[i][j];
    Cell& gp2 = gauss_path2_d[i][j];

    for (int itr=0; itr<6; itr++){
        c.nearestDirFromNestList[itr] = NONE;
        gp1.nearestDirFromNestList[itr] = NONE;
        gp2.nearestDirFromNestList[itr] = NONE;
    }

    enum Direction dir = UP;
    for(int itr=0; dir<=UPLEFT; itr++) {
        if ( c.criticalAngle&dir == NONE ){
            continue;
        }

        c.nearestDirFromNestList[itr] = selectNextDir(c, dir);
        gp1.nearestDirFromNestList[itr] = c.nearestDirFromNestList[itr];
        gp2.nearestDirFromNestList[itr] = c.nearestDirFromNestList[itr];
        dir<<=1;
    }
}

__global__ void setNestDirs(){
    const int i = threadIdx.x;
    const int j = blockIdx.x;

    Cell *c;

    double d = cells_d[i][j].distFromNest;
    double tmp;
    for(enum Direction dir = UP; dir<=UPLEFT; (dir<<=1) ){

        c = getCell(cells_d,i,j,dir);


        tmp=c->distFromNest;
        if( fabs(tmp-d)<MACRO_EPS ){
            cells_d[i][j].nestDir |= dir;
        }
        else if(tmp<d) {
            cells_d[i][j].nestDir = dir;
            d = tmp;
        }
        gauss_path1_d[i][j].nestDir = cells_d[i][j].nestDir;
        gauss_path2_d[i][j].nestDir = cells_d[i][j].nestDir;
    }
}

__global__ void foodsReset(){
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    foods_d[i].vol = MACRO_FOODSOURCE;
}

__global__ void setFoodsDir(){
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    const double dtheta = degToRad(MACRO_FOOD_ANGLE);

    Cell *nearCell=NULL;
    Cell *nearCellgp1=NULL;
    Cell *nearCellgp2=NULL;

    double x,y;
    x=MACRO_FOOD_DIST * cos(i*dtheta);
    y=MACRO_FOOD_DIST * sin(i*dtheta);
    for(int j=0; j<MACRO_MAX; j++){
        for(int k=0; k<MACRO_MAX; k++){
            if(distCandP(cells_d[j][k],x,y)<=sqrt(3.0)/3.0+MACRO_EPS){
                nearCell = &cells_d[j][k];
                nearCellgp1 = &gauss_path1_d[j][k];
                nearCellgp2 = &gauss_path2_d[j][k];
                break;
            }
        }
    }
    if(nearCell==NULL){
    }
    else{
        Cell *c=NULL;
        Cell *gp1=NULL;
        Cell *gp2=NULL;
        double d = distCandP(*nearCell,x,y);
        int j = nearCell->i;
        int k = nearCell->j;

        for(enum Direction dir = UP; dir<=UPLEFT; (dir<<=1) ){
            c = getCell(cells_d,j,k,dir);
            gp1 = getCell(gauss_path1_d,j,k,dir);
            gp2 = getCell(gauss_path2_d,j,k,dir);
            if( distCandP(*c,x,y)<d ){
                nearCell = c;
                nearCellgp1 = gp1;
                nearCellgp2 = gp2;
                d = distCandP(*nearCell,x,y);
            }
        }
        foods_d[i].i = nearCell->i;
        foods_d[i].j = nearCell->j;

        nearCell->foodNo = i;
        nearCell->status |= FOOD_CELL;
        nearCellgp1->foodNo = nearCell->foodNo;
        nearCellgp1->status = nearCell->status;
        nearCellgp2->foodNo = nearCell->foodNo;
        nearCellgp2->status = nearCell->status;


        for(enum Direction dir = UP; dir<=UPLEFT; (dir<<=1) ){
            c = getCell(cells_d,foods_d[i].i,foods_d[i].j,dir);
            gp1 = getCell(gauss_path1_d,foods_d[i].i,foods_d[i].j,dir);
            gp2 = getCell(gauss_path2_d,foods_d[i].i,foods_d[i].j,dir);
            c->foodNo = i;
            c->status |= FOOD_NEIGHBOUR_CELL;
            gp1->foodNo = c->foodNo;
            gp1->status = c->status;
            gp2->foodNo = c->foodNo;
            gp2->status = c->status;
        }
    }

}


//Calculation

__global__ void selectAnts(){
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    int rnd =  genAntNumRand(id);
    atomicAdd(&selectedCounts[rnd], 1);
}

__global__ void sortKeyInit(){
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    sort_key_d[id] = curand(&rnd_state[id]);
    //printf("id:%d,%u\n",id,sort_key_d[id]);
}

__global__ void diffusion(){
    const int i = blockIdx.x;
    const int j = threadIdx.x;

    double tmp = 0.0;
    for (enum Direction dir = UP; dir<=UPLEFT; (dir<<=1) ){
        tmp += getCell(cells_d,i,j,dir)->phero;
    }
    tmpPhero_d[i][j] = cells_d[i][j].phero+MACRO_DIFFE*(tmp/6.0-cells_d[i][j].phero);
}

__global__ void pheroUpdate(){
    const int i = blockIdx.x;
    const int j = threadIdx.x;

    cells_d[i][j].phero = tmpPhero_d[i][j];
}

__global__ void naturalFoodDecrease(){
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    // 常に一定量で増加する。
    foods_d[id].vol = foods_d[id].vol + MACRO_REC;
    if (foods_d[id].vol > 100.0)
        foods_d[id].vol = 100.0;
}

__global__ void evapolation(){
    const int i = blockIdx.x;
    const int j = threadIdx.x;
    cells_d[i][j].phero *= (1.0-MACRO_EVAPOLATION_CONST);
}


__global__ void chemotaxis(){
    const int id = threadIdx.x + blockIdx.x * blockDim.x;
    Ant *ant = &(ants_d[id]);

    for(int dummy=0; dummy<selectedCounts[id]; dummy++){
        ant->searchTime++;

        int i = ant->i;
        int j = ant->j;
        enum Direction dir = ant->dir;
        enum Direction nestDir = cells_d[i][j].nestDir;

        double leftPhero, frontPhero, rightPhero;

        Cell *leftCell  = getCell(cells_d,i,j,left(dir));
        Cell *frontCell = getCell(cells_d,i,j,dir);
        Cell *rightCell = getCell(cells_d,i,j,right(dir));

        if(
                ant->searchTime>=MACRO_MAX_SEARCH_TIME
                && ant->status!=EMERGENCY
          ){
            ant->status = EMERGENCY;
        }

        if(ant->status==GOHOME){
            atomicAddDouble(&(cells_d[i][j].phero),MACRO_EMI*MACRO_ENEST);
        }
        __threadfence();
        if(ant->status==RANDOM_SEARCH){
            leftPhero = 1.0;
            frontPhero = 1.0;
            rightPhero = 1.0;
        }
        else {
            leftPhero = leftCell->phero;
            frontPhero = frontCell->phero;
            rightPhero = rightCell->phero;
        }

        if( (ant->status==GOHOME || ant->status==EMERGENCY) && isOppositeDir(cells_d[i][j], dir)){

            enum Direction nextDir = cells_d[i][j].nearestDirFromNestList[dirToNum(dir)];

            if( nextDir == left(dir) ){
                ant->dir = left(dir);
                frontCell = leftCell;
            }
            else if( nextDir == right(dir) ){
                ant->dir = right(dir);
                frontCell = rightCell;
            }
            else{
                if(genProbRand(id)<=0.5){
                    ant->dir = right(dir);
                    frontCell = rightCell;
                }
                else{
                    ant->dir = left(dir);
                    frontCell = leftCell;
                }
            }
            ant->i = frontCell->i;
            ant->j = frontCell->j;
        }
        else{
            double s1,s2,s3,s12,t,tot,rand;
            if(ant->ch == NORMAL_CH){
                t = MACRO_HIL_CONST;
            }
            else{
                t = SENSOR2*MACRO_HIL_CONST;
            }

            s1=hilFunc(leftPhero,t);
            s2=hilFunc(frontPhero,t);
            s3=hilFunc(rightPhero,t);
            /*
               if(s1<MACRO_EPS && s2<MACRO_EPS && s3<MACRO_EPS){
               s1=1.0;
               s2=1.0;
               s3=1.0;
               }
               */
            tot = s1+s2+s3;
            s1/=tot;
            s2/=tot;

            s12=s1+s2;

            rand=genProbRand(id);

            if(rand<=s1){
                ant->dir = left(dir);
                ant->i   = leftCell->i;
                ant->j   = leftCell->j;
            }
            else if(rand<=s12){
                ant->i   = frontCell->i;
                ant->j   = frontCell->j;
            }
            else{
                ant->dir = right(dir);
                ant->i   = rightCell->i;
                ant->j   = rightCell->j;
            }

        }

        if( (cells_d[ant->i][ant->j].status&NEAR_FOOD)!=NORMAL_CELL
                &&  foods_d[  cells_d[ant->i][ant->j].foodNo  ].vol>=0.1
                &&  (ant->status != GOHOME && ant->status != EMERGENCY) ){
            //atomicAddDouble(&(foods_d[  cells_d[ant->i][ant->j].foodNo  ].vol),-MACRO_UNIT);
            //ant->status = GOHOME;
            //ant->searchTime = 0;
            int fNo = cells_d[ant->i][ant->j].foodNo;

            if(isGotFood(foods_d[fNo])){
                ant->status = GOHOME;
                ant->searchTime = 0;
                ant->_foodNo = fNo;
                ant->dir = left(left(left(dir)));
            }
        }
        __threadfence();

        if( (cells_d[ant->i][ant->j].status&NEAR_NEST)!=NORMAL_CELL
                &&  (ant->status == GOHOME || ant->status == EMERGENCY)){
            if(ant->status == GOHOME){
                ant->homing[ant->_foodNo]++;
                //atomicAddDouble(&(cells_d[i][j].phero),MACRO_EMI*MACRO_ENEST);
            }
            ant->status = FORAGE;
            ant->searchTime = 0;
            ant->dir = genDirRand(id);
            ant->i   = MACRO_NEST_Y;
            ant->j   = MACRO_NEST_X;
        }
    }
    selectedCounts[id] = 0;
}


//DataHandler
__device__ __host__  enum Direction operator<<(enum Direction d, int i){
    return static_cast<enum Direction>(static_cast<int>(d)<<i);
}

__device__ __host__  enum Direction operator>>(enum Direction d, int i){
    return static_cast<enum Direction>(static_cast<int>(d)>>i);
}

__device__ __host__  enum Direction operator|(enum Direction d1, enum Direction d2){
    return static_cast<enum Direction>(static_cast<int>(d1)|static_cast<int>(d2));
}
__device__ __host__  enum Direction operator&(enum Direction d1, enum Direction d2){
    return static_cast<enum Direction>(static_cast<int>(d1)&static_cast<int>(d2));
}

__device__ __host__  enum Direction& operator|=(enum Direction& d1, enum Direction d2){
    d1 = (d1 | d2);
    return d1;
}

__device__ __host__  enum Direction& operator&=(enum Direction& d1, enum Direction d2){
    d1 = (d1 & d2);
    return d1;
}

__device__ __host__  enum Direction& operator<<=(enum Direction& d1, int i){
    d1 = (d1 << i);
    return d1;
}

__device__ __host__  enum Direction& operator>>=(enum Direction& d1, int i){
    d1 = (d1 >> i);
    return d1;
}

__device__ __host__  bool operator<=(enum Direction d1, enum Direction d2){
    return (static_cast<int>(d1) <= static_cast<int>(d2));
}







__device__ __host__  enum CELLStatus operator<<(enum CELLStatus d, int i){
    return static_cast<enum CELLStatus>(static_cast<int>(d)<<i);
}

__device__ __host__  enum CELLStatus operator>>(enum CELLStatus d, int i){
    return static_cast<enum CELLStatus>(static_cast<int>(d)>>i);
}

__device__ __host__  enum CELLStatus operator|(enum CELLStatus d1, enum CELLStatus d2){
    return static_cast<enum CELLStatus>(static_cast<int>(d1)|static_cast<int>(d2));
}
__device__ __host__  enum CELLStatus operator&(enum CELLStatus d1, enum CELLStatus d2){
    return static_cast<enum CELLStatus>(static_cast<int>(d1)&static_cast<int>(d2));
}

__device__ __host__  enum CELLStatus& operator|=(enum CELLStatus& d1, enum CELLStatus d2){
    d1 = (d1 | d2);
    return d1;
}

__device__ __host__  enum CELLStatus& operator&=(enum CELLStatus& d1, enum CELLStatus d2){
    d1 = (d1 & d2);
    return d1;
}



__device__ __host__ __forceinline__ enum Direction left(enum Direction dir){
    if(dir == UP){
        return UPLEFT;
    }
    else{
        return (dir >> 1)&ALL_DIR;
    }
}

__device__ __host__ __forceinline__ enum Direction right(enum Direction dir){
    if(dir == UPLEFT){
        return UP;
    }
    else{
        return (dir << 1)&ALL_DIR;
    }
}

__device__ __host__ __forceinline__ Cell* up(Cell cells[MACRO_MAX][MACRO_MAX],int i,int j){
    if( (cells[i][j].edge&UP)!=NONE ){
        return &cells[0][j];
    }
    else{
        return &cells[i+1][j];
    }
}

__device__ __host__ __forceinline__ Cell* upright(Cell cells[MACRO_MAX][MACRO_MAX],int i,int j){
    int ii,jj;
    if( (cells[i][j].edge&UPRIGHT)!=NONE ){
        jj = 0;
        if(abs(j-MACRO_CART_X_ZERO)%2==0){
            ii = i;
        }
        else{
            ii = i+1;
            if(ii==MACRO_MAX){
                ii = 0;
            }
        }
    }
    else{
        jj = j+1;
        if(abs(j-MACRO_CART_X_ZERO)%2==0){
            ii = i;
        }
        else{
            ii = i+1;
        }
    }
    return &cells[ii][jj];
}

__device__ __host__ __forceinline__ Cell* lowright(Cell cells[MACRO_MAX][MACRO_MAX],int i,int j){

    int ii,jj;

    if( (cells[i][j].edge&LOWRIGHT)!=NONE ){
        jj = 0;
        if(abs(j-MACRO_CART_X_ZERO)%2==0){
            ii = i-1;
            if(ii<0){
                ii=MACRO_MAX-1;
            }
        }
        else{
            ii = i;
        }
    }
    else{
        jj = j+1;
        if(abs(j-MACRO_CART_X_ZERO)%2==0){
            ii = i-1;
        }
        else{
            ii = i;
        }
    }
    return &cells[ii][jj];
}

__device__ __host__ __forceinline__ Cell* low(Cell cells[MACRO_MAX][MACRO_MAX],int i,int j){
    if( (cells[i][j].edge&LOW)!=NONE ){
        return &cells[MACRO_MAX-1][j];
    }
    else{
        return &cells[i-1][j];
    }
}

__device__ __host__ __forceinline__ Cell* lowleft(Cell cells[MACRO_MAX][MACRO_MAX],int i,int j){
    int ii,jj;

    if( (cells[i][j].edge&LOWLEFT)!=NONE ){
        jj = MACRO_MAX-1;
        if(abs(j-MACRO_CART_X_ZERO)%2==0){
            ii = i-1;
            if(ii<0){
                ii = MACRO_MAX-1;
            }
        }
        else{
            ii = i;
        }
    }
    else{
        jj = j-1;
        if(abs(j-MACRO_CART_X_ZERO)%2==0){
            ii = i-1;
        }
        else{
            ii=i;
        }
    }
    return &cells[ii][jj];
}

__device__ __host__ __forceinline__ Cell* upleft(Cell cells[MACRO_MAX][MACRO_MAX],int i,int j){
    int ii,jj;
    if( (cells[i][j].edge&UPLEFT)!=NONE ){
        jj = MACRO_MAX-1;

        if(abs(j-MACRO_CART_X_ZERO)%2==0){
            ii = i;
        }
        else{
            ii= i+1;
            if(ii==MACRO_MAX){
                ii=0;
            }
        }
    }
    else{
        jj = j-1;
        if(abs(j-MACRO_CART_X_ZERO)%2==0){
            ii = i;
        }
        else{
            ii = i+1;
        }
    }
    return &cells[ii][jj];
}

__device__ __host__ Cell* getCell(Cell cells[MACRO_MAX][MACRO_MAX],int i,int j, enum Direction dir){

    switch (dir){
        case UP:
            return up(cells,i,j);
        case UPRIGHT:
            return upright(cells,i,j);
        case LOWRIGHT:
            return lowright(cells,i,j);
        case LOW:
            return low(cells,i,j);
        case LOWLEFT:
            return lowleft(cells,i,j);
        case UPLEFT:
            return upleft(cells,i,j);
        default:
            return NULL;
    }
}

__device__ __host__ int dirToNum(enum Direction dir){
    switch (dir){
        case UP:
            return 0;
        case UPRIGHT:
            return 1;
        case LOWRIGHT:
            return 2;
        case LOW:
            return 3;
        case LOWLEFT:
            return 4;
        case UPLEFT:
            return 5;
        default:
            return -1;
    }
}



//Misc
__device__ __forceinline__ bool isGotFood(Food& food){
    unsigned long long int* address_as_ull =
        (unsigned long long int*)(&(food.vol));
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        if(__longlong_as_double(assumed)<0.1){
            return false;
        }
        old = atomicCAS(address_as_ull, assumed,__double_as_longlong(-MACRO_UNIT + __longlong_as_double(assumed)));
    } while (assumed != old);
    return true;
}

__device__ __forceinline__ double atomicAddDouble(double* address, double val){
    unsigned long long int* address_as_ull =
        (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                __double_as_longlong(val +
                    __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__device__ __forceinline__ enum Direction genDirRand(int id){
    return static_cast<enum Direction>(1 << (curand(&rnd_state[id])%6));
}

__device__ __forceinline__ double genProbRand(int id){
    return curand_uniform_double(&rnd_state[id]);
}

__device__ __forceinline__ int genAntNumRand(int id){
    return curand(&rnd_state[id])%MACRO_NMAX;
}

__device__ __forceinline__ double degToRad(double a) {
    return a * M_PI / 180.0;
}

__device__ __forceinline__ double dist(Cell a,Cell b){
    return sqrt( (a.cart.x - b.cart.x)*(a.cart.x - b.cart.x)
            + (a.cart.y - b.cart.y)*(a.cart.y - b.cart.y) );
}

__device__ __forceinline__ double distCandP(Cell a,double x,double y){
    return sqrt( (a.cart.x - x)*(a.cart.x - x)
            + (a.cart.y - y)*(a.cart.y - y) );
}

__device__ __forceinline__ bool isOppositeDir(enum Direction nestDir,enum Direction dir){
    //If theta = 60 deg., this is OK.
    if( (dir&nestDir)        !=NONE
            ||  (left(dir)&nestDir)  !=NONE
            ||  (right(dir)&nestDir) !=NONE){
        return false;
    }
    else{
        return true;
    }
}

__device__ __forceinline__ bool isOppositeDir(Cell& cell, enum Direction dir){
    if ( (cell.criticalAngle & dir)==dir ){
        return true;
    }
    else{
        return false;
    }
}

__device__ __forceinline__ enum Direction selectNextDir(Cell& cell, enum Direction dir){
    int rightCount = 0;
    int leftCount  = 0;
    for (enum Direction currentDir=right(dir); currentDir!=dir; currentDir=right(currentDir)){
        if( (cell.criticalAngle & currentDir)!=currentDir ){
            break;
        }
        rightCount++;
    }

    for (enum Direction currentDir=left(dir); currentDir!=dir; currentDir=left(currentDir)){
        if( (cell.criticalAngle & currentDir)!=currentDir ){
            break;
        }
        leftCount++;
    }

    if ( rightCount < leftCount ){
        return right(dir);
    }
    else if ( rightCount > leftCount ){
        return left(dir);
    }
    else{
        return NONE;
    }
}

__device__ __forceinline__ double hilFunc(double x,double alpha){
    return pow(alpha*x+0.05,10);
}

__host__ void initialize(){
    getDevicePtrs();

    //antsInit<<<MACRO_NMAX,1>>>();
    cellsInit<<<MACRO_MAX,MACRO_MAX>>>();

    setEdges<<<MACRO_MAX,MACRO_MAX>>>();
    setNest<<<MACRO_MAX,MACRO_MAX>>>();
    setDistFromNest<<<MACRO_MAX,MACRO_MAX>>>();

    setCriticalAngle<<<MACRO_MAX,MACRO_MAX>>>();
    setNearestDirFromNest<<<MACRO_MAX,MACRO_MAX>>>();

    setNestDirs<<<MACRO_MAX,MACRO_MAX>>>();
    setFoodsDir<<<MACRO_NUM_FOODS,1>>>();

    // addition by Shiraishi
    for (int i=0; i<MACRO_NUM_FOODS; i++){
        x_Nest2Food_d[i] = cells_d[foods_d[i].i][foods_d[i].j].cart.x - cells_d[MACRO_NEST_Y][MACRO_NEST_X].cart.x;
        y_Nest2Food_d[i] = cells_d[foods_d[i].i][foods_d[i].j].cart.y - cells_d[MACRO_NEST_Y][MACRO_NEST_X].cart.y;
        double norm_inv = 1.0 / ( x_Nest2Food_d[i]*x_Nest2Food_d[i] + y_Nest2Food_d[i]*y_Nest2Food_d[i] );
        x_Nest2Food_d[i] = x_Nest2Food_d[i] * norm_inv;
        y_Nest2Food_d[i] = y_Nest2Food_d[i] * norm_inv;
    }

    gauss_path<<<MACRO_MAX,MACRO_MAX>>>();
}

__host__ void reset(double sensor,int naho,unsigned long long int step){
    cudaMemcpyToSymbol(SENSOR2,&sensor,sizeof(double),0);
    cudaMemcpyToSymbol(NAHO,&naho,sizeof(int),0);

    //initialize();
    //antsInit<<<MACRO_NMAX,1>>>();
    //cellsInit<<<MACRO_MAX,MACRO_MAX>>>();

    //setEdges<<<MACRO_MAX,MACRO_MAX>>>();
    //setNest<<<MACRO_MAX,MACRO_MAX>>>();
    //setDistFromNest<<<MACRO_MAX,MACRO_MAX>>>();

    //setNestDirs<<<MACRO_MAX,MACRO_MAX>>>();
    //setFoodsDir<<<MACRO_NUM_FOODS,1>>>();

    srand(MACRO_RND_SEED+step);

    thrust::host_vector<unsigned long long int> seeds_vec_h(MACRO_NMAX);
    std::generate(seeds_vec_h.begin(), seeds_vec_h.end(), rand);
    thrust::copy(seeds_vec_h.begin(), seeds_vec_h.end(), seeds_d_ptr);
    randInit<<<MACRO_NMAX,1>>>();

    antsReset<<<MACRO_NMAX,1>>>();
    cellsReset<<<MACRO_MAX,MACRO_MAX>>>();
    foodsReset<<<MACRO_NUM_FOODS,1>>>();
}
