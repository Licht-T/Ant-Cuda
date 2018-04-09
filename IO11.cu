#include "IO11.h"

int pw_old;

const std::string constHeader("Constants.h");
std::string path;
std::string fool_num_plot;
std::ofstream *ofs;
std::string fool_num_plot_prob;
std::ofstream *ofsProb;
std::string fool_num_food;
std::ofstream *ofsfood;

template <int max> struct getHomingFood{
    int operator()(int i){
        static const int n = max-1;
        static thrust::plus<int> binary_op;
        if(n==i){
            return thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, getHomingWithFoodNum<n>(), 0, binary_op);;
        }
        else{
            return getHomingFood<n>()(i);
        }
    }
};

template <> struct getHomingFood<0>{
    int operator()(int i){
        return -1;
    }
};

template <AntCharacter ch,int max> struct getHomingTypeAndFood{
    int operator()(int i){
        static const int n = max-1;
        static thrust::plus<int> binary_op;
        if(n==i){
            return thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, getHomingWithTypeAndFoodNum<ch,n>(), 0, binary_op);;
        }
        else{
            return getHomingTypeAndFood<ch,n>()(i);
        }
    }
};

template <AntCharacter ch> struct getHomingTypeAndFood<ch,0>{
    int operator()(int i){
        return -1;
    }
};

struct getPheroAroundFood{
    __device__ double operator()(const Food food) const{
        int i = food.i;
        int j = food.j;
        double sum = 0.0;
        sum += cells_d[i][j].phero;

        Cell *c=NULL;
        for(enum Direction dir = UP; dir<=UPLEFT; (dir<<=1) ){
            c = getCell(cells_d,i,j,dir);
            sum += c->phero;
        }
        return sum/7.0;
    }
};

// template <class T> std::string toString(T t){
//     static std::ostringstream ss;
//     ss << (T)t;

//     std::string str(ss.str());
//     ss.str("");
//     ss.clear(std::stringstream::goodbit);

//     return str;
// }

void IOInit(){	

    struct stat st;

    pw_old = -1;

    std::string nmax = toString(MACRO_NMAX);
    std::string max = toString(MACRO_MAX);

    std::string initDir(nmax+"ants_"+max+"x"+max+"cells");
    if(stat(initDir.c_str(), &st) != 0){
#ifdef _WIN32
        _mkdir(initDir.c_str());
#else
        mkdir(initDir.c_str(), 0775);
#endif
    }

    std::string fnum = toString(MACRO_NUM_FOODS);

    std::string fNumDir(initDir+"/"+fnum+"foodnum");
    if(stat(fNumDir.c_str(), &st) != 0){
#ifdef _WIN32
        _mkdir(fNumDir.c_str());
#else
        mkdir(fNumDir.c_str(), 0775);
#endif
    }

    std::string fsource = toString(MACRO_FOODSOURCE);
    std::string fdist = toString(MACRO_FOOD_DIST);

    std::string fCondDir(fNumDir+"/"+fsource+"initfvol_"+fdist+"fdist");
    if(stat(fCondDir.c_str(), &st) != 0){
#ifdef _WIN32
        _mkdir(fCondDir.c_str());
#else
        mkdir(fCondDir.c_str(), 0775);
#endif
    }

    std::string step = toString(MACRO_MAX_STEP);
    std::string angle = toString(MACRO_FOOD_ANGLE);

    std::string stepAngleDir(fCondDir+"/"+step+"steps_"+angle+"deg");
    if(stat(stepAngleDir.c_str(), &st) != 0){
#ifdef _WIN32
        _mkdir(stepAngleDir.c_str());
#else
        mkdir(stepAngleDir.c_str(), 0775);
#endif
    }

    path = std::string(stepAngleDir+"/");
    std::string rec = toString(MACRO_REC);
    fool_num_plot = std::string(path+step+"steps_"+angle+"deg_"+rec+"rec.dat");

    ofs = new std::ofstream(fool_num_plot.c_str());

    fool_num_plot_prob = std::string(path+step+"steps_"+angle+"deg_"+rec+"rec_prob.dat");
    ofsProb = new std::ofstream(fool_num_plot_prob.c_str());

    std::ifstream constifs(constHeader.c_str());
    std::ofstream constofs((path+constHeader).c_str());
    constofs << constifs.rdbuf() << std::flush;


    // fool_num_food = std::string(path+step+"steps_"+angle+"deg_food.dat");
    // ofsfood = new std::ofstream(fool_num_food.c_str());
}

void IOEffWrite(int pw, int n, double sum){
    if(pw_old<pw){
        pw_old = pw;
        (*ofs) << std::endl;
    }
    (*ofs)  << pw << " "
        << n  << " "
        << (sum/(MACRO_MAX_STEP))/(MACRO_MAX_TIME-1000)
        << std::endl;
}

void IOProbWrite(int pw, int n, double prob[]){
    if(pw_old<pw){
        pw_old = pw;
        (*ofsProb) << std::endl;
    }
    (*ofsProb)  << pw << " "
        << n  << " "
        << ((double)prob[0]/(MACRO_MAX_STEP))/(MACRO_MAX_TIME) << " "
        << ((double)prob[1]/(MACRO_MAX_STEP))/(MACRO_MAX_TIME) << " "
        << ((double)prob[2]/(MACRO_MAX_STEP))/(MACRO_MAX_TIME) << " "
        << ((double)prob[3]/(MACRO_MAX_STEP))/(MACRO_MAX_TIME) << " "
        << ((double)prob[4]/(MACRO_MAX_STEP))/(MACRO_MAX_TIME)
        << std::endl;
}

void IOFinWrite(int pw, int n, double sum){
    if(pw_old<pw){
        pw_old = pw;
        (*ofs) << std::endl;
    }
    (*ofs)  << pw << " "
        << n  << " "
        << (sum/(MACRO_MAX_STEP))
        << std::endl;
}

void IOFoodWrite(int pw, int n, double ft[], double et[]){
    if(pw_old<pw){
        pw_old = pw;
        (*ofsfood) << std::endl;
    }
    (*ofsfood)  << pw << " " << n  << " ";

    for(int id=0; id<MACRO_NUM_FOODS; id++)
        (*ofsfood) << ft[id]/(MACRO_MAX_STEP) << " ";

    for(int id=0; id<MACRO_NUM_FOODS; id++)
        for(int id2=id+1; id2<MACRO_NUM_FOODS; id2++)
            (*ofsfood) << (ft[id] - ft[id2])/(MACRO_MAX_STEP) << " ";

    for(int id=0; id<MACRO_NUM_FOODS; id++)
        (*ofsfood) << et[id]/(MACRO_MAX_STEP) << " ";

    for(int id=0; id<MACRO_NUM_FOODS; id++)
        for(int id2=id+1; id2<MACRO_NUM_FOODS; id2++)
            (*ofsfood) << (et[id] - et[id2])/(MACRO_MAX_STEP) << " ";

    (*ofsfood) << std::endl;
}

void IOFoodAmountWrite(std::ofstream & ofs, std::vector< std::vector<double> > & ft){
    // if(pw_old<pw){
    //     pw_old = pw;
    //     (*ofs) << std::endl;
    // }

    for(int t=0; t<MACRO_MAX_TIME; t++){
        (ofs)  << t << "\t" ;
        for(int id=0; id<MACRO_NUM_FOODS; id++)
            (ofs) << ft[t][id] << "\t";
        (ofs) << std::endl;
    }
}


// void IOCellWrite(int pw, int n){
//     std::string pwstr = toString(pw);
//     std::string nstr = toString(n);
//     std::string anglestr = toString(MACRO_FOOD_ANGLE);

//     std::string celldata(path+"cell_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+".dat");
//     std::ofstream cellfs(celldata.c_str());

//     cudaMemcpyFromSymbol(cells,cells_d,MACRO_MAX*MACRO_MAX*sizeof(Cell),0);
//     for(int i=0; i<MACRO_MAX; i++){
//         for(int j=0; j<MACRO_MAX; j++){
//             cellfs << cells[j][i].cart.x << " "
//                 << cells[j][i].cart.y << " "
//                 << cells[j][i].phero
//                 << std::endl;
//         }
//         cellfs << std::endl;
//     }
// }

void IOCellWrite(int pw, int n, int sample, int t){
    std::string pwstr = toString(pw);
    std::string nstr = toString(n);
    std::string anglestr = toString(MACRO_FOOD_ANGLE);
    std::string samplestr = toString(sample);
    std::string stepstr = toString(t);

    std::string celldata(path+"cell_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_sampleNo"+samplestr+"_of_"+stepstr+".dat");
    std::ofstream cellfs(celldata.c_str());

    cudaMemcpyFromSymbol(cells,cells_d,MACRO_MAX*MACRO_MAX*sizeof(Cell),0);
    for(int i=0; i<MACRO_MAX; i++){
        for(int j=0; j<MACRO_MAX; j++){
            cellfs << cells[j][i].cart.x << " "
                << cells[j][i].cart.y << " "
                << cells[j][i].phero
                << std::endl;
        }
        cellfs << std::endl;
    }

    std::string min1data(path+"min1_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_sampleNo"+samplestr+"_of_"+stepstr+".dat");
    std::ofstream min1fs(min1data.c_str());

    cudaMemcpyFromSymbol(min_path1,min_path1_d,MACRO_MAX*MACRO_MAX*sizeof(Cell),0);
    for(int i=0; i<MACRO_MAX; i++){
        for(int j=0; j<MACRO_MAX; j++){
            min1fs << min_path1[j][i].cart.x << " "
                << min_path1[j][i].cart.y << " "
                << min_path1[j][i].phero
                << std::endl;
        }
        min1fs << std::endl;
    }

    std::string min2data(path+"min2_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_sampleNo"+samplestr+"_of_"+stepstr+".dat");
    std::ofstream min2fs(min2data.c_str());

    cudaMemcpyFromSymbol(min_path2,min_path2_d,MACRO_MAX*MACRO_MAX*sizeof(Cell),0);
    for(int i=0; i<MACRO_MAX; i++){
        for(int j=0; j<MACRO_MAX; j++){
            min2fs << min_path2[j][i].cart.x << " "
                << min_path2[j][i].cart.y << " "
                << min_path2[j][i].phero
                << std::endl;
        }
        min2fs << std::endl;
    }
}

void IOPheroStateWrite(int pw, int n, int sample, double conv [][3]){
    std::string pwstr = toString(pw);
    std::string nstr = toString(n);
    std::string anglestr = toString(MACRO_FOOD_ANGLE);
    std::string samplestr = toString(sample);

    std::string pherostatedata(path+"cell_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_sampleNo"+samplestr+"_phero_state.dat");
    std::ofstream pherostate_fs(pherostatedata.c_str(), std::ios::app);

    for (int i=0; i<10; i++)
        pherostate_fs << i*1000 << " " << conv[i][0] << " " << conv[i][1] << " " << conv[i][2] << std::endl;
}

void IOPheroPathWrite(int pw, int n, int sample, int t, std::vector< std::vector< std::vector<double> > > & pathxyp){
    std::string pwstr = toString(pw);
    std::string nstr = toString(n);
    std::string anglestr = toString(MACRO_FOOD_ANGLE);
    std::string samplestr = toString(sample);
    std::string stepstr = toString(t);

    std::string pheropathdata(path+"cell_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_sampleNo"+samplestr+"_of_"+stepstr+"_phero_path.dat");
    std::ofstream pheropath_fs(pheropathdata.c_str(), std::ios::app);

    for (int fi=0; fi<MACRO_NUM_FOODS; fi++){
        for (int i=0; i<int(pathxyp[fi].size()); i++)
            pheropath_fs << i*1000 << " " << pathxyp[fi][i][0] << " " << pathxyp[fi][i][1] << " " << pathxyp[fi][i][2] << " " << pathxyp[fi][i][3] << std::endl;
        pheropath_fs << std::endl << std::endl;
    }
}

void IOPheroPathLengthWrite(int pw, int n, std::vector< std::vector< std::vector<double> > > & pathlen){
    std::string pwstr = toString(pw);
    std::string nstr = toString(n);
    std::string anglestr = toString(MACRO_FOOD_ANGLE);

    std::string pheropathlendata(path+"cell_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_phero_pathlen_len.dat");
    std::ofstream pheropathlen_fs(pheropathlendata.c_str(), std::ios::app);

    for (int i=0; i<MACRO_MAX_STEP; i++){
        pheropathlen_fs << i ;
        for (int fi=0; fi<int(pathlen.size()); fi++){
            for (int t=0; t<int(pathlen[fi].size()); t++){
                pheropathlen_fs << " " << pathlen[fi][t][i];
            }
        }
        pheropathlen_fs << std::endl;
    }
}

void IOEffPoll(int pw, int n, int sample, int t){

    static getHomingWithType<NORMAL_CH> normalOp;
    static getHomingWithType<FOOL_CH> ahoOp;
    static thrust::plus<int> binary_op;
    static int homingFoods[MACRO_NUM_FOODS];
    static int foolHomingFoods[MACRO_NUM_FOODS];
    static int normalHomingFoods[MACRO_NUM_FOODS];
    static getHomingFood<MACRO_NUM_FOODS> homingFoodFunctor;
    static getHomingTypeAndFood<FOOL_CH,MACRO_NUM_FOODS> foolHomingFoodFunctor;
    static getHomingTypeAndFood<NORMAL_CH,MACRO_NUM_FOODS> normalHomingFoodFunctor;
    static thrust::host_vector<double> phero_h(MACRO_NUM_FOODS);
    static thrust::device_vector<double> phero_d(MACRO_NUM_FOODS);

    std::string stepstr = toString(MACRO_MAX_STEP);
    std::string anglestr = toString(MACRO_FOOD_ANGLE);

    std::string pwstr = toString(pw);
    std::string nstr = toString(n);
    std::string samplestr = toString(sample);

    std::string celldata(path+"food_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_sampleNo"+samplestr+"_of_"+stepstr+".dat");
    std::ofstream pollfs(celldata.c_str(),std::ios::out | std::ios::app);

    int nor = thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, normalOp, 0, binary_op);
    int aho = thrust::transform_reduce(ants_d_ptr, ants_d_ptr+MACRO_NMAX, ahoOp, 0, binary_op);

    thrust::transform(foods_d_ptr, foods_d_ptr+MACRO_NUM_FOODS, phero_d.begin(), getPheroAroundFood());
    thrust::copy(phero_d.begin(), phero_d.end(), phero_h.begin());

    pollfs  << t << " ";

    for (int i=0; i<MACRO_NUM_FOODS; i++){
        //homingFoods[i] = homingFoodFunctor(i);
        foolHomingFoods[i] = foolHomingFoodFunctor(i);
        normalHomingFoods[i] = normalHomingFoodFunctor(i);
        homingFoods[i]=foolHomingFoods[i]+normalHomingFoods[i];
    }

    for (int i=0; i<MACRO_NUM_FOODS; i++){
        pollfs  << homingFoods[i]
            << " ";
    }
    pollfs  << nor      << " "
        << aho      << " "
        << nor/(double)n << " "
        << aho/(double)(MACRO_NMAX-n) << " "
        << (nor+aho)<< " ";
    for (int i=0; i<MACRO_NUM_FOODS; i++){
        pollfs << normalHomingFoods[i] << " ";
    }
    for (int i=0; i<MACRO_NUM_FOODS; i++){
        pollfs << foolHomingFoods[i] << " ";
    }
    for (int i=0; i<MACRO_NUM_FOODS; i++){
        pollfs << normalHomingFoods[i]/(double)n << " ";
    }
    for (int i=0; i<MACRO_NUM_FOODS; i++){
        pollfs << foolHomingFoods[i]/(double)(MACRO_NMAX-n) << " ";
    }
    pollfs << (nor+aho) << " ";
    for (int i=0; i<MACRO_NUM_FOODS; i++){
        pollfs << phero_h[i] << " ";
    }
    pollfs << std::endl;
}
