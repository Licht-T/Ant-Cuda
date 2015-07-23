#include "IO.h"

int pw_old;

const std::string constHeader("Constants.h");
std::string path;
std::string fool_num_plot;
std::ofstream *ofs;

template <int max> struct getHomingFood{
	int operator()(int i){
		static const int n = max-1;
		static thrust::plus<int> binary_op;
		if(n==i){
			return thrust::transform_reduce(ants_d_ptr, ants_d_ptr+NMAX, getHomingWithFoodNum<n>(), 0, binary_op);;
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
			return thrust::transform_reduce(ants_d_ptr, ants_d_ptr+NMAX, getHomingWithTypeAndFoodNum<ch,n>(), 0, binary_op);;
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

template <class T> std::string toString(T t){
	static std::ostringstream ss;
	ss << (T)t;

	std::string str(ss.str());
	ss.str("");
	ss.clear(std::stringstream::goodbit);

	return str;
}

void IOInit(){	

	struct stat st;

	pw_old = -1;

	std::string nmax = toString(NMAX);
	std::string max = toString(MAX);

	std::string initDir(nmax+"ants_"+max+"x"+max+"cells");
	if(stat(initDir.c_str(), &st) != 0){
		mkdir(initDir.c_str(), 0775);
	}

	std::string fnum = toString(NUM_FOODS);

	std::string fNumDir(initDir+"/"+fnum+"foodnum");
	if(stat(fNumDir.c_str(), &st) != 0){
		mkdir(fNumDir.c_str(), 0775);
	}

	std::string fsource = toString(FOODSOURCE);
	std::string fdist = toString(FOOD_DIST);

	std::string fCondDir(fNumDir+"/"+fsource+"initfvol_"+fdist+"fdist");
	if(stat(fCondDir.c_str(), &st) != 0){
		mkdir(fCondDir.c_str(), 0775);
	}

	std::string step = toString(MAX_STEP);
	std::string angle = toString(FOOD_ANGLE);

	std::string stepAngleDir(fCondDir+"/"+step+"steps_"+angle+"deg");
	if(stat(stepAngleDir.c_str(), &st) != 0){
		mkdir(stepAngleDir.c_str(), 0775);
	}

	path = std::string(stepAngleDir+"/");
	fool_num_plot = std::string(path+step+"steps_"+angle+"deg.dat");

	ofs = new std::ofstream(fool_num_plot.c_str());

	std::ifstream constifs(constHeader.c_str());
	std::ofstream constofs((path+constHeader).c_str());
	constofs << constifs.rdbuf() << std::flush;
}

void IOEffWrite(int pw, int n, double sum){
	if(pw_old<pw){
		pw_old = pw;
		(*ofs) << std::endl;
	}
	(*ofs)  << pw << " "
		<< n  << " "
		<< (sum/(MAX_STEP))/(MAX_TIME-1000)
		<< std::endl;
}

void IOCellWrite(int pw, int n){
	std::string pwstr = toString(pw);
	std::string nstr = toString(n);
	std::string anglestr = toString(FOOD_ANGLE);

	std::string celldata(path+"cell_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+".dat");
	std::ofstream cellfs(celldata.c_str());

	cudaMemcpyFromSymbol(cells,cells_d,MAX*MAX*sizeof(Cell),0);
	for(int i=0; i<MAX; i++){
		for(int j=0; j<MAX; j++){
			cellfs << cells[j][i].cart.x << " "
				<< cells[j][i].cart.y << " "
				<< cells[j][i].phero
				<< std::endl;
		}
		cellfs << std::endl;
	}
}

void IOEffPoll(int pw, int n, int sample, int t){

	static getHomingWithType<NORMAL_CH> normalOp;
	static getHomingWithType<FOOL_CH> ahoOp;
	static thrust::plus<int> binary_op;
	static int homingFoods[NUM_FOODS];
	static int foolHomingFoods[NUM_FOODS];
	static int normalHomingFoods[NUM_FOODS];
	static getHomingFood<NUM_FOODS> homingFoodFunctor;
	static getHomingTypeAndFood<FOOL_CH,NUM_FOODS> foolHomingFoodFunctor;
	static getHomingTypeAndFood<NORMAL_CH,NUM_FOODS> normalHomingFoodFunctor;
	static thrust::host_vector<double> phero_h(NUM_FOODS);
	static thrust::device_vector<double> phero_d(NUM_FOODS);	

	std::string stepstr = toString(MAX_STEP);
	std::string anglestr = toString(FOOD_ANGLE);

	std::string pwstr = toString(pw);
	std::string nstr = toString(n);
	std::string samplestr = toString(sample);

	std::string celldata(path+"food_"+anglestr+"deg_10e-"+pwstr+"_"+nstr+"normal"+"_sampleNo"+samplestr+"_of_"+stepstr+".dat");
	std::ofstream pollfs(celldata.c_str(),std::ios::out | std::ios::app);

	int nor = thrust::transform_reduce(ants_d_ptr, ants_d_ptr+NMAX, normalOp, 0, binary_op);
	int aho = thrust::transform_reduce(ants_d_ptr, ants_d_ptr+NMAX, ahoOp, 0, binary_op);

	thrust::transform(foods_d_ptr, foods_d_ptr+NUM_FOODS, phero_d.begin(), getPheroAroundFood());
	thrust::copy(phero_d.begin(), phero_d.end(), phero_h.begin());
	
	pollfs  << t << " ";

	for (int i=0; i<NUM_FOODS; i++){
		//homingFoods[i] = homingFoodFunctor(i);
		foolHomingFoods[i] = foolHomingFoodFunctor(i);
		normalHomingFoods[i] = normalHomingFoodFunctor(i);
		homingFoods[i]=foolHomingFoods[i]+normalHomingFoods[i];
	}

	for (int i=0; i<NUM_FOODS; i++){
		pollfs  << homingFoods[i]
			<< " ";
	}
	pollfs  << nor      << " "
		<< aho      << " "
		<< nor/(double)n << " "
		<< aho/(double)(NMAX-n) << " "
		<< (nor+aho)<< " ";
	for (int i=0; i<NUM_FOODS; i++){
		pollfs << normalHomingFoods[i] << " ";
	}
	for (int i=0; i<NUM_FOODS; i++){
		pollfs << foolHomingFoods[i] << " ";
	}
	for (int i=0; i<NUM_FOODS; i++){
		pollfs << normalHomingFoods[i]/(double)n << " ";
	}
	for (int i=0; i<NUM_FOODS; i++){
		pollfs << foolHomingFoods[i]/(double)(NMAX-n) << " ";
	}
	pollfs << (nor+aho) << " ";
	for (int i=0; i<NUM_FOODS; i++){
		pollfs << phero_h[i] << " ";
	}
	pollfs << std::endl;
}
