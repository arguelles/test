#include "Verosimilitud.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdexcept>
#include <Python.h>
#include "Flux.h"
#include "Tensor.h"
#include "H5Cpp.h"
#include "EffectiveArea.h"
#include "ConventionalFlux.h"
#include "ICData.h"
#include <cmath>
#include <cfloat>
#include <math.h>
#include <string>

#include <dlib/optimization.h>
#include <dlib/member_function_pointer.h>

int main()
{

std::shared_ptr<Verosimilitud> Vp = std::make_shared<Verosimilitud>(3,0,2);

//std::vector<double> Verosimilitud::MinLLH(std::vector<double> param, std::vector<double> low_bound, std::vector<double> high_bound, std::vector<bool> param_to_minimize)

double max = 4.0;

for (double prof_param = 0.01; prof_param < max; prof_param += 0.01){

	std::vector<double> param = {prof_param, 0.01, 1.0, 1.0};
	// std::vector<double> low_bound = {0.01, -0.05, 0.5, 0.875};
// 	std::vector<double> high_bound = {3.0, 0.05, 1.5, 1.125};
	
	std::vector<double> low_bound = {0.01, -0.15, 0.5, 0.875};
	std::vector<double> high_bound = {3.0, 10.00, 10.0, 10.0};

	std::vector<bool> param_to_minimize = {false, true, true, true};

	//std::cout << "================ BEGIN CODE =================" <<std::endl;

	dlib::matrix<double,0,1> r(4);
	for(size_t i=0; i<param.size(); i++) r(i)=param[i];
	// std::cout << "Chi2Initial: "<< std::endl << Vp->Chi2(r) <<std::endl;
	// std::cout << "Chi2GradientInitial: "<< std::endl << Vp->Chi2Gradient(r) <<std::endl;
	// 
	// std::cout << "================ BEGIN MINIMIZATION =================" <<std::endl;

	std::vector<double> result = Vp->MinLLH(param, low_bound, high_bound, param_to_minimize);

	// std::cout << "================ END MINIMIZATION =================" <<std::endl;

	// for(unsigned int i=0; i<result.size(); i++){
	// std::cout << result[i] << " ";
	// }

	std::cout<< prof_param << "," << result[result.size()-1] << "," ;


}
std::cout<<std::endl;













/*
std::shared_ptr<Verosimilitud> Vp = std::make_shared<Verosimilitud>(3,0,2);

//std::vector<double> Verosimilitud::MinLLH(std::vector<double> param, std::vector<double> low_bound, std::vector<double> high_bound, std::vector<bool> param_to_minimize)

std::vector<double> param = {1.0, 0.01, 1.0, 1.0};
std::vector<double> low_bound = {0.01, -0.05, 0.5, 0.875};
std::vector<double> high_bound = {3.0, 0.05, 1.5, 1.125};
std::vector<bool> param_to_minimize = {true, true, true, true};

std::cout << "================ BEGIN CODE =================" <<std::endl;

dlib::matrix<double,0,1> r(4);
for(size_t i=0; i<param.size(); i++) r(i)=param[i];
std::cout << "Chi2Initial: "<< std::endl << Vp->Chi2(r) <<std::endl;
std::cout << "Chi2GradientInitial: "<< std::endl << Vp->Chi2Gradient(r) <<std::endl;

std::cout << "================ BEGIN MINIMIZATION =================" <<std::endl;

std::vector<double> result = Vp->MinLLH(param, low_bound, high_bound, param_to_minimize);

std::cout << "================ END MINIMIZATION =================" <<std::endl;

for(unsigned int i=0; i<result.size(); i++){
std::cout << result[i] << " ";
}
std::cout << std::endl;

std::cout<<"I am done."<<std::endl;
*/



















/*
std::vector<double> LLH_values;
std::vector<double> norm_values;
std::vector<double> gamma_values;
std::vector<double> r_kpi_values;
std::vector<double> r_nubarnu_values;
std::vector<double> param;
*/


/*
//	loop over overall normalization values
LLH_values.clear();
norm_values.clear();


for( double	norm = 1.0; norm < 4.0; norm += 0.1)
{
	param.clear();

	param.push_back(norm);	// overall norm
	param.push_back(0.0);	// change to spectral index
	param.push_back(1.0);	// r_kpi
	param.push_back(1.0);	// r_nubarnu
	
	LLH_values.push_back(Vp->LLH(param));
	norm_values.push_back(norm);
	

}

for( int i = 0; i<LLH_values.size(); i++){
	std::cout<< norm_values[i] << "," << LLH_values[i]/550 << ",";
}
*/

/*
//  loop over spectral index changes

LLH_values.clear();
gamma_values.clear();

for( double	gamma = -1.0; gamma < 1.0; gamma += 0.05)
{

	param.clear();
	
	param.push_back(2.0);	// overall norm
	param.push_back(gamma);	// change to spectral index
	param.push_back(1.0);	// r_kpi
	param.push_back(1.0);	// r_nubarnu
	
	LLH_values.push_back(Vp->LLH(param));
	gamma_values.push_back(gamma);
		
}

for( int i = 0; i<LLH_values.size(); i++){
	std::cout<< gamma_values[i] << "," << LLH_values[i]/550 << ",";
}
*/

/*
//	loop over k/pi ratio
LLH_values.clear();
r_kpi_values.clear();


for( double	r_kpi = 0.5; r_kpi < 1.5; r_kpi += 0.05)
{
	param.clear();

	param.push_back(2.0);	// overall norm
	param.push_back(0.0);	// change to spectral index
	param.push_back(r_kpi);	// r_kpi
	param.push_back(1.0);	// r_nubarnu
	
	LLH_values.push_back(Vp->LLH(param));
	r_kpi_values.push_back(r_kpi);
	

}

for( int i = 0; i<LLH_values.size(); i++){
	std::cout<< r_kpi_values[i] << "," << LLH_values[i]/550 << ",";
}
*/

/*
loop over nubar/nu ratio
LLH_values.clear();
r_nubarnu_values.clear();


for( double	r_nubarnu = 0.875; r_nubarnu < 1.25; r_nubarnu += 0.01)
{
	param.clear();

	param.push_back(2.0);	// overall norm
	param.push_back(0.0);	// change to spectral index
	param.push_back(1.0);	// r_kpi
	param.push_back(r_nubarnu);	// r_nubarnu
	
	LLH_values.push_back(Vp->LLH(param));
	r_nubarnu_values.push_back(r_nubarnu);
	

}

for( int i = 0; i<LLH_values.size(); i++){
	std::cout<< r_nubarnu_values[i] << "," << LLH_values[i]/550 << ",";
}
*/



/*
	std::vector<double> param;
	param.push_back(1.0);	// overall norm
	param.push_back(0.0);	// change to spectral index
	param.push_back(1.0);	// r_kpi
	param.push_back(1.0);	// r_nubarnu


 std::cout << "LLH: "<< Vp->LLH(param) <<std::endl;
std::cout << "LLH: above!!!" <<std::endl;

std::vector<double> mc = Vp->GetPertExpectationVec(param);
std::vector<double> data = Vp->GetDataVec();

const int datasize = data.size();
const int mcsize = mc.size();

std::cout<< "MC: "<<std::endl;
for (int i = 0; i<mcsize; i++)
{
	std::cout<<mc[i]<<","<<std::endl;
}

std::cout<< "data: "<<std::endl;
for (int i = 0; i<datasize; i++)
{
	std::cout<<data[i]<<","<<std::endl;
}
*/


return 0;
	
}
