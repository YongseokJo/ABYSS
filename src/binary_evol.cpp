//
// Created by iorio on 12/7/22.
//
#include <sevn.h>
#include <vector>
#include <string>


int binary (int argc, char **argv){
    //IO class
    IO sevnio;
    sevnio.load(argc,argv);

    std::vector<std::string> init_params{"40.52","0.002","0.0","rapid","tams","12.34","0.002","0.0","rapid","zams","1492.799","0.15","30000","all"};
    size_t ID=0;

    std::vector<double> time;
    std::vector<double> mass0;
    std::vector<double> mass1;
    std::vector<double> a;

    Binstar binary(&sevnio, init_params, ID);

    while(binary.getp(BWorldtime::ID)<1){
        binary.evolve();
        time.push_back(binary.getp(BWorldtime::ID));
        mass0.push_back(binary.getstar(0)->getp(Mass::ID));
        mass1.push_back(binary.getstar(1)->getp(Mass::ID));
        a.push_back(binary.getp(Semimajor::ID));
    }

    for (auto i=0; i< time.size(); i++){
        std::cout<<"Time [Myr]: "<<time[i]<<" Mass0 [Msun]: "<<mass0[i]<<" Mass1 [Msun]: "<<mass1[i]<<" Semiamjor [Rsun]: "<<a[i]<<std::endl;
    }

    std::cout<<SEVNinfo::get_full_info()<<std::endl;

    return EXIT_SUCCESS;
}
