/*
 Created by iorio on 7/26/22.
 Base header for Kick models, it includes:
- Class Kicks: Base class for Kick models
*/

//TODO Add Documentation

#ifndef SEVN_KICKS_H
#define SEVN_KICKS_H

#include <string>
#include <map>
#include <vector>
#include <sevnlog.h>
#include <random>
#include <utilities.h>

class Star;
using sevnstd::SevnLogging;

//TODO Kicks here should be refactored, Kicks should be general enough to not contain gaussian_265
//then we can create a son class Maxwellian with methods that can be used by Hobbs Unified etc
class Kicks {

public:

    Kicks() : standard_gaussian(0,1), uniformRealDistribution(0,1){
    }

    virtual ~Kicks(){
    }

    static Kicks *Instance(std::string const &name);
    virtual Kicks *instance(){ svlog.critical("This is supposed to be pure virtual, do you have an uninitialized module?", __FILE__, __LINE__); return nullptr; }
    virtual inline std::string name() { return "Generic SN-kick model"; }

    /**
     * Wrapper for specified _apply functions
     * @param s
     */
    void apply(_UNUSED Star *s);

    virtual void _apply(_UNUSED Star *s) = 0;

    inline double get_random_kick(){
        if (random_velocity_kick==utilities::NULL_DOUBLE)
            svlog.critical("Trying to get random_kick before it is initiliased",__FILE__,__LINE__);

        return random_velocity_kick;
    }

    
protected:

    void Register(Kicks *ptr, const std::string &_name) {
        GetStaticMap().insert(std::make_pair(_name, ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size()-1] = 0;
    }

    std::normal_distribution<double> standard_gaussian;
    std::uniform_real_distribution<double> uniformRealDistribution;



    inline void set_random_kick(const double& a){
        random_velocity_kick =a;
        return;
    }
    void kick_initializer();
    SevnLogging svlog;
    double random_velocity_kick=utilities::NULL_DOUBLE;

    /**
     * Check if we have to make correction to the final Vkick (after all the fallback and similar correction)
     * It checks:
     *      - If Mremant=0 (e.g. after a PPISN)
     * So far it just check that the final Vkick is not lower than the parameter sn_min_vkick. If this is the case
     * it just sets the final velocity to the minimum value and rescale all the components by the factor min_vkick/old_vkick
     * if old_vkick is 0, new isotropic velocity components are randomly drawn
     * @param s
     */
    virtual void check_and_correct_vkick(Star* s);


    inline double draw_from_gaussian(double std, double mean=0.){
        return std*standard_gaussian(utilities::mtrand) + mean;
    }

private:

    static std::map<std::string, Kicks *> & GetStaticMap(){
        static std::map<std::string, Kicks *> _locmap;
        return _locmap;
    }
    static std::vector<int> & GetUsed(){
        static std::vector<int> _used;
        return _used;
    }


};



#endif //SEVN_KICKS_H
