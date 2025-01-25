//
// Created by iorio on 7/25/22.
//

#include <pairinstability.h>
#include <star.h>

///Base class
PairInstability *PairInstability::Instance(const std::string &name) {
    auto it = GetStaticMap().find(name);
    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;

    if (it != GetStaticMap().end()){
        GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;
    }
    //otherwise return an error
    else{
        SevnLogging svlog;
        std::string available_options;
        for (auto const& option : GetStaticMap()) available_options+=option.first+", ";
        svlog.critical("Option "+name+" not available for Process PairInstability. The available options are: \n"+available_options,
                       __FILE__,__LINE__,sevnstd::notimplemented_error(" "));
    }


    return (it->second)->instance(); //when I call the instance create a new object with the non trivial constructor
}

void PairInstability::Register(PairInstability *ptr){
    //Register only if the name is not already in the map,
    //The check is necessary because, the ctor call a register by default, therefore
    //without this check each call to a PairInstability ctor will add a new record on the map
    if (GetStaticMap().find(ptr->name())==GetStaticMap().end()) {
        GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size() - 1] = 0;
    }
}

//Simple standard beforeSN, just copy the mass values
utilities::MassContainer PairInstability::apply_beforeSN(Star *s)  const {
    return utilities::MassContainer{s->getp(Mass::ID),s->getp(MHE::ID), s->getp(MCO::ID)};
}

//Just set Mremnant_after_pi=mremnant and SNType=CoreCollapse
PISNreturn PairInstability::apply_afterSN(_UNUSED Star *s, double mremnant) const {
    return  PISNreturn{mremnant, Lookup::SNExplosionType::CoreCollapse};
}


///Disabled
PIDisabled::PIDisabled(bool reg) : PairInstability(false) {
    if (reg){
        Register(this);
    }
}


