//
// Created by Giuliano Iorio on 26/07/2022.
//

#include <neutrinomassloss.h>

void NeutrinoMassLoss::Register(NeutrinoMassLoss *ptr) {
    //Register only if the name is not already in the map,
    //The check is necessary because, the ctor call a register by default, therefore
    //without this check each call to a NeutrinoMassLoss ctor will add a new record on the map
    if (GetStaticMap().find(ptr->name())==GetStaticMap().end()) {
        GetStaticMap().insert(std::make_pair(ptr->name(), ptr));
        GetUsed().resize(GetStaticMap().size());
        GetUsed()[GetUsed().size() - 1] = 0;
    }
}

NeutrinoMassLoss *NeutrinoMassLoss::Instance(const std::string &name)  {
    auto it = GetStaticMap().find(name);
    if(it!= GetStaticMap().end()) GetUsed()[std::distance(GetStaticMap().begin(), it)] = 1;

    //If name in hte map set the GetUsed vector
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

NMLDisabled::NMLDisabled(bool reg) : NeutrinoMassLoss(false) {
    if (reg){
        Register(this);
    }
}
