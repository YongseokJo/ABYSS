//
// Created by spera on 11/02/19.
//
#include <lookup_and_phases.h>
#include <string>


std::string Lookup::literal(Lookup::Phases A){

    switch(A){
        case Lookup::PreMainSequence: return "PreMainSequence"; break;
        case Lookup::MainSequence: return "MainSequence"; break;
        case Lookup::TerminalMainSequence: return "TerminalMainSequence"; break;
        case Lookup::ShellHBurning: return "ShellHBurning"; break;
        case Lookup::CoreHeBurning: return "CoreHeBurning"; break;
        case Lookup::TerminalCoreHeBurning: return "TerminalCoreHeBurning"; break;
        case Lookup::ShellHeBurning: return "ShellHeBurning"; break;
        case Lookup::Remnant: return "Remnant"; break;
        default: return "Empty"; break;

    }

}

std::string Lookup::literal(Lookup::Remnants A){

    if(A == Lookup::NotARemnant) return "NotARemnant";
    else if (A == Lookup::HeWD) return "HeWD";
    else if (A == Lookup::COWD) return "COWD";
    else if (A == Lookup::ONeWD) return "ONeWD";
    else if (A == Lookup::NS_ECSN) return "NS_ECSN";
    else if (A == Lookup::NS_CCSN) return "NS_CCSN";
    else if (A == Lookup::BH) return "BH";
    else if (A == Lookup::Empty) return "Empty";
    else if (A == Lookup::Zombie) return "Zombie";
    else return "Unknown_state";

}

std::string Lookup::literal(Lookup::EventsList A){

    if(A == Lookup::NoEvent) return "";
    else if (A == Lookup::ChangePhase) return "CPhase";
    else if (A == Lookup::ChangeRemnant) return "CRem";
    else if (A == Lookup::QHE) return "QHE";
    else if (A == Lookup::GWBegin) return "GWB";
    else if (A == Lookup::RLOBegin) return "RLOB";
    else if (A == Lookup::RLOEnd) return "RLOE";
    else if (A == Lookup::Collision) return "Collision";
    else if (A == Lookup::Merger) return "Merger";
    else if (A == Lookup::CE) return "CE";
    else if (A == Lookup::CE_Merger) return "CE_Merger";
    else if (A == Lookup::RLOB_Merger) return "RLOB_Merger";
    else if (A == Lookup::RLOB_CE) return "RLOB_CE";
    else if (A == Lookup::RLOB_CE_Merger) return "RLOB_CE_Merger";
    else if (A == Lookup::Collision_Merger) return "Collision_Merger";
    else if (A == Lookup::Collision_CE) return "Collision_CE";
    else if (A == Lookup::Collision_CE_Merger) return "Collision_CE_Merger";
    else if (A == Lookup::Swallowed) return "Swallowed";
    else if (A == Lookup::RLOB_Swallowed) return "RLOB_Swallowed";
    else if (A == Lookup::GW_Merger) return "GW_Merger";
    else if (A == Lookup::SNBroken) return "SNBroken";
    else if (A == Lookup::SNIa) return "SNIa";
    else return "Unknown_state";

}


const Lookup::FILEMAP Lookup::filemap = {
        {"mass.dat", Lookup::Tables::_Mass},
        {"lumi.dat", Lookup::Tables::_Lumi},
        {"mco.dat", Lookup::Tables::_MCO},
        {"mhe.dat", Lookup::Tables::_MHE},
        {"phase.dat", Lookup::Tables::_Phase},
        {"radius.dat", Lookup::Tables::_Radius},
        {"time.dat", Lookup::Tables::_Time}};

const Lookup::FILEMAP Lookup::filemap_optional = {
        {"rco.dat", Lookup::Tables::_RCO},
        {"rhe.dat", Lookup::Tables::_RHE},
        {"inertia.dat", Lookup::Tables::_Inertia},
        {"hsup.dat", Lookup::Tables::_Hsup},
        {"hesup.dat", Lookup::Tables::_HEsup},
        {"csup.dat", Lookup::Tables::_Csup},
        {"nsup.dat", Lookup::Tables::_Nsup},
        {"osup.dat", Lookup::Tables::_Osup},
        {"qconv.dat", Lookup::Tables::_Qconv},
        {"depthconv.dat", Lookup::Tables::_Depthconv},
        {"tconv.dat", Lookup::Tables::_Tconv}};

//GI 81219: Map to manage the output options
const Lookup::OUTPUTMAP Lookup::outputmap{
        {"h5", OutputOption::_hdf5},
        {"hdf5", OutputOption::_hdf5},
        {"hf5", OutputOption::_hdf5},
        {"hd5", OutputOption::_hdf5},
        {"txt", OutputOption::_ascii},
        {"dat", OutputOption::_ascii},
        {"text", OutputOption::_ascii},
        {"ascii", OutputOption::_ascii},
        {"asci", OutputOption::_ascii},
        {"csv", OutputOption::_csv},
        {"binary", OutputOption::_binary}};


const Lookup::INPUTMAPBIN_NPARAM Lookup::inputmapbin_nparam{
        {InputBinaryOption::_new, 14},
        {InputBinaryOption::_legacy, 14}
};

//GI 81219: Map to manage the output options
const Lookup::INPUTMAPBIN Lookup::inputmapbin{
        {"new",  std::make_pair (InputBinaryOption::_new, inputmapbin_nparam.at(InputBinaryOption::_new)) },
        {"legacy", std::make_pair (InputBinaryOption::_legacy, inputmapbin_nparam.at(InputBinaryOption::_legacy))},
        {"sevn1", std::make_pair (InputBinaryOption::_legacy, inputmapbin_nparam.at(InputBinaryOption::_legacy))},
        {"old", std::make_pair (InputBinaryOption::_legacy, inputmapbin_nparam.at(InputBinaryOption::_legacy))}
        };



const Lookup::RLMAP Lookup::rlmap{
        {"hurley", RLMode::_RLHurley},
        {"hurley_rl", RLMode::_RLHurley},
        {"h", RLMode::_RLHurley},
        {"hurley_rl_prop", RLMode::_RLHurleyprop},
        {"hurley_rl_propeller", RLMode::_RLHurleyprop},
        {"hurley_bse", RLMode::_RLHurleyBSE},
        {"hurley_rl_bse", RLMode::_RLHurleyBSE},
        {"hbse", RLMode::_RLHurleyBSE},
        {"hurleymod", RLMode::_RLHurleymod},
        {"hurleymod_rl", RLMode::_RLHurleymod},
        {"hm", RLMode::_RLHurleymod},
        {"disabled", RLMode::_RLdisabled},
        {"disabled_rl", RLMode::_RLdisabled},
        {"d", RLMode::_RLdisabled}

};


const Lookup::TIDESMAP Lookup::tidesmap{
        {"tides_simple", TidesMode::_Tsimple},
        {"simple", TidesMode::_Tsimple},
        {"s", TidesMode::_Tsimple},
        {"ts", TidesMode::_Tsimple},
        {"tides_simple_notab", TidesMode::_Tsimple_notab},
        {"simple_notab", TidesMode::_Tsimple_notab},
        {"snt", TidesMode::_Tsimple_notab},
        {"tsnt", TidesMode::_Tsimple_notab},
        {"disabled", TidesMode::_Tdisabled},
        {"d", TidesMode::_Tdisabled}

};


const Lookup::GWMAP Lookup::gwmap{
        {"disabled", GWMode::_GWdisabled},
        {"d", GWMode::_GWdisabled},
        {"standard",GWMode::_GWPeters},
        {"s",GWMode::_GWPeters},
        {"peters",GWMode::_GWPeters},
        {"p",GWMode::_GWPeters}
};

const Lookup::MIXMAP Lookup::mixmap{
        {"disabled", MixMode::_Mixdisabled},
        {"d", MixMode::_Mixdisabled},
        {"simple",MixMode::_Mixsimple},
        {"s",MixMode::_Mixsimple},
};

const Lookup::SNKMAP Lookup::snkmap{
        {"disabled", SNKickMode::_SNKickdisabled},
        {"d", SNKickMode::_SNKickdisabled},
        {"hurley",SNKickMode::_SNKickHurley},
        {"h",SNKickMode::_SNKickHurley},
        {"hurley_kick",SNKickMode::_SNKickHurley},
        {"hk",SNKickMode::_SNKickHurley},
};

const Lookup::CEMAP Lookup::cemap{
        {"disabled", CEMode::_CEdisabled},
        {"d", CEMode::_CEdisabled},
        {"hurley",CEMode::_CEEnergy},
        {"h",CEMode::_CEEnergy},
        {"hurley_ce",CEMode::_CEEnergy},
        {"hce",CEMode::_CEEnergy},
        {"default",CEMode::_CEEnergy},
        {"energy",CEMode::_CEEnergy},
};


