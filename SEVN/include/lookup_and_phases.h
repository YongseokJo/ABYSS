//
// Created by spera on 11/02/19.
//

#ifndef SEVN_LOOKUP_AND_PHASES_H
#define SEVN_LOOKUP_AND_PHASES_H

#include <utility>
#include <vector>
#include <map>
#include <string>

class Orbital_change;



namespace Lookup{

    enum Tables {
        _Time = 0,
        _Mass,
        _Radius,
        _MHE,
        _MCO,
        _RHE,
        _RCO,
        _Lumi,
        _Phase,
        _Inertia,
        _Hsup,
        _HEsup,
        _Csup,
        _Nsup,
        _Osup,
        _Qconv,
        _Depthconv,
        _Tconv,
        _Ntables
    };

    ///NOTICE: The order is important from less evolved to more evolved
    enum Phases {
        PreMainSequence = 0,
        MainSequence,
        TerminalMainSequence,
        ShellHBurning,
        CoreHeBurning,
        TerminalCoreHeBurning,
        ShellHeBurning,
        Remnant,
        Nphases
    };

    enum  SNExplosionType{
        Unknown = 0,
        ElectronCapture,
        CoreCollapse,
        PPISN,
        PISN,
        Ia,
        WDformation,
        NSNexplosiontypes
    };

    enum Material {
        H = 0,
        He,
        CO,
        ONe,
        Neutron,
        BHmaterial,
        Unknownmaterial,
        Nmaterial
    };

    //Different remnants have different way to handle property::remnant
    ///NOTICE, the order is important (mass ordered): WD, NS, BH
    enum Remnants {
        Zombie=-2,
        Empty=-1,
        NotARemnant,
        HeWD,
        COWD,
        ONeWD,
        NS_ECSN,
        NS_CCSN,
        BH,
        Nremnants
    };

    enum OutputOption {
        _hdf5 = 0,
        _ascii,
        _csv,
        _binary,
        _Noption
    };

    enum WindsMode {
        _WHurley = 0,
        _WFaniAD,
        _WFaniDE,
        _Wdisabled
    };

    enum RLMode {
        _RLHurley = 0,
        _RLHurleyBSE,
        _RLHurleymod,
        _RLHurleyprop,
        _RLdisabled,
    };

    enum TidesMode {
        _Tsimple = 0,
        _Tsimple_notab,
        _Tdisabled
    };

    enum GWMode {
        _GWPeters = 0,
        _GWdisabled
    };

    enum MixMode {
        _Mixsimple = 0,
        _Mixdisabled,
    };

    enum SNKickMode {
        _SNKickHurley =0,
        _SNKickdisabled,
    };

    enum CEMode {
        _CEEnergy = 0,
        _CEdisabled
    };

    /*Vector containing available xspin modes*/
    const std::vector<std::string> xspinmodes = {"geneva", "mesa", "fuller", "maxwellian", "zeros", "accretion", "disabled"};

    /**
     * This Enum match the kind of input for binaries with the expected length of the input param vectors:
     *  - Option standard: new formalism in SEVN2 (len 14) -> Mass1, Z1, Spin1, t_1 initial, Mass2, Z2, Spin2, t_2 initial, sn_type, t final, dt_out, bin_separation, bin_eccentricity
     *  - Standard formalism in SEVN1 (len 13) ->   Mass1, Mass2, Z1, Z2, Spin1, Spin2, bin_separation, bin_eccentricity, t_final, t_ini, tstep  sn_type1, sn_type2, dt_out
     */
    enum InputBinaryOption {
        _new = 0,
        _legacy,
        _Niop
    };

    /**
     * This enum stores the Events code
     */
    enum EventsList {
        NoEvent = -1,
        ChangePhase,
        ChangeRemnant,
        QHE,
        GWBegin,
        RLOBegin,
        RLOEnd,
        Collision,
        CE,
        Merger,
        CE_Merger,
        RLOB_Merger,
        RLOB_CE,
        RLOB_CE_Merger,
        Collision_Merger,
        Collision_CE,
        Collision_CE_Merger,
        Swallowed,
        RLOB_Swallowed,
        GW_Merger,
        SNBroken,
        SNIa
    };

    /** Store the DTout type:
     * */
    enum DTout_type : unsigned int {
        _DTUNKNOWN =0, /*!< not know */
        _DTVAL, /*!< just a number */
        _DTLIST, /*!< just of values */
        _DTGENINTERVAL, /*!< three numbers setting an interval */
        _DTALL, /*!< print all */
        _DTEND, /*!< prìnt only end  */
        _DTEVNT, /*!< prìnt events  */
        _DTEVNTRLO /*!< print events + all steps during a RLO */
    };


    //TODO Make an input map also for single stellar evolution?


    //GI 81219: Map to handle the output options
    typedef std::map<std::string, OutputOption > OUTPUTMAP;
    extern const OUTPUTMAP outputmap;


    std::string literal(Phases A); //converts phases to strings
    std::string literal(Remnants A); //converts phases to strings
    std::string literal(EventsList A); //convert events to strings


    /***************************************************
     *    GI 230120: Map to handle the input options
     *    This part handle the input params for the binary evolution.
     *    All the possibile modality are reported in the enum InputBinaryOption.
     *    I attached to the modalities the expected number of parameters to read from the file
     *    to be able to properly handle the possible errors and exceptions.
     *    This info is stored in the map INPUTMAPBIN_PARAM.
     *    Finally the map INPUTMAPBIN maps a string that can be given as option in io.load to
     *    a pair containing the Input modality and the associated number of expected parameters.
     **************************************************/
    //GI 230120: Map to manage the input options
    typedef std::pair<InputBinaryOption,int> _INPAIR; //PAIR  INPUT OPTION - NUMBER OF EXPECTED INPUT PARAMETERS
    typedef std::map<InputBinaryOption, int> INPUTMAPBIN_NPARAM; //Map between Input option  and Number of expected input parameters
    typedef std::map<std::string, _INPAIR> INPUTMAPBIN;  //MAP string - PAIR (Input Option, N expected)
    //They are defined in cpp
    extern const INPUTMAPBIN inputmapbin;
    extern const INPUTMAPBIN_NPARAM inputmapbin_nparam;
    /***************************************************/


    //GI 210220: Map to handle the binary process options
    typedef std::map<std::string, WindsMode > WINDSMAP;
    typedef std::map<std::string, RLMode > RLMAP;
    typedef std::map<std::string, TidesMode > TIDESMAP;
    typedef std::map<std::string, GWMode > GWMAP;
    typedef std::map<std::string, RLMode > RLMAP;
    typedef std::map<std::string, MixMode > MIXMAP;
    typedef std::map<std::string, SNKickMode > SNKMAP;
    typedef std::map<std::string, CEMode > CEMAP;
    //These are defined in cpp
    extern const WINDSMAP windsmap;
    extern const RLMAP rlmap;
    extern const TIDESMAP tidesmap;
    extern const GWMAP gwmap;
    extern const RLMAP rlmap;
    extern const MIXMAP mixmap;
    extern const SNKMAP snkmap;
    extern const CEMAP cemap;
    //These are defined in static_main.h
    typedef std::map<WindsMode, std::string> WINDSMAP_NAME;
    extern const WINDSMAP_NAME windsmap_name;
    typedef std::map<TidesMode, std::string> TIDESMAP_NAME;
    extern const TIDESMAP_NAME tidesmap_name;
    typedef std::map<GWMode, std::string> GW_NAME;
    extern const GW_NAME gwmap_name;
    typedef std::map<RLMode, std::string> RL_NAME;
    extern const RL_NAME rlmap_name;
    typedef std::map<MixMode, std::string> MIX_NAME;
    extern const MIX_NAME mixmap_name;
    typedef std::map<SNKickMode, std::string> SNK_NAME;
    extern const SNK_NAME snkmap_name;
    typedef std::map<CEMode , std::string> CE_NAME;
    extern const CE_NAME cemap_name;


    typedef std::map<std::string, Tables> FILEMAP;
    extern const FILEMAP filemap; //filemap containing all the required tables
    extern const FILEMAP filemap_optional; //Filemap containing all the tables that are optional


}


#endif //SEVN_LOOKUP_AND_PHASES_H
