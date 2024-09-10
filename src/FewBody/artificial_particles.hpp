#pragma once

#include "Common/binary_tree.h"
#include "tidal_tensor.hpp"
#include "pseudoparticle_multipole.hpp"
#include <cassert>
#include "Particle.h" // Eunwoo debug
#include "Group.h" // Eunwoo debug

typedef PseudoParticleMultipoleManager OrbitManager;


//! class to organize artificial particles
class ArtificialParticleManager{
public:
    // be careful to make consistent read/writeBinary and operator = if modified
    REAL r_tidal_tensor;  ///> tidal tensor maximum distance of particles
    int id_offset;       ///> offset to generate ID for artificial particles
    REAL gravitational_constant; ///> gravitational constant
    OrbitManager orbit_manager; ///> orbit particle method

    //! initializer, require setParticleSplitN later to initialize the size
    ArtificialParticleManager(): r_tidal_tensor(-1.0), id_offset(-1), gravitational_constant(1.0), orbit_manager() {}

    //! check paramters
    bool checkParams() {
        assert(TidalTensor::getParticleN()%2==0);
        assert(r_tidal_tensor>=0.0);
        assert(id_offset>0);
        assert(gravitational_constant>0);
        assert(orbit_manager.checkParams());
        return true;
    }

    //! create artificial particles 
    /*! First TidalTensor::getParticleN() are tidal tensor particles; others-1 are orbitial sample particles; last is c.m.  
      id: 
      tt/orb: id_offset + abs(member->id)*(n_artificial-1)/2 + index/2;
      cm: - abs(_bin.id)

      status:
      Tidal tensors/orbital
      0: left member N
      1: right member N
      others: index+1 (_data_to_store>0.0)

      c.m.: n_members

      mass_backup:
      Tidial tensors/orbital: 0.0 (_data_to_store <=0.0)
      c.m.:  mass(cm)

      @param[out]    _ptcl_new: artificial particles that will be added
      @param[in]     _bin: binary tree root
      @param[in]     _data_to_store: array of data to be stored in the status of artificial particles
      @param[in]     _n_data: number of data
    */
    // template <class Tptcl> // Eunwoo debug
    void createArtificialParticles(Particle* _ptcl_artificial,
                                   COMM::BinaryTree<Particle,COMM::Binary> &_bin, 
                                   const REAL* _data_to_store,
                                   const int _n_data) {
        ASSERT(checkParams());
        // set id and status.d except cm particle
        int n_artificial = getArtificialParticleN();
        for (int i=0; i<n_artificial-1; i++) {
            Particle* pi = &_ptcl_artificial[i];
            Particle* binary_member_i = _bin.getMember(i%2);
            // pi->PID = id_offset + abs(binary_member_i->PID)*n_artificial +i;
            auto& pi_artificial = pi->GroupInfo->artificial;
            pi_artificial.setParticleTypeToArtificial(REAL(i+1));
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(pi_artificial.isArtificial());
#endif
        }

        // First TidalTensor::getParticleN() is used for tidal tensor points
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(r_tidal_tensor<=_bin.changeover.getRin());
#endif
        TidalTensor::createTidalTensorMeasureParticles(_ptcl_artificial, *((Particle*)&_bin), r_tidal_tensor);

        // remaining is for orbital sample particles
        orbit_manager.createSampleParticles(&(_ptcl_artificial[getIndexOffsetOrb()]), _bin);
        
        // store the component member number 
        for (int j=0; j<2; j++) {
            int n_members = _bin.isMemberTree(j) ? ((COMM::BinaryTree<Particle,COMM::Binary>*)(_bin.getMember(j)))->getMemberN() : 1;
            _ptcl_artificial[j].GroupInfo->artificial.storeData(n_members); 
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(n_members>0);
#endif
        }

        // store the additional data (should be positive) 
        // ensure the data size is not overflow
        assert(_n_data+2<n_artificial-1);
        for (int j=0; j<_n_data; j++) {
            _ptcl_artificial[j+2].GroupInfo->artificial.storeData(_data_to_store[j]);
#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(_data_to_store[j]>0);
#endif
        }

        // last member is the c.m. particle
        Particle* pcm;
        pcm = &_ptcl_artificial[getIndexOffsetCM()];
        pcm->GroupInfo->artificial.setParticleTypeToCM(_bin.Mass, _bin.getMemberN());
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(pcm->GroupInfo->artificial.isCM());
#endif        
        if (getOrbitalParticleN()>0) pcm->Mass = 0.0;
        else pcm->Mass = _bin.Mass;
        pcm->Position[0] = _bin.Position[0];
        pcm->Position[1] = _bin.Position[1];
        pcm->Position[2] = _bin.Position[2];
        pcm->Velocity[0] = _bin.Velocity[0];
        pcm->Velocity[1] = _bin.Velocity[1];
        pcm->Velocity[2] = _bin.Velocity[2];
        pcm->PID  = - std::abs(_bin.PID);
    }

    //! correct orbit-samping/pseudo particles force
    /*!
      replace c.m. force by the averaged force on sample/pseudo particles
      @param[in,out] _ptcl_artificial: one group of artificial particles 
    */
    // template <class Tptcl>
    void correctOrbitalParticleForce(Particle* _ptcl_artificial) {
        auto* porb = getOrbitalParticles(_ptcl_artificial);
        const int n_orb = getOrbitalParticleN();

        if (n_orb>0) {
            auto* pcm = getCMParticles(_ptcl_artificial);
;
            REAL& acc_cm_0 = pcm->a_irr[0][0];  // First row, first column
            REAL& acc_cm_1 = pcm->a_irr[1][0];  // Second row, first column
            REAL& acc_cm_2 = pcm->a_irr[2][0];
        
            acc_cm_0=0.0;
            acc_cm_1=0.0;
            acc_cm_2=0.0;
            REAL m_ob_tot = 0.0;

            for (int k=0; k<n_orb; k++) {
                acc_cm_0 += porb[k].Mass*porb[k].a_irr[0][0]; 
                acc_cm_1 += porb[k].Mass*porb[k].a_irr[1][0]; 
                acc_cm_2 += porb[k].Mass*porb[k].a_irr[2][0]; 
                m_ob_tot += porb[k].Mass;
            }
            acc_cm_0 /= m_ob_tot;
            acc_cm_1 /= m_ob_tot;
            acc_cm_2 /= m_ob_tot;

#ifdef ARTIFICIAL_PARTICLE_DEBUG
            assert(abs(m_ob_tot-pcm->GroupInfo->artificial.getMassBackup())<1e-10);
#endif
        }
    }

    //! correct artificial particles force for furture use
    /*! 
      substract c.m. force (acc) from tidal tensor force (acc)\n
      replace c.m. force by the averaged force on orbital particles
      @param[in,out] _ptcl_artificial: one group of artificial particles 
    */
    // template <class Tptcl>
    void correctArtficialParticleForce(Particle* _ptcl_artificial) {
        auto* pcm = getCMParticles(_ptcl_artificial);
        // substract c.m. force (acc) from tidal tensor force (acc)
        auto* ptt = getTidalTensorParticles(_ptcl_artificial);
        TidalTensor::subtractCMForce(ptt, *pcm);

        // not consistent
        // After c.m. force used, it can be replaced by the averaged force on orbital particles
        // correctOrbitalParticleForce(_ptcl_artificial);
    }

    //! get oribit/pseudo particle list address from a artificial particle array
    // template <class Tptcl>
    Particle* getOrbitalParticles(Particle* _ptcl_list)  {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetOrb()].GroupInfo.artificial.isArtificial());
#endif
        return &_ptcl_list[getIndexOffsetOrb()];
    }

    //! get tidal tensor particle list address from a artificial particle array
    // template <class Tptcl>
    Particle* getTidalTensorParticles(Particle* _ptcl_list) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetTT()].GroupInfo.artificial.isArtificial());
#endif
        return &_ptcl_list[getIndexOffsetTT()];
    }

    //! get c.m. particle list address from a artificial particle array
    // template <class Tptcl>
    Particle* getCMParticles(Particle* _ptcl_list)  {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetCM()].GroupInfo.artificial.isCM());
#endif
        return &_ptcl_list[getIndexOffsetCM()];
    }

    //! tidal tensor particle index offset
    int getIndexOffsetTT() const {
        return 0;
    }

    //! orbital/pseudo particle index offset
    int getIndexOffsetOrb() const {
        return TidalTensor::getParticleN();
    }

    //! CM particle index offset
    int getIndexOffsetCM() const {
        return TidalTensor::getParticleN() + orbit_manager.getParticleN();
    }

    //! get artificial particle total number
    int getArtificialParticleN() const {
        return TidalTensor::getParticleN() + orbit_manager.getParticleN() + 1;
    }

    //! get artificial particle total number
    int getTidalTensorParticleN() const {
        return TidalTensor::getParticleN();
    }

    //! get orbitial particle number 
    int getOrbitalParticleN() const {
        return orbit_manager.getParticleN();
    }

    //! get left member number
    // template <class Tptcl>
    int getLeftMemberN(const Particle* _ptcl_list[]) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[0].GroupInfo.artificial.isArtificial());
#endif
        return int(_ptcl_list[0]->GroupInfo->artificial.getData(true));
    }

    //! get member number
    // template <class Tptcl>
    int getMemberN(const Particle* _ptcl_list[]) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetCM()].GroupInfo.artificial.isCM());
#endif
        return int(_ptcl_list[getIndexOffsetCM()]->GroupInfo->artificial.getData(true));
    }

    //! get right member number
    // template <class Tptcl>
    int getRightMemberN(const Particle* _ptcl_list[]) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[1].GroupInfo.artificial.isArtificial());
#endif
        return int(_ptcl_list[1]->GroupInfo->artificial.getData(true));
    }

    //! get center of mass id
    // template <class Tptcl>
    int getCMID(const Particle* _ptcl_list) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[getIndexOffsetCM()].GroupInfo.artificial.isCM());
#endif
        return -int(_ptcl_list[getIndexOffsetCM()].PID);
    }

    //! get stored data 
    // template <class Tptcl>
    REAL getStoredData(const Particle* _ptcl_list[], const int _index, const bool _is_positive) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_ptcl_list[_index+2].GroupInfo.artificial.isArtificial());
#endif
        return _ptcl_list[_index+2]->GroupInfo->artificial.getData(_is_positive);
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) {
        fwrite(&r_tidal_tensor, sizeof(REAL), 1, _fp);
        fwrite(&id_offset,      sizeof(int), 1, _fp);
        fwrite(&gravitational_constant, sizeof(REAL), 1, _fp);
        orbit_manager.writeBinary(_fp);
    }    

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
        size_t rcount = 0;
        rcount += fread(&r_tidal_tensor, sizeof(REAL), 1, _fin);
        rcount += fread(&id_offset,      sizeof(int), 1, _fin);
        rcount += fread(&gravitational_constant, sizeof(REAL), 1, _fin);
        if (rcount<3) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
        orbit_manager.readBinary(_fin);
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"r_tidal_tensor   : "<<r_tidal_tensor<<std::endl
             <<"id_offset        : "<<id_offset<<std::endl
             <<"G:               : "<<gravitational_constant<<std::endl;
        orbit_manager.print(_fout);
    }    


};
