#pragma once
#include "Common/binary_tree.h"
#include <iostream>
#include "../defs.h"
#include "Particle.h"

//! PseudoParticle Multipole method for counter force from binaries
class PseudoParticleMultipoleManager{
public:
    //! check paramters
    bool checkParams() {
        return true;
    }

    //! create pseudoparticle multipole particles
    /*! Create three particles that have the same quadrupole moment as orbit-averaged binaries.
        @param[in] _ptcl_artificial: particle array to store the sample particles, 2*n_split_ will be used
        @param[in] _bin: binary orbit 
     */
    // template <class Tptcl>
    void createSampleParticles(Particle* _ptcl_artificial,
                               COMM::BinaryTree<Particle,COMM::Binary> &_bin) {
        REAL m12 = _bin.Mass;
        REAL mu = _bin.m1*_bin.m2/m12;
        REAL prefactor = _bin.semi*std::sqrt(mu/m12);
        REAL pmass = m12/3.0;
        REAL ecc2 = _bin.ecc*_bin.ecc;

        REAL beta = prefactor*std::sqrt((1-ecc2)/4.0);
        REAL alpha= prefactor*std::sqrt(3*ecc2+0.75);

        Particle* pi = &(_ptcl_artificial[0]);
        pi->Mass = pmass;
        // pi->pos = PS::F64vec(0, 2*beta, 0 );
        pi->Position[0] = 0;
        pi->Position[1] = 2*beta;
        pi->Position[2] = 0;
        _bin.rotateToOriginalFrame(&(pi->Position[0]));
        for (int i; i<Dim; i++) {
            pi->Position[i] += _bin.Position[i];
            pi->Velocity[i] = _bin.Velocity[i];
        }
        
        pi = &(_ptcl_artificial[1]);
        pi->Mass = pmass;
        // pi->pos = PS::F64vec(alpha, -beta, 0 );
        pi->Position[0] = alpha;
        pi->Position[1] = -beta;
        pi->Position[2] = 0;
        _bin.rotateToOriginalFrame(&(pi->Position[0]));
        for (int i; i<Dim; i++) {
            pi->Position[i] += _bin.Position[i];
            pi->Velocity[i] = _bin.Velocity[i];
        }

        pi = &(_ptcl_artificial[2]);
        pi->Mass = pmass;
        // pi->pos = PS::F64vec(-alpha, -beta, 0 );
        pi->Position[0] = -alpha;
        pi->Position[1] = -beta;
        pi->Position[2] = 0;
        _bin.rotateToOriginalFrame(&(pi->Position[0]));
        for (int i; i<Dim; i++) {
            pi->Position[i] += _bin.Position[i];
            pi->Velocity[i] = _bin.Velocity[i];
        }

#ifdef ARTIFICIAL_PARTICLE_DEBUG
        PS::F64vec dv = (_ptcl_artificial[0].pos + _ptcl_artificial[1].pos + _ptcl_artificial[2].pos)/3 - _bin.pos;
        assert(dv*dv<1e-10);
#endif
    }

    //! get particle number 
    static int getParticleN() {
        return 3;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) {
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
    }
};
