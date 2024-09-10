#pragma once

#include "Common/list.h"
// #include "Hermite/hermite_particle.h"
// #include "hard_ptcl.hpp"
#include "Common/binary_tree.h"
#include "tidal_tensor.hpp"
// #include "information.h"
#include "Particle.h"
#include "FB_defs.h"
// #include "Group.h"


//! Perturber class for AR integration
class Perturber {
public:

    TidalTensor* soft_pert;  ///> soft perturbation 
    Float soft_pert_min; ///> minimum soft perturbation

    Perturber(): soft_pert(NULL), soft_pert_min(Float(0.0)) {}

    //! clear function
    void clear() {
        if (soft_pert!=NULL) {
            assert(soft_pert->GroupInfo != nullptr);
            soft_pert->GroupInfo = nullptr;
            soft_pert = NULL;
        }
    }

    //! check parameters status
    bool checkParams() {
        //assert(soft_pert_min>=0.0);
        //assert(soft_pert!=NULL);
        return true;
    }

    //! calculate soft_pert_min for slowdown pert_out
    /*! \Delta F = G m_cm m_p (apo) / rp^3
        Pert_out = \Delta F /(G apo)
        @param[in] _bin: binary parameters
        @param[in] _G: gravitatioal constant
     */
    // template <class Tptcl>
    void calcSoftPertMin(const AR::BinaryTree<Particle>& _bin, const Float _G) {
        soft_pert_min = 0.0;
#ifdef SOFT_PERT
        if(soft_pert!=NULL&&_bin.semi>0.0) {
            Particle p[2];
            _bin.calcParticlesEcca(p[0], p[1], COMM::PI, _G);
            Float dacc_soft = 0.0;
            Float acc_p1[3] = {0.0, 0.0, 0.0};
            Float acc_p2[3] = {0.0, 0.0, 0.0};
            soft_pert->eval(acc_p1, p[0].Position);
            soft_pert->eval(acc_p2, p[1].Position);
            Float dacc[3] = {acc_p1[0]-acc_p2[0], 
                             acc_p1[1]-acc_p2[1],
                             acc_p1[2]-acc_p2[2]};
            dacc_soft = std::sqrt(dacc[0]*dacc[0] + dacc[1]*dacc[1] + dacc[2]*dacc[2]);
            Float apo = _bin.semi*(1.0+_bin.ecc);
            soft_pert_min = _bin.Mass*dacc_soft/(_G*apo);
        }
#endif
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumnTitle(std::ostream & _fout, const int _width=20) {
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
    */
    void printColumn(std::ostream & _fout, const int _width=20){
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) const {
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
    } 

};
