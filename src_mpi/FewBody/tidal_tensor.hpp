#pragma once
#include <iostream>
#include <iomanip>
#include "../def.h"
#include "../particle.h"
#include "Group.h"

//! Tidal tensor perterbation for AR
class TidalTensor{
private:
    double T1[3];     // 0 constant 
    double T2[9];  // 1st (9)    general tensor
    //PS::F64 T2[6];  // 1st (6)  symmetry tensor
#ifdef TIDAL_TENSOR_3RD
    double T3[10]; // 2nd Tensor (10)
#endif
public:
    double Position[Dim];  // position of c.m.
    // REAL group_id; // indicate which group use the tensor // Eunwoo: or GroupInfo???
    Group* GroupInfo;


    TidalTensor(): T1{0.0, 0.0, 0.0}, 
                   T2{0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0}, 
#ifdef TIDAL_TENSOR_3RD
                   T3{0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.0,  0.0}, 
#endif
                   Position{0.0, 0.0, 0.0}, GroupInfo(nullptr) {}

    void dump(FILE *fp){
        fwrite(this, sizeof(TidalTensor),1,fp);
    }

    void read(FILE *fp){
        size_t rcount = fread(this, sizeof(TidalTensor),1,fp);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    void clear(){
        T1[0] = T1[1] = T1[2] = 0;
        for(int i=0; i<9; i++) T2[i] = 0;
#ifdef TIDAL_TENSOR_3RD
        for(int i=0; i<10; i++) T3[i] = 0;
#endif
        Position[0] = Position[1] = Position[2] = 0.0;
        GroupInfo = nullptr;
    }

    //! create tidal tensor measurement particles 
    /*! 
       2nd order: creat 4 zero-mass particles at the corners of regular tentrahedron with edge size of 0.16*_size. the cente is c.m. particle
       3rd order: creat 8 zero-mass particles at the corners of cube with edge size of 0.16*_size. the cente is c.m. particle
     */
    // template<class Tptcl>
    static void createTidalTensorMeasureParticles(Particle* _ptcl_tt, const Particle& _ptcl_cm, const double _size) {
#ifdef TIDAL_TENSOR_3RD
        ///* Assume _size is the maximum length inside box
        //   Then the edge length=_size/(2*sqrt(2))
        // */

        double lscale = 0.16*_size; // half edge size

        // set box 
        // _ptcl_tt[0].pos = PS::F64vec(lscale,   0,       -lscale) + _ptcl_cm.pos;
        // _ptcl_tt[1].pos = PS::F64vec(0,        lscale,  -lscale) + _ptcl_cm.pos;
        // _ptcl_tt[2].pos = PS::F64vec(-lscale,  0,       -lscale) + _ptcl_cm.pos;
        // _ptcl_tt[3].pos = PS::F64vec(0,       -lscale,  -lscale) + _ptcl_cm.pos;
        // _ptcl_tt[4].pos = PS::F64vec(lscale,   0,        lscale) + _ptcl_cm.pos;
        // _ptcl_tt[5].pos = PS::F64vec(0,        lscale,   lscale) + _ptcl_cm.pos;
        // _ptcl_tt[6].pos = PS::F64vec(-lscale,  0,        lscale) + _ptcl_cm.pos;
        // _ptcl_tt[7].pos = PS::F64vec(0,       -lscale,   lscale) + _ptcl_cm.pos;

        _ptcl_tt[0].Position[0] = lscale + _ptcl_cm.Position[0];
        _ptcl_tt[0].Position[1] = 0 + _ptcl_cm.Position[1];
        _ptcl_tt[0].Position[2] = -lscale + _ptcl_cm.Position[2];
        _ptcl_tt[1].Position[0] = 0 + _ptcl_cm.Position[0];
        _ptcl_tt[1].Position[1] = lscale + _ptcl_cm.Position[1];
        _ptcl_tt[1].Position[2] = -lscale + _ptcl_cm.Position[2];
        _ptcl_tt[2].Position[0] = -lscale + _ptcl_cm.Position[0];
        _ptcl_tt[2].Position[1] = 0 + _ptcl_cm.Position[1];
        _ptcl_tt[2].Position[2] = -lscale + _ptcl_cm.Position[2];
        _ptcl_tt[3].Position[0] = 0 + _ptcl_cm.Position[0];
        _ptcl_tt[3].Position[1] = -lscale + _ptcl_cm.Position[1];
        _ptcl_tt[3].Position[2] = -lscale + _ptcl_cm.Position[2];
        _ptcl_tt[4].Position[0] = lscale + _ptcl_cm.Position[0];
        _ptcl_tt[4].Position[1] = 0 + _ptcl_cm.Position[1];
        _ptcl_tt[4].Position[2] = lscale + _ptcl_cm.Position[2];
        _ptcl_tt[5].Position[0] = 0 + _ptcl_cm.Position[0];
        _ptcl_tt[5].Position[1] = lscale + _ptcl_cm.Position[1];
        _ptcl_tt[5].Position[2] = lscale + _ptcl_cm.Position[2];
        _ptcl_tt[6].Position[0] = -lscale + _ptcl_cm.Position[0];
        _ptcl_tt[6].Position[1] = 0 + _ptcl_cm.Position[1];
        _ptcl_tt[6].Position[2] = lscale + _ptcl_cm.Position[2];
        _ptcl_tt[7].Position[0] = 0 + _ptcl_cm.Position[0];
        _ptcl_tt[7].Position[1] = -lscale + _ptcl_cm.Position[1];
        _ptcl_tt[7].Position[2] = lscale + _ptcl_cm.Position[2];

#else
        ///* Assume _size is the maximum length, 
        //   Then the edge length=_size
        // */

        double lscale = 0.5*_size; // half edge size

        // set tetrahedron
        // _ptcl_tt[0].pos = PS::F64vec( lscale,  0,       -lscale*0.707106781186548) + _ptcl_cm.pos;
        // _ptcl_tt[1].pos = PS::F64vec(-lscale,  0,       -lscale*0.707106781186548) + _ptcl_cm.pos;
        // _ptcl_tt[2].pos = PS::F64vec(0,        lscale,   lscale*0.707106781186548) + _ptcl_cm.pos;
        // _ptcl_tt[3].pos = PS::F64vec(0,       -lscale,   lscale*0.707106781186548) + _ptcl_cm.pos;

        _ptcl_tt[0].Position[0] = lscale + _ptcl_cm.Position[0];
        _ptcl_tt[0].Position[1] = 0 + _ptcl_cm.Position[1];
        _ptcl_tt[0].Position[2] = -lscale*0.707106781186548 + _ptcl_cm.Position[2];
        _ptcl_tt[1].Position[0] = -lscale + _ptcl_cm.Position[0];
        _ptcl_tt[1].Position[1] = 0 + _ptcl_cm.Position[1];
        _ptcl_tt[1].Position[2] = -lscale*0.707106781186548 + _ptcl_cm.Position[2];
        _ptcl_tt[2].Position[0] = 0 + _ptcl_cm.Position[0];
        _ptcl_tt[2].Position[1] = lscale + _ptcl_cm.Position[1];
        _ptcl_tt[2].Position[2] = lscale*0.707106781186548 + _ptcl_cm.Position[2];
        _ptcl_tt[3].Position[0] = 0 + _ptcl_cm.Position[0];
        _ptcl_tt[3].Position[1] = -lscale + _ptcl_cm.Position[1];
        _ptcl_tt[3].Position[2] = lscale*0.707106781186548 + _ptcl_cm.Position[2];

#endif

        for (int i=0; i<getParticleN(); i++) {
            // co-moving velocity
            _ptcl_tt[i].Velocity[0]  = _ptcl_cm.Velocity[0];
            _ptcl_tt[i].Velocity[1]  = _ptcl_cm.Velocity[1];
            _ptcl_tt[i].Velocity[2]  = _ptcl_cm.Velocity[2];
            // no mass
            _ptcl_tt[i].Mass = 0.0;
        }
    }

    //! subtract c.m. force from measure points
    // template<class Tptcl>
    static void subtractCMForce(Particle* _ptcl_tt, const Particle& _ptcl_cm) {
        for (int k=0; k<getParticleN(); k++) {
            // _ptcl_tt[k].acc -= _ptcl_cm.acc;

            _ptcl_tt[k].a_irr[0][0] -= _ptcl_cm.a_irr[0][0];
            _ptcl_tt[k].a_irr[1][0] -= _ptcl_cm.a_irr[1][0];
            _ptcl_tt[k].a_irr[2][0] -= _ptcl_cm.a_irr[2][0];
        }
    }

    //! get particle number 
    static int getParticleN() {
#ifdef TIDAL_TENSOR_3RD
        return 8;
#else
        return 4;
#endif
    }

    //! tidal tensor fitting function,
    /*! 
       Symmetry T2:
       xx xy xz 0     1     2
       yx yy yz 3(1)  4     5
       zx zy zz 6(2)  7(5)  8

       General T2:
       xx xy xz 0  1  2
       yx yy yz 3  4  5
       zx zy zz 6  7  8

       Symmetry T3:
       xxx xxy xxz  0  1  2
       xyx xyy xyz  1  3  4
       xzx xzy xzz  2  4  5

       yxx yxy yxz  1  3  4
       yyx yyy yyz  3  6  7  
       yzx yzy yzz  4  7  8
      
       zxx zxy zxz  2  4  5
       zyx zyy zyz  4  7  8
       zzx zzy zzz  5  8  9
       
       @param[in] _ptcl_tt: tidal tensor measure particles
       @param[in] _ptcl_cm: tidal tensor measure particle c.m.
       @param[in] _size: particle box size
    */
    // template<class Tptcl>
    void fit(Particle* _ptcl_tt, Particle& _ptcl_cm,  const double _size) {
        // get c.m. position
        Position[0] = _ptcl_cm.Position[0];
        Position[1] = _ptcl_cm.Position[1];
        Position[2] = _ptcl_cm.Position[2];

        int n_point = getParticleN();

        double fi[n_point][Dim];

        // get acceleration
        for (int i=0; i<n_point; i++) {
            // fi[i] = _ptcl_tt[i].acc;
            fi[i][0] = _ptcl_tt[i].a_irr[0][0];
            fi[i][1] = _ptcl_tt[i].a_irr[1][0];
            fi[i][2] = _ptcl_tt[i].a_irr[2][0];
        }

        // get cofficients
        // T1, assume input force already remove the c.m.
        T1[0] = T1[1] = T1[2] = 0.0;
        
#ifdef TIDAL_TENSOR_3RD
        // T2, general form
        // 0 1 2
        T2[0] =  0.250000000000000*fi[0][0] + -0.250000000000000*fi[2][0] +  0.250000000000000*fi[4][0] + -0.250000000000000*fi[6][0];
        T2[1] =  0.125000000000000*fi[0][1] +  0.125000000000000*fi[1][0] + -0.125000000000000*fi[2][1] + -0.125000000000000*fi[3][0] 
            +    0.125000000000000*fi[4][1] +  0.125000000000000*fi[5][0] + -0.125000000000000*fi[6][1] + -0.125000000000000*fi[7][0];
        T2[2] = -0.083333333333333*fi[0][0] +  0.083333333333333*fi[0][2] + -0.083333333333333*fi[1][0] + -0.083333333333333*fi[2][0]
            +   -0.083333333333333*fi[2][2] + -0.083333333333333*fi[3][0] +  0.083333333333333*fi[4][0] +  0.083333333333333*fi[4][2]
            +    0.083333333333333*fi[5][0] +  0.083333333333333*fi[6][0] + -0.083333333333333*fi[6][2] +  0.083333333333333*fi[7][0];

        // 3 4 5
        T2[3] =  T2[1];
        T2[4] =  0.250000000000000*fi[1][1] + -0.250000000000000*fi[3][1] +  0.250000000000000*fi[5][1] + -0.250000000000000*fi[7][1];
        T2[5] = -0.083333333333333*fi[0][1] + -0.083333333333333*fi[1][1] +  0.083333333333333*fi[1][2] + -0.083333333333333*fi[2][1]
            +   -0.083333333333333*fi[3][1] + -0.083333333333333*fi[3][2] +  0.083333333333333*fi[4][1] +  0.083333333333333*fi[5][1]
            +    0.083333333333333*fi[5][2] +  0.083333333333333*fi[6][1] +  0.083333333333333*fi[7][1] + -0.083333333333333*fi[7][2];

        // 6 7 8
        T2[6] =  T2[2];
        T2[7] =  T2[5];
        T2[8] = -0.125000000000000*fi[0][2] + -0.125000000000000*fi[1][2] + -0.125000000000000*fi[2][2] + -0.125000000000000*fi[3][2]
            +    0.125000000000000*fi[4][2] +  0.125000000000000*fi[5][2] +  0.125000000000000*fi[6][2] +  0.125000000000000*fi[7][2];

        // T3, symmetry form
        T3[0] =  0.250000000000000*fi[0][0] +  0.125000000000000*fi[0][2] +  0.250000000000000*fi[2][0] + -0.125000000000000*fi[2][2]
            +    0.250000000000000*fi[4][0] + -0.125000000000000*fi[4][2] +  0.250000000000000*fi[6][0] +  0.125000000000000*fi[6][2];
        T3[1] =  0.250000000000000*fi[0][1] +  0.125000000000000*fi[1][2] +  0.250000000000000*fi[2][1] + -0.125000000000000*fi[3][2]
            +    0.250000000000000*fi[4][1] + -0.125000000000000*fi[5][2] +  0.250000000000000*fi[6][1] +  0.125000000000000*fi[7][2];
        T3[2] = -0.112500000000000*fi[0][0] +  0.025000000000000*fi[0][2] + -0.012500000000000*fi[1][1] + -0.025000000000000*fi[1][2]
            +    0.112500000000000*fi[2][0] +  0.025000000000000*fi[2][2] +  0.012500000000000*fi[3][1] + -0.025000000000000*fi[3][2]
            +    0.112500000000000*fi[4][0] +  0.025000000000000*fi[4][2] +  0.012500000000000*fi[5][1] + -0.025000000000000*fi[5][2]
            +   -0.112500000000000*fi[6][0] +  0.025000000000000*fi[6][2] + -0.012500000000000*fi[7][1] + -0.025000000000000*fi[7][2];
        T3[3] =  0.125000000000000*fi[0][2] +  0.250000000000000*fi[1][0] + -0.125000000000000*fi[2][2] +  0.250000000000000*fi[3][0]
            +   -0.125000000000000*fi[4][2] +  0.250000000000000*fi[5][0] +  0.125000000000000*fi[6][2] +  0.250000000000000*fi[7][0];
        T3[4] = -0.062500000000000*fi[0][1] + -0.062500000000000*fi[1][0] +  0.062500000000000*fi[2][1] +  0.062500000000000*fi[3][0]
            +    0.062500000000000*fi[4][1] +  0.062500000000000*fi[5][0] + -0.062500000000000*fi[6][1] + -0.062500000000000*fi[7][0];
        T3[5] = -0.125000000000000*fi[0][2] +  0.125000000000000*fi[2][2] +  0.125000000000000*fi[4][2] + -0.125000000000000*fi[6][2];
        T3[6] =  0.250000000000000*fi[1][1] +  0.125000000000000*fi[1][2] +  0.250000000000000*fi[3][1] + -0.125000000000000*fi[3][2]
            +    0.250000000000000*fi[5][1] + -0.125000000000000*fi[5][2] +  0.250000000000000*fi[7][1] +  0.125000000000000*fi[7][2];
        T3[7] = -0.012500000000000*fi[0][0] + -0.025000000000000*fi[0][2] + -0.112500000000000*fi[1][1] +  0.025000000000000*fi[1][2]
            +    0.012500000000000*fi[2][0] + -0.025000000000000*fi[2][2] +  0.112500000000000*fi[3][1] +  0.025000000000000*fi[3][2]
            +    0.012500000000000*fi[4][0] + -0.025000000000000*fi[4][2] +  0.112500000000000*fi[5][1] +  0.025000000000000*fi[5][2]
            +   -0.012500000000000*fi[6][0] + -0.025000000000000*fi[6][2] + -0.112500000000000*fi[7][1] +  0.025000000000000*fi[7][2];
        T3[8] = -0.125000000000000*fi[1][2] +  0.125000000000000*fi[3][2] +  0.125000000000000*fi[5][2] + -0.125000000000000*fi[7][2];
        T3[9] =  0.062500000000000*fi[0][0] +  0.125000000000000*fi[0][2] +  0.062500000000000*fi[1][1] +  0.125000000000000*fi[1][2]
            +   -0.062500000000000*fi[2][0] +  0.125000000000000*fi[2][2] + -0.062500000000000*fi[3][1] +  0.125000000000000*fi[3][2]
            +   -0.062500000000000*fi[4][0] +  0.125000000000000*fi[4][2] + -0.062500000000000*fi[5][1] +  0.125000000000000*fi[5][2]
            +    0.062500000000000*fi[6][0] +  0.125000000000000*fi[6][2] +  0.062500000000000*fi[7][1] +  0.125000000000000*fi[7][2];

        // Rescale
        //PS::F64 T2S = 1.0/(_bin.semi*(1+_bin.ecc)*0.35);
        double T2S = 1.0/(_size*0.16);
        double T3S = T2S*T2S;
        for (int i=0; i<9;  i++) T2[i] *= T2S;
        for (int i=0; i<10; i++) T3[i] *= T3S;

#else

        // T2, general form
        // 0 1 2
        T2[0] =  0.500000000000000*fi[0][0]+   -0.500000000000000*fi[1][0];
        T2[1] =  0.250000000000000*fi[0][1]+   -0.250000000000000*fi[1][1]+    0.250000000000000*fi[2][0]+   -0.250000000000000*fi[3][0];
        T2[2] = -0.176776695296637*fi[0][0]+    0.250000000000000*fi[0][2]+   -0.176776695296637*fi[1][0]+   -0.250000000000000*fi[1][2]+    0.176776695296637*fi[2][0]+    0.176776695296637*fi[3][0];
        
        // 3 4 5
        T2[3] = T2[1];
        T2[4] =  0.500000000000000*fi[2][1]+   -0.500000000000000*fi[3][1];
        T2[5] = -0.176776695296637*fi[0][1]+   -0.176776695296637*fi[1][1]+    0.176776695296637*fi[2][1]+    0.250000000000000*fi[2][2]+    0.176776695296637*fi[3][1]+   -0.250000000000000*fi[3][2];

        // 6 7 8
        T2[6] =  T2[2];
        T2[7] =  T2[5];
        T2[8] = -0.353553390593274*fi[0][2]+   -0.353553390593274*fi[1][2]+    0.353553390593274*fi[2][2]+    0.353553390593274*fi[3][2];

        // Rescale
        double T2S = 2.0/_size;
        for (int i=0; i<9;  i++) T2[i] *= T2S;
#endif
    }

    //! Shift c.m. to new reference position
    /*! Only the 1st order tensor need a correction from 2nd order 
      
      ### 1st order:
      T2:
      xx xy xz 0  1  2
      yx yy yz 3  4  5
      zx zy zz 6  7  8

      ### 2nd order:
      xxx xxy xxz  0  1  2
      xyx xyy xyz  1  3  4
      xzx xzy xzz  2  4  5

      yxx yxy yxz  1  3  4
      yyx yyy yyz  3  6  7  
      yzx yzy yzz  4  7  8
      
      zxx zxy zxz  2  4  5
      zyx zyy zyz  4  7  8
      zzx zzy zzz  5  8  9
      
      @param[in] _pos: new c.m. position
     */
    void shiftCM(const double (&_pos)[3]) {
#ifdef TIDAL_TENSOR_3RD
        double dr[0] = _pos[0]-Position[0];
        double dr[1] = _pos[1]-Position[1];
        double dr[2] = _pos[2]-Position[2];

        double x = dr[0];
        double y = dr[1];
        double z = dr[2];
        //PS::F64 x2 = x*x;
        //PS::F64 xy = x*y;
        //PS::F64 xz = x*z;
        //PS::F64 y2 = y*y;
        //PS::F64 yz = y*z;
        //PS::F64 z2 = z*z;

        // T1 += T2^dr + dr^T3^dr
        /*! c.m. force need to be removed, otherwise the perturbation force sum of all particles are not zero 
        T1[0] += T2[0]*x + T2[1]*y + T2[2]*z 
            +    T3[0]*x2 + 2*T3[1]*xy + 2*T3[2]*xz + T3[3]*y2 + 2*T3[4]*yz + T3[5]*z2;
        T1[1] += T2[3]*x + T2[4]*y + T2[5]*z
            +    T3[1]*x2 + 2*T3[3]*xy + 2*T3[4]*xz + T3[6]*y2 + 2*T3[7]*yz + T3[8]*z2;
        T1[2] += T2[6]*x + T2[7]*y + T2[8]*z
            +    T3[2]*x2 + 2*T3[4]*xy + 2*T3[5]*xz + T3[7]*y2 + 2*T3[8]*yz + T3[9]*z2;
        */
        
        // T2 += 2*dr^T3
        T2[0] += 2.0*(T3[0]*x + T3[1]*y + T3[2]*z); // xx: xxx*x + xyx*y + xzx*z
        T2[1] += 2.0*(T3[1]*x + T3[3]*y + T3[4]*z); // xy: xxy*x + xyy*y + xzy*z
        T2[2] += 2.0*(T3[2]*x + T3[4]*y + T3[5]*z); // xy: xxz*x + xyz*y + xzz*z

        T2[3] += 2.0*(T3[1]*x + T3[3]*y + T3[4]*z); // yx: yxx*x + yyx*y + yzx*z
        T2[4] += 2.0*(T3[3]*x + T3[6]*y + T3[7]*z); // yy: yxy*x + yyy*y + yzy*z
        T2[5] += 2.0*(T3[4]*x + T3[7]*y + T3[8]*z); // yy: yxz*x + yyz*y + yzz*z

        T2[6] += 2.0*(T3[2]*x + T3[4]*y + T3[5]*z); // zx: zxx*x + zyx*y + zzx*z
        T2[7] += 2.0*(T3[4]*x + T3[7]*y + T3[8]*z); // zy: zxy*x + zyy*y + zzy*z
        T2[8] += 2.0*(T3[5]*x + T3[8]*y + T3[9]*z); // zy: zxz*x + zyz*y + zzz*z

#endif
        // update c.m.
        Position[0] = _pos[0];
        Position[1] = _pos[1];
        Position[2] = _pos[2];
    }

    void eval(double (&acc)[3], const double (&pos)[3]) const {
        double acc0=acc[0];
        double acc1=acc[1];
        double acc2=acc[2];

        double x = pos[0];
        double y = pos[1];
        double z = pos[2];

#ifdef TIDAL_TENSOR_3RD
        /*
          T2:
          [[0 1 2]
           [3 4 5]
           [6 7 8]]

          T3:
          [[[6  7  8 ]
            [7  9  10]
            [8  10 11]]

           [[7  9  10]
            [9  12 13]
            [10 13 14]]

           [[8  10 11]
            [10 13 14]
            [11 14 15]]]

        */
        double x2 = x*x;
        double xy = x*y;
        double xz = x*z;
        double y2 = y*y;
        double yz = y*z;
        double z2 = z*z;

        acc0 +=  T1[0] + T2[0]*x + T2[1]*y + T2[2]*z 
            +      T3[0]*x2 + 2*T3[1]*xy + 2*T3[2]*xz + T3[3]*y2 + 2*T3[4]*yz + T3[5]*z2;
        acc1 +=  T1[1] + T2[3]*x + T2[4]*y + T2[5]*z
            +      T3[1]*x2 + 2*T3[3]*xy + 2*T3[4]*xz + T3[6]*y2 + 2*T3[7]*yz + T3[8]*z2;
        acc2 +=  T1[2] + T2[6]*x + T2[7]*y + T2[8]*z
            +      T3[2]*x2 + 2*T3[4]*xy + 2*T3[5]*xz + T3[7]*y2 + 2*T3[8]*yz + T3[9]*z2;

#else
        acc0 +=  T1[0] + T2[0]*x + T2[1]*y + T2[2]*z;
        acc1 +=  T1[1] + T2[3]*x + T2[4]*y + T2[5]*z;
        acc2 +=  T1[2] + T2[6]*x + T2[7]*y + T2[8]*z;
#endif

        acc[0] = acc0;
        acc[1] = acc1;
        acc[2] = acc2;
    }

    double evalPot(const double (&pos)[3]) const {
        double x = pos[0];
        double y = pos[1];
        double z = pos[2];
#ifdef TIDAL_TENSOR_3RD
        double x2 = x*x;
        double xy = x*y;
        double xz = x*z;
        double y2 = y*y;
        double yz = y*z;
        double z2 = z*z;

        double acc0 =  T1[0] + 0.5*(T2[0]*x + T2[1]*y + T2[2]*z) 
            +      (T3[0]*x2 + 2*T3[1]*xy + 2*T3[2]*xz + T3[3]*y2 + 2*T3[4]*yz + T3[5]*z2)/3.0;
        double acc1 =  T1[1] + 0.5*(T2[3]*x + T2[4]*y + T2[5]*z)
            +      (T3[1]*x2 + 2*T3[3]*xy + 2*T3[4]*xz + T3[6]*y2 + 2*T3[7]*yz + T3[8]*z2)/3.0;
        double acc2 =  T1[2] + 0.5*(T2[6]*x + T2[7]*y + T2[8]*z)
            +      (T3[2]*x2 + 2*T3[4]*xy + 2*T3[5]*xz + T3[7]*y2 + 2*T3[8]*yz + T3[9]*z2)/3.0;
#else
        double acc0 =  T1[0] + 0.5*(T2[0]*x + T2[1]*y + T2[2]*z);
        double acc1 =  T1[1] + 0.5*(T2[3]*x + T2[4]*y + T2[5]*z);
        double acc2 =  T1[2] + 0.5*(T2[6]*x + T2[7]*y + T2[8]*z);
#endif        

        return - x*acc0 - y*acc1 - z*acc2;
    }

    void print(std::ostream & _fout, const int _width) const{
        _fout<<"T1: \n"
             <<std::setw(_width)<<T1[0]
             <<std::setw(_width)<<T1[1]
             <<std::setw(_width)<<T1[2]
             <<std::endl
             <<"T2: \n"
             <<std::setw(_width)<<T2[0]
             <<std::setw(_width)<<T2[1]
             <<std::setw(_width)<<T2[2]
             <<std::endl
             <<std::setw(_width)<<T2[3]
             <<std::setw(_width)<<T2[4]
             <<std::setw(_width)<<T2[5]
             <<std::endl
             <<std::setw(_width)<<T2[6]
             <<std::setw(_width)<<T2[7]
             <<std::setw(_width)<<T2[8]
#ifdef TIDAL_TENSOR_3RD
             <<std::endl
             <<"T3: \n"
             <<"x: \n"
             <<std::setw(_width)<<T3[0]
             <<std::setw(_width)<<T3[1]
             <<std::setw(_width)<<T3[2]
             <<std::endl
             <<std::setw(_width)<<T3[1]
             <<std::setw(_width)<<T3[3]
             <<std::setw(_width)<<T3[4]
             <<std::endl
             <<std::setw(_width)<<T3[2]
             <<std::setw(_width)<<T3[4]
             <<std::setw(_width)<<T3[5]
             <<std::endl
             <<"y: \n"
             <<std::setw(_width)<<T3[1]
             <<std::setw(_width)<<T3[3]
             <<std::setw(_width)<<T3[4]
             <<std::endl
             <<std::setw(_width)<<T3[3]
             <<std::setw(_width)<<T3[6]
             <<std::setw(_width)<<T3[7]
             <<std::endl
             <<std::setw(_width)<<T3[4]
             <<std::setw(_width)<<T3[7]
             <<std::setw(_width)<<T3[8]
             <<std::endl
             <<"x: \n"
             <<std::setw(_width)<<T3[2]
             <<std::setw(_width)<<T3[4]
             <<std::setw(_width)<<T3[5]
             <<std::endl
             <<std::setw(_width)<<T3[4]
             <<std::setw(_width)<<T3[7]
             <<std::setw(_width)<<T3[8]
             <<std::endl
             <<std::setw(_width)<<T3[5]
             <<std::setw(_width)<<T3[8]
             <<std::setw(_width)<<T3[9]
#endif
             <<std::endl;
    }
};