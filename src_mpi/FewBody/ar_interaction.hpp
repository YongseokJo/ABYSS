#pragma once
#include <cmath>

#include "FB_defs.h"

extern REAL EnzoTimeStep;
extern FILE* mergerout;

#include "Common/Float.h"
#include "Common/binary_tree.h"
#include "AR/force.h"
#include "ar_perturber.hpp"
#include <cassert>


#define ASSERT(x) assert(x)

// #define SEVN
#ifdef SEVN
extern FILE* SEVNout;
// extern IO* sevnio;
void UpdateEvolution(Particle* ptcl);
void Mix(Star* star1, Star* star2);
#endif


//! AR interaction clas
class Interaction{
public:

    Float gravitational_constant;
    // TwoBodyTide tide;

    Interaction(): gravitational_constant(Float(1.0)) {}


    //! (Necessary) check whether publicly initialized parameters are correctly set
    /*! \return true: all parmeters are correct. In this case no parameters, return true;
     */
    bool checkParams() {
        ASSERT(gravitational_constant>0.0);
        return true;
    }        

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<"G      : "<<gravitational_constant<<std::endl;
    }    

    //! (Necessary) calculate inner member acceleration, potential and inverse time transformation function gradient and factor for kick (two-body case)
    /*!
      @param[out] _f1: force for particle 1 to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard are overwritten, not accummulating old values)
      @param[out] _f2: force for particle 2
      @param[out] _epot: total inner potential energy
      @param[in] _p1: particle 1
      @param[in] _p2: particle 2
      \return the inverse time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInvTwo(AR::Force& _f1, AR::Force& _f2, Float& _epot, const Particle& _p1, const Particle& _p2) {
        // acceleration
        const Float mass1 = _p1.Mass;
        const Float* pos1 = _p1.Position;

        const Float mass2 = _p2.Mass;
        const Float* pos2 = _p2.Position;

        Float gm1 = gravitational_constant*mass1;
        Float gm2 = gravitational_constant*mass2;
        Float gm1m2 = gm1*mass2;
        
        Float dr[3] = {pos2[0] -pos1[0],
                       pos2[1] -pos1[1],
                       pos2[2] -pos1[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float inv_r = 1.0/sqrt(r2);
        Float inv_r3 = inv_r*inv_r*inv_r;

        Float* acc1 = _f1.acc_in;
        Float* acc2 = _f2.acc_in;


        Float gmor3_1 = gm2*inv_r3;
        Float gmor3_2 = gm1*inv_r3;

        Float gm1or =  gm1*inv_r;
        Float gm2or =  gm2*inv_r;
        Float gm1m2or = gm1m2*inv_r;


        acc1[0] = gmor3_1 * dr[0];
        acc1[1] = gmor3_1 * dr[1];
        acc1[2] = gmor3_1 * dr[2];

        _f1.pot_in = -gm2or;

        acc2[0] = - gmor3_2 * dr[0];
        acc2[1] = - gmor3_2 * dr[1];
        acc2[2] = - gmor3_2 * dr[2];

        _f2.pot_in = -gm1or;


#ifdef AR_TTL 

        Float gm1m2or3 = gm1m2*inv_r3;
        Float* gtgrad1 = _f1.gtgrad;
        Float* gtgrad2 = _f2.gtgrad;
        gtgrad1[0] = gm1m2or3 * dr[0];
        gtgrad1[1] = gm1m2or3 * dr[1];
        gtgrad1[2] = gm1m2or3 * dr[2];

        gtgrad2[0] = - gtgrad1[0];
        gtgrad2[1] = - gtgrad1[1];
        gtgrad2[2] = - gtgrad1[2];
#endif

        // potential energy
        _epot = - gm1m2or;

        // transformation factor for kick
        Float gt_kick_inv = gm1m2or;
        // fprintf(stderr, "gt_kick_inv: %e\n", gt_kick_inv); // Eunwoo debug

        return gt_kick_inv;
    }

    //! calculate inner member acceleration, potential and inverse time transformation function gradient and factor for kick
    /*!
      @param[out] _force: force array to store the calculation results (in acc_in[3] for acceleration and gtgrad[3] for gradient, notice acc/gtgard may need to reset zero to avoid accummulating old values)
      @param[out] _epot: total inner potential energy
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      \return the inverse time transformation factor (gt_kick_inv) for kick step
    */
    inline Float calcInnerAccPotAndGTKickInv(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle) {
        _epot = Float(0.0);
        Float gt_kick_inv = Float(0.0);

        for (int i=0; i<_n_particle; i++) {
            const Float massi = _particles[i].Mass;
            const Float* posi = _particles[i].Position;
            Float* acci = _force[i].acc_in;
            acci[0] = acci[1] = acci[2] = Float(0.0);

#ifdef AR_TTL 
            Float* gtgradi = _force[i].gtgrad;
            gtgradi[0] = gtgradi[1] = gtgradi[2] = Float(0.0);
#endif

            Float poti = Float(0.0);
            Float gtki = Float(0.0);

            for (int j=0; j<_n_particle; j++) {
                if (i==j) continue;
                const Float massj = _particles[j].Mass;
                const Float* posj = _particles[j].Position; 
                Float dr[3] = {posj[0] -posi[0],
                               posj[1] -posi[1],
                               posj[2] -posi[2]};
                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float inv_r = 1.0/sqrt(r2);
                Float inv_r3 = inv_r*inv_r*inv_r;

                Float gmor3 = gravitational_constant*massj*inv_r3;
                Float gmor = gravitational_constant*massj*inv_r;

                acci[0] += gmor3 * dr[0];
                acci[1] += gmor3 * dr[1];
                acci[2] += gmor3 * dr[2];

#ifdef AR_TTL                     
                Float gmimjor3 = massi*gmor3;
                gtgradi[0] += gmimjor3 * dr[0];
                gtgradi[1] += gmimjor3 * dr[1];
                gtgradi[2] += gmimjor3 * dr[2];
#endif

                poti -= gmor;
                gtki += gmor;
                    
            }
            _epot += poti * massi;
            gt_kick_inv += gtki * massi;
        }
        _epot   *= 0.5;
        gt_kick_inv *= 0.5;
        // fprintf(stderr, "gt_kick_inv: %e\n", gt_kick_inv); // Eunwoo debug

        return gt_kick_inv;
    }

    //! (Necessary) calculate acceleration from perturber and the perturbation factor for slowdown calculation
    /*!@param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
    */
    void calcAccPert(AR::Force* _force, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber, const Float _time) {
        static const Float inv3 = 1.0 / 3.0;

        // perturber force
        // const int n_pert = _perturber.neighbor_address.getSize();
        const int n_pert = _particle_cm.NumberOfAC;
        // const int n_pert_single = _perturber.n_neighbor_single;
        // const int n_pert_group = _perturber.n_neighbor_group;

        if (n_pert>0) {

            Float time = _time;

            // auto* pert_adr = _perturber.neighbor_address.getDataAddress();
            auto pert_adr = _particle_cm.ACList;

            Float xp[n_pert][3], xcm[3], m[n_pert];
            // ChangeOver* changeover[n_pert_single];
            // H4::NBAdr<Particle>::Group* ptclgroup[n_pert_group];

            // int n_single_count=0;
            // int n_group_count=0;
            for (int j=0; j<n_pert; j++) {
                // H4::NBAdr<Particle>::Single* pertj;
                Particle* pertj;
                pertj = pert_adr[j];
                // int k; // index of predicted data
                // if (pert_adr[j].type==H4::NBType::group) {
                //     pertj = &(((H4::NBAdr<Particle>::Group*)pert_adr[j].adr)->cm);
                //     k = n_group_count + n_pert_single;
                //     ptclgroup[n_group_count] = (H4::NBAdr<Particle>::Group*)pert_adr[j].adr;
                //     n_group_count++;
                // }
                // else {
                //     pertj = (H4::NBAdr<Particle>::Single*)pert_adr[j].adr;
                //     k = n_single_count;
                //     changeover[n_single_count] = &pertj->changeover;
                //     n_single_count++;
                // }

                Float dt = time - pertj->CurrentTimeIrr*EnzoTimeStep;
                // ASSERT(dt>=0.0); // Eunwoo debug // Is this right?
                //ASSERT(dt>=-1e-7);
                xp[j][0] = pertj->Position[0] + dt*(pertj->Velocity[0] + 0.5*dt*(pertj->a_irr[0][0] + inv3*dt*pertj->a_irr[0][1]));
                xp[j][1] = pertj->Position[1] + dt*(pertj->Velocity[1] + 0.5*dt*(pertj->a_irr[1][0] + inv3*dt*pertj->a_irr[1][1]));
                xp[j][2] = pertj->Position[2] + dt*(pertj->Velocity[2] + 0.5*dt*(pertj->a_irr[2][0] + inv3*dt*pertj->a_irr[2][1]));


                m[j] = pertj->Mass;
            }
            // ASSERT(n_single_count == n_pert_single);
            // ASSERT(n_group_count == n_pert_group);

            Float dt = time - _particle_cm.CurrentTimeIrr*EnzoTimeStep;
            // ASSERT(dt>=0.0); // Eunwoo debug // Is this right?

            xcm[0] = _particle_cm.Position[0] + dt*(_particle_cm.Velocity[0] + 0.5*dt*(_particle_cm.a_irr[0][0] + inv3*dt*_particle_cm.a_irr[0][1]));
            xcm[1] = _particle_cm.Position[1] + dt*(_particle_cm.Velocity[1] + 0.5*dt*(_particle_cm.a_irr[1][0] + inv3*dt*_particle_cm.a_irr[1][1]));
            xcm[2] = _particle_cm.Position[2] + dt*(_particle_cm.Velocity[2] + 0.5*dt*(_particle_cm.a_irr[2][0] + inv3*dt*_particle_cm.a_irr[2][1]));


            Float acc_pert_cm[3]={0.0, 0.0, 0.0};
            Float mcm = 0.0;
            // if (_perturber.need_resolve_flag) {
            // calculate component perturbation
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                Float& pot_pert = _force[i].pot_pert;
                const auto& pi = _particles[i];
                // auto& chi = pi.changeover;
                acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
                pot_pert = 0.0;

                Float xi[3];
                xi[0] = pi.Position[0] + xcm[0];
                xi[1] = pi.Position[1] + xcm[1];
                xi[2] = pi.Position[2] + xcm[2];

                // single perturber
                for (int j=0; j<n_pert; j++) {
                    Float dr[3] = {xp[j][0] - xi[0],
                                   xp[j][1] - xi[1],
                                   xp[j][2] - xi[2]};
                    Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                    Float r  = sqrt(r2);
                    Float r3 = r*r2;
                    Float gm = gravitational_constant*m[j];
                    Float gmor3 = gm/r3;

                    acc_pert[0] += gmor3 * dr[0];
                    acc_pert[1] += gmor3 * dr[1];
                    acc_pert[2] += gmor3 * dr[2];

                }

                acc_pert_cm[0] += pi.Mass *acc_pert[0];
                acc_pert_cm[1] += pi.Mass *acc_pert[1];
                acc_pert_cm[2] += pi.Mass *acc_pert[2];

                mcm += pi.Mass;

            }
//#ifdef AR_DEBUG
//            ASSERT(abs(mcm-_particle_cm.mass)<1e-10);
//#endif
                
            // get cm perturbation (exclude soft pert)
            acc_pert_cm[0] /= mcm;
            acc_pert_cm[1] /= mcm;
            acc_pert_cm[2] /= mcm;

            // remove cm. perturbation
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                Float& pot_pert = _force[i].pot_pert;
                const auto& pi = _particles[i];
                acc_pert[0] -= acc_pert_cm[0]; 
                acc_pert[1] -= acc_pert_cm[1];        
                acc_pert[2] -= acc_pert_cm[2]; 
                
                pot_pert -= acc_pert[0]*pi.Position[0] + acc_pert[1]*pi.Position[1] + acc_pert[2]*pi.Position[2];

            }

        }
        else {
            for (int i=0; i<_n_particle; i++) {
                Float* acc_pert = _force[i].acc_pert;
                acc_pert[0] = acc_pert[1] = acc_pert[2] = Float(0.0);
                _force[i].pot_pert = 0.0;
            }  
        }
    }
    
    //! (Necessary) calculate acceleration from internal members and perturbers
    /*! The Force class acc_pert should be updated
      @param[out] _force: force array to store the calculation results (in acc_pert[3], notice acc_pert may need to reset zero to avoid accummulating old values)
      @param[out] _epot: potential 
      @param[in] _particles: member particle array
      @param[in] _n_particle: number of member particles
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
      @param[in] _time: current time
      \return perturbation energy to calculate slowdown factor
    */
    Float calcAccPotAndGTKickInv(AR::Force* _force, Float& _epot, const Particle* _particles, const int _n_particle, const Particle& _particle_cm, const Perturber& _perturber, const Float _time) {
        // inner force
        Float gt_kick_inv;
        if (_n_particle==2) gt_kick_inv = calcInnerAccPotAndGTKickInvTwo(_force[0], _force[1], _epot, _particles[0], _particles[1]);
        else gt_kick_inv = calcInnerAccPotAndGTKickInv(_force, _epot, _particles, _n_particle);

        calcAccPert(_force, _particles, _n_particle, _particle_cm, _perturber, _time);

        return gt_kick_inv;
    }


    //! calculate perturbation from c.m. acceleration
    Float calcPertFromForce(const Float* _force, const Float _mp, const Float _mpert) {
        Float force2 = _force[0]*_force[0]+_force[1]*_force[1]+_force[2]*_force[2];
#ifdef AR_SLOWDOWN_PERT_R4
        return force2/(gravitational_constant*_mp*_mpert);
#else
        Float force = sqrt(force2)/gravitational_constant;
        return sqrt(force/(_mp*_mpert))*force;
#endif
    }

    //! calculate perturbation from binary tree
    static Float calcPertFromBinary(const COMM::Binary& _bin) {
        Float apo = _bin.semi*(1.0+_bin.ecc);
        Float apo2 = apo*apo;
#ifdef AR_SLOWDOWN_PERT_R4
        return (_bin.m1*_bin.m2)/(apo2*apo2);
#else
        return (_bin.m1*_bin.m2)/(apo2*apo);
#endif
    }

    //! calculate perturbation from distance to perturber and masses of particle and perturber 
    static Float calcPertFromMR(const Float _r, const Float _mp, const Float _mpert) {
        Float r2 = _r*_r;
#ifdef AR_SLOWDOWN_PERT_R4
        return _mp*_mpert/(r2*r2);
#else
        return (_mp*_mpert)/(r2*_r);
#endif
    }

#if (defined AR_SLOWDOWN_ARRAY) || (defined AR_SLOWDOWN_TREE)

    //! calculate slowdown timescale
    void calcSlowDownTimeScale(Float& _t_min_sq, const Float dv[3], const Float dr[3], const Float& r, const Float& gm) {

        Float r2 = r*r;
        Float v2 = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
        Float drdv = dr[0]*dv[0] + dr[1]*dv[1] + dr[2]*dv[2];

        Float semi = 1.0/(2.0/r - v2/gm);
        //hyperbolic, directly use velocity v
        if (semi<0) 
            _t_min_sq = std::min(_t_min_sq, r2/v2);
        else {
            Float ra_fact = (1 - r/semi); 
            Float e2 = drdv*drdv/(gm*semi) + ra_fact*ra_fact; // ecc^2
            Float r_vrmax = semi*(1-e2);
            if (r<r_vrmax) {
                // avoid decrese of vr once the orbit pass, calculate vr max at cos(E)=e (r==semi*(1-e^2))
                // vr_max = sqrt(er*(drdv^2*er + r*vcr2^2))/(G(m1+m2)r)
                //        = e*sqrt[G(m1+m2)/(a*(1-e^2)]
                Float vrmax_sq = e2*gm/r_vrmax;
                //Float rv2 = r*v2;
                //Float er = 2*gm - rv2;
                //Float vcr2 = gm - rv2;
                //Float vrmax_sq = er*(drdv*drdv*er + r*vcr2*vcr2)/(gm*gm*r2);
                _t_min_sq = std::min(_t_min_sq, semi*semi/vrmax_sq);
            }
            else {
                // r/vr
                Float rovr = r2/abs(drdv);
                _t_min_sq = std::min(_t_min_sq, rovr*rovr);
            }
        }
    }

    //! calculate slowdown perturbation and timescale from particle j to particle i
    /*! 
      @param[out] _pert_out: perturbation from particle j
      @param[out] _t_min_sq: timescale limit from particle j
      @param[in] _pi: particle i (cm of binary)
      @param[in] _pj: particle j 
     */
    void calcSlowDownPertOne(Float& _pert_out, Float& _t_min_sq, const Particle& pi, const Particle& pj) {
        Float dr[3] = {pj.Position[0] - pi.Position[0],
                       pj.Position[1] - pi.Position[1],
                       pj.Position[2] - pi.Position[2]};
        Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
        Float r = sqrt(r2);
        _pert_out += calcPertFromMR(r, pi.Mass, pj.Mass);
            
#ifdef AR_SLOWDOWN_TIMESCALE
        Float dv[3] = {pj.Velocity[0] - pi.Velocity[0],
                       pj.Velocity[1] - pi.Velocity[1],
                       pj.Velocity[2] - pi.Velocity[2]};

        // identify whether hyperbolic or closed orbit
        Float gm = gravitational_constant*(pi.Mass+pj.Mass);

        calcSlowDownTimeScale(_t_min_sq, dv, dr, r, gm);
        // force dependent method
        // min sqrt(r^3/(G m))
        //Float gmor3 = (mp+mcm)*r*r2/(sdt->G*mp*mcm);
        //sdt->trf2_min =  std::min(sdt->trf2_min, gmor3);
#endif
    }

    //! (Necessary) calculate slowdown perturbation and timescale
    /*!
      @param[out] _pert_out: perturbation 
      @param[out] _t_min_sq: timescale limit 
      @param[in] _time: physical time for prediction
      @param[in] _particle_cm: center-of-mass particle
      @param[in] _perturber: pertuber container
    */
    void calcSlowDownPert(Float& _pert_out, Float& _t_min_sq, const Float& _time, const Particle& _particle_cm, const Perturber& _perturber) {
        static const Float inv3 = 1.0 / 3.0;

        // const int n_pert = _perturber.neighbor_address.getSize();
        const int n_pert = _particle_cm.NumberOfAC;

        if (n_pert>0) {

            auto pert_adr = _particle_cm.ACList;

            Float xp[3], xcm[3];
            Float dt = _time - _particle_cm.CurrentTimeIrr*EnzoTimeStep;
            // ASSERT(dt>=0.0); // Eunwoo debug // Is this necessary?
            xcm[0] = _particle_cm.Position[0] + dt*(_particle_cm.Velocity[0] + 0.5*dt*(_particle_cm.a_irr[0][0] + inv3*dt*_particle_cm.a_irr[0][1]));
            xcm[1] = _particle_cm.Position[1] + dt*(_particle_cm.Velocity[1] + 0.5*dt*(_particle_cm.a_irr[1][0] + inv3*dt*_particle_cm.a_irr[1][1]));
            xcm[2] = _particle_cm.Position[2] + dt*(_particle_cm.Velocity[2] + 0.5*dt*(_particle_cm.a_irr[2][0] + inv3*dt*_particle_cm.a_irr[2][1]));

            Float mcm = _particle_cm.Mass;
            // auto& chi = _particle_cm.changeover;

#ifdef AR_SLOWDOWN_TIMESCALE
            // velocity dependent method 
            Float vp[3], vcm[3];

            vcm[0] = _particle_cm.Velocity[0] + dt*(_particle_cm.a_irr[0][0] + 0.5*dt*_particle_cm.a_irr[0][1]);
            vcm[1] = _particle_cm.Velocity[1] + dt*(_particle_cm.a_irr[1][0] + 0.5*dt*_particle_cm.a_irr[1][1]);
            vcm[2] = _particle_cm.Velocity[2] + dt*(_particle_cm.a_irr[2][0] + 0.5*dt*_particle_cm.a_irr[2][1]);
#endif

            for (int j=0; j<n_pert; j++) {
                Particle* pertj;
                pertj = pert_adr[j];

                Float dt = _time - pertj->CurrentTimeIrr*EnzoTimeStep;
                // ASSERT(dt>=0.0); // Eunwoo debug // Is this necessary?
                xp[0] = pertj->Position[0] + dt*(pertj->Velocity[0] + 0.5*dt*(pertj->a_irr[0][0] + inv3*dt*pertj->a_irr[0][1]));
                xp[1] = pertj->Position[1] + dt*(pertj->Velocity[1] + 0.5*dt*(pertj->a_irr[1][0] + inv3*dt*pertj->a_irr[1][1]));
                xp[2] = pertj->Position[2] + dt*(pertj->Velocity[2] + 0.5*dt*(pertj->a_irr[2][0] + inv3*dt*pertj->a_irr[2][1]));

                Float mj = pertj->Mass;

                // auto& chj = pertj->changeover;

                Float dr[3] = {xp[0] - xcm[0],
                               xp[1] - xcm[1],
                               xp[2] - xcm[2]};

                Float r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
                Float r = sqrt(r2);
                // Float k  = ChangeOver::calcAcc0WTwo(chi, chj, r);
                // _pert_out += calcPertFromMR(r, mcm, k*mj);
                _pert_out += calcPertFromMR(r, mcm, mj);

#ifdef AR_SLOWDOWN_TIMESCALE
                // velocity dependent method 
                vp[0] = pertj->Velocity[0] + dt*(pertj->a_irr[0][0] + 0.5*dt*pertj->a_irr[0][1]);
                vp[1] = pertj->Velocity[1] + dt*(pertj->a_irr[1][0] + 0.5*dt*pertj->a_irr[1][1]);
                vp[2] = pertj->Velocity[2] + dt*(pertj->a_irr[2][0] + 0.5*dt*pertj->a_irr[2][1]);

                Float dv[3] = {vp[0] - vcm[0],
                               vp[1] - vcm[1],
                               vp[2] - vcm[2]};

                // identify whether hyperbolic or closed orbit
                Float gm = gravitational_constant*(mcm+mj);

                calcSlowDownTimeScale(_t_min_sq, dv, dr, r, gm);
#endif
            }
        }
        else {
            _pert_out = 0.0;
            _t_min_sq = 0.0;
        }

        // add soft perturbation
        _pert_out += _perturber.soft_pert_min;

    }
#endif

    //! (Necessary) modify one particle function
    /*!
      @param[in] _p: particle
      @param[in] _time_now: current time (not physical time, NB unit)
      @param[in] _time_end: required evolved time (not physical time, do not directly use, NB unit)
      \return 0: no modification; 1: modify mass; 2: modify mass and velocity; 3: mass become zero
     */
    // template <class Tparticle>
    int modifyOneParticle(Particle& _p, const Float& _time_now, const Float& _time_end) {
        return 0;
    }

  
// This function is used if manager->interrupt_detection_option > 0
//! (Necessary) modify the orbits and interrupt check 
    /*! check the inner left binary whether their separation is smaller than particle radius sum and become close, if true, set one component stauts to merger with cm mass and the other unused with zero mass. Return the binary tree address 
      @param[in] _bin_interrupt: interrupt binary information: adr: binary tree address; time_now: current physical time; time_end: integration finishing time; status: interrupt status: change, merge,none
      @param[in] _bin: binarytree to check iteratively
     */
    // void modifyAndInterruptIter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, REAL dt) {
    void modifyAndInterruptIter(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin) {

        if (_bin_interrupt.status==AR::InterruptStatus::none) {
            auto* p1 = _bin.getLeftMember();
            auto* p2 = _bin.getRightMember();

            if (p1->Mass < p2->Mass) 
                std::swap(p1, p2); // p1 should have the larger mass than p2 (p1->Mass > p2->Mass)

            // if (p1->PID == 75) {
            //     fprintf(mergerout, "APID: %d. time: %e Myr\n", p1->PID, _bin_interrupt.time_now*1e4);
            //     fprintf(mergerout, "APID: %d. x: %e, y: %e, z: %e\n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            //     fflush(mergerout);
            // }


            if(_bin.getMemberN()==2) {

                Float radius = mergerRadius(p1, p2); // p1->radius + p2->radius;

                if (p1->getBinaryInterruptState()== BinaryInterruptState::collision && 
                    p2->getBinaryInterruptState()== BinaryInterruptState::collision &&
                    (p1->time_check<_bin_interrupt.time_end || p2->time_check<_bin_interrupt.time_end) &&
                    (p1->getBinaryPairID()==p2->PID||p2->getBinaryPairID()==p1->PID)) {

                        merge(_bin_interrupt, _bin);
                    }
                else {
                    Float semi, ecc, dr, drdv;
                    _bin.particleToSemiEcc(semi, ecc, dr, drdv, *_bin.getLeftMember(), *_bin.getRightMember(), gravitational_constant);
                    Float peri = semi*(1 - ecc);

                    // if (peri<radius && p1->getBinaryPairID()!=p2->PID&&p2->getBinaryPairID()!=p1->PID) { // original
                    if (peri<radius) { // Eunwoo modified
                        Float ecc_anomaly  = _bin.calcEccAnomaly(dr);
                        Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, ecc);
                        Float mean_motion  = sqrt(gravitational_constant*_bin.Mass/(fabs(_bin.semi*_bin.semi*_bin.semi))); 
                        Float t_peri = mean_anomaly/mean_motion;
                        if (drdv<0 && t_peri<_bin_interrupt.time_end-_bin_interrupt.time_now) {
                            fprintf(mergerout, "Merger1. peri: %e pc, radius: %e pc\n", peri*position_unit, radius*position_unit);
                            merge(_bin_interrupt, _bin);
                        }
                            
                        else if (semi>0||(semi<0&&drdv<0)) {
                            p1->setBinaryPairID(p2->PID);
                            p2->setBinaryPairID(p1->PID);
                            p1->setBinaryInterruptState(BinaryInterruptState::collision);
                            p2->setBinaryInterruptState(BinaryInterruptState::collision);
                            p1->time_check = std::min(p1->time_check, _bin_interrupt.time_now + (drdv<0 ? t_peri : (_bin.period - t_peri)));
                            p2->time_check = std::min(p1->time_check, p2->time_check);
                            fprintf(mergerout, "Merger2. PID: %d and %d might merge soon!\n", p1->PID, p2->PID);
                            fprintf(mergerout, "peri: %e pc, radius: %e pc\n", peri*position_unit, radius*position_unit);
                        }
                    }
                }
            }
            else {
                for (int k=0; k<2; k++) 
                    if (_bin.isMemberTree(k)) modifyAndInterruptIter(_bin_interrupt, *_bin.getMemberAsTree(k)); // Eunwoo added dt for orbit shrinking
                    // if (_bin.isMemberTree(k)) modifyAndInterruptIter(_bin_interrupt, *_bin.getMemberAsTree(k), dt); // Eunwoo added dt for orbit shrinking
            }
        }
    }

// This function is used if manager->interrupt_detection_option > 0
//! Modify the orbits and interrupt check for Kepler solver 
    /*! check the inner left binary whether their separation is smaller than particle radius sum and become close, if true, set one component stauts to merger with cm mass and the other unused with zero mass. Return the binary tree address 
      @param[in] _bin_interrupt: interrupt binary information: adr: binary tree address; time_now: current physical time; time_end: integration finishing time; status: interrupt status: change, merge,none
      @param[in] _bin: binarytree to check iteratively
      @param[in] _dt: integration time
     */
    void modifyAndInterruptKepler(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin, REAL _dt) {

        if (_bin_interrupt.status==AR::InterruptStatus::none) {
            auto* p1 = _bin.getLeftMember();
            auto* p2 = _bin.getRightMember();

            if (p1->Mass < p2->Mass) 
                std::swap(p1, p2); // p1 should have the larger mass than p2 (p1->Mass > p2->Mass)

            Float radius = mergerRadius(p1, p2);

            Float semi, ecc, dr, drdv;
            _bin.particleToSemiEcc(semi, ecc, dr, drdv, *_bin.getLeftMember(), *_bin.getRightMember(), gravitational_constant);
            Float peri = semi*(1 - ecc);
            Float ecc_anomaly  = _bin.calcEccAnomaly(dr); // 0 ~ pi
            Float mean_anomaly = _bin.calcMeanAnomaly(ecc_anomaly, ecc); // 0 ~ pi
            Float mean_motion  = sqrt(gravitational_constant*_bin.Mass/(fabs(_bin.semi*_bin.semi*_bin.semi))); 
            Float t_peri = abs(mean_anomaly/mean_motion); // always smaller than period / 2
            Float period = 2*M_PI/mean_motion;

            if (peri < radius && _dt > t_peri) {
                fprintf(mergerout, "peri: %e pc, radius: %e pc\n", peri*position_unit, radius*position_unit);
                merge(_bin_interrupt, _bin);
            }
        }
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

    Float mergerRadius(Particle *p1, Particle *p2) {

        Float radius = 0.0;
#ifdef SEVN
        if (p1->ParticleType == (Blackhole+SingleParticle) && p2->ParticleType == (Blackhole+SingleParticle)) {
            radius = (p1->radius > p2->radius) ? 3*p1->radius : 3*p2->radius; // r_ISCO == 3 * Schwartzschild raiuds
        }
        else if (p1->ParticleType == (Blackhole+SingleParticle) && p2->ParticleType == (NormalStar+SingleParticle)) {
            radius = 1.3*pow((p1->Mass + p2->Mass)/p2->Mass, 1./3)*p2->radius; // TDE radius
        }
        else if (p1->ParticleType == (NormalStar+SingleParticle) && p2->ParticleType == (Blackhole+SingleParticle)) {
            radius = 1.3*pow((p1->Mass + p2->Mass)/p1->Mass, 1./3)*p1->radius; // TDE radius
        }
        else if (p1->ParticleType == (NormalStar+SingleParticle) && p2->ParticleType == (NormalStar+SingleParticle)) {
            radius = p1->radius + p2->radius; // Sum of two stellar radius
        }
#else
        if (p1->ParticleType == (Blackhole+SingleParticle) && p2->ParticleType == (Blackhole+SingleParticle)) {
            radius = (p1->radius > p2->radius) ? p1->radius : p2->radius; // r_ISCO == 3 * Schwartzschild raiuds
        }
        else if (p1->ParticleType == (Blackhole+SingleParticle) && p2->ParticleType == (NormalStar+SingleParticle)) {
            radius = 1.3*pow((p1->Mass + p2->Mass)/p2->Mass, 1./3)*p2->radius; // TDE radius
        }
        else if (p1->ParticleType == (NormalStar+SingleParticle) && p2->ParticleType == (NormalStar+SingleParticle)) {
            radius = p1->radius + p2->radius; // Sum of two stellar radius
        }
#endif
        return radius;
    }

    void merge(AR::InterruptBinary<Particle>& _bin_interrupt, AR::BinaryTree<Particle>& _bin) { // Stellar merger

        _bin_interrupt.adr = &_bin;
        auto* p1 = _bin.getLeftMember();
        auto* p2 = _bin.getRightMember();

        if (p1->Mass < p2->Mass) 
            std::swap(p1, p2); // p1 should have the larger mass than p2 (p1->Mass > p2->Mass)

        Float radius = mergerRadius(p1, p2);

        if (p1->ParticleType == (Blackhole+SingleParticle) && p2->ParticleType == (Blackhole+SingleParticle)) {
            // radius = (p1->radius > p2->radius) ? p1->radius : p2->radius; // r_ISCO == 3 * Schwartzschild raiuds
            // fprintf(mergerout, "Separation: %e pc\n", dist(p1->Position, p2->Position)*position_unit);
            // fprintf(mergerout, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
            fprintf(mergerout, "r_ISCO: %e pc\n", radius*position_unit);

            fprintf(mergerout, "GW driven merger happens!!! (PID: %d, PID: %d)\n", p1->PID, p2->PID);
            fprintf(mergerout, "Time: %e Myr\n", _bin_interrupt.time_now*1e4);
            fprintf(mergerout, "In center-of-mass frame...\n");
            fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->PID, p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p1->PID, p1->Mass*mass_unit);
            fprintf(mergerout, "PID: %d. Dimensionless spin - %e, %e, %e\n", p1->PID, p1->a_spin[0], p1->a_spin[1], p1->a_spin[2]);
            fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
            fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->PID, p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p2->PID, p2->Mass*mass_unit);
            fprintf(mergerout, "PID: %d. Dimensionless spin - %e, %e, %e\n", p2->PID, p2->a_spin[0], p2->a_spin[1], p2->a_spin[2]);

            _bin_interrupt.status = AR::InterruptStatus::GWmerge;

            Float mcm = p1->Mass + p2->Mass;
            for (int k=0; k<3; k++) {
                p1->Position[k] = (p1->Mass*p1->Position[k] + p2->Mass*p2->Position[k])/mcm;
                p1->Velocity[k] = (p1->Mass*p1->Velocity[k] + p2->Mass*p2->Velocity[k])/mcm;
                p2->Position[k] = 0.0;
                p2->Velocity[k] = 0.0;
            }
            recoilKick(_bin, p1, p2);
            remnantSpinMass(_bin, p1, p2);

            p1->setBinaryInterruptState(BinaryInterruptState::none);
            p2->setBinaryInterruptState(BinaryInterruptState::none);
            // p1->dm = mcm - p1->Mass;
            // p2->dm = -p2->Mass;
            // p1->Mass = mcm;
            p2->Mass = 0.0;
            fprintf(mergerout, "---------------Merger remnant properties---------------\n");
            fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
            fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
        }
        else if ((p1->ParticleType == (Blackhole+SingleParticle) && p2->ParticleType == (NormalStar+SingleParticle)) ||
                (p1->ParticleType == (NormalStar+SingleParticle) && p2->ParticleType == (Blackhole+SingleParticle))) {
            // radius = 1.3*pow((p1->Mass + p2->Mass)/p2->Mass, 1./3)*p2->radius; // TDE radius

            if (p2->ParticleType == (Blackhole+SingleParticle)) 
                std::swap(p1, p2); // p1 should be BH

            // fprintf(mergerout, "Separation: %e pc\n", dist(p1->Position, p2->Position)*position_unit);
            // fprintf(mergerout, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
            fprintf(mergerout, "r_TDE: %e pc\n", radius*position_unit);

            fprintf(mergerout, "TDE happens!!! (PID: %d, PID: %d)\n", p1->PID, p2->PID);
            fprintf(mergerout, "Time: %e Myr\n", _bin_interrupt.time_now*1e4);
            fprintf(mergerout, "In center-of-mass frame...\n");
            fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->PID, p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p1->PID, p1->Mass*mass_unit);
            fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
            fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->PID, p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p2->PID, p2->Mass*mass_unit);

            _bin_interrupt.status = AR::InterruptStatus::TDE;
            Float mcm = p1->Mass + p2->Mass * 0.5; // If TDE happens, the half of the mass of star is accreted to a BH.
            for (int k=0; k<3; k++) {
                p1->Position[k] = (p1->Mass*p1->Position[k] + p2->Mass*p2->Position[k])/mcm;
                p1->Velocity[k] = (p1->Mass*p1->Velocity[k] + p2->Mass*p2->Velocity[k])/mcm;
            }
            p1->setBinaryInterruptState(BinaryInterruptState::none);
            p2->setBinaryInterruptState(BinaryInterruptState::none);
            // p1->dm = mcm - p1->Mass;
            p1->dm += 0.5 * p2->Mass;
            p1->Mass = mcm;
            p2->Mass = 0.0;
            fprintf(mergerout, "---------------Merger remnant properties---------------\n");
            fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
            fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
        }
        else if (p1->ParticleType == (NormalStar+SingleParticle) && p2->ParticleType == (NormalStar+SingleParticle)) {
            // radius = p1->radius + p2->radius; // Sum of two stellar radius
            // fprintf(mergerout, "Separation: %e pc\n", dist(p1->Position, p2->Position)*position_unit);
            // fprintf(mergerout, "peri: %e pc\n", _bin.semi*(1 - _bin.ecc)*position_unit);
            fprintf(mergerout, "r1 + r2: %e pc\n", radius*position_unit);

            fprintf(mergerout, "Stellar merger happens!!! (PID: %d, PID: %d)\n", p1->PID, p2->PID);
            fprintf(mergerout, "Time: %e Myr\n", _bin_interrupt.time_now*1e4);
            fprintf(mergerout, "In center-of-mass frame...\n");
            fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p1->PID, p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
            fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->PID, p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p1->PID, p1->Mass*mass_unit);
            fprintf(mergerout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", p2->PID, p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
            fprintf(mergerout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->PID, p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(mergerout, "PID: %d. Mass (Msol) - %e, \n", p2->PID, p2->Mass*mass_unit);

            _bin_interrupt.status = AR::InterruptStatus::merge;
            Float mcm = p1->Mass + p2->Mass;
            for (int k=0; k<3; k++) {
                p1->Position[k] = (p1->Mass*p1->Position[k] + p2->Mass*p2->Position[k])/mcm;
                p1->Velocity[k] = (p1->Mass*p1->Velocity[k] + p2->Mass*p2->Velocity[k])/mcm;
                p2->Position[k] = p1->Position[k];
                p2->Velocity[k] = p1->Velocity[k];
                
            }
            p1->setBinaryInterruptState(BinaryInterruptState::none);
            p2->setBinaryInterruptState(BinaryInterruptState::none);
#ifdef SEVN
            if (p1->star == nullptr && p2->star == nullptr) {

                // p1->dm = mcm - p1->Mass;
                // p2->dm = -p2->Mass;
                p1->Mass = mcm;
                p2->Mass = 0.0;

                if (mcm*mass_unit > 2.2 && mcm*mass_unit < 600) {

                    std::vector<std::string> args = {"empty", // Not used
                                                    // "-myself", "/data/vinicius/NbodyPlus/SEVN",
                                                    "-tables", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_parsec_ov04_AGB", 
                                                    //  "-tables", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_MIST_AGB",
                                                    // "-tables_HE", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_parsec_pureHe36",
                                                    // "-turn_WR_to_pureHe", "false",
                                                    "-snmode", "delayed",
                                                    "-Z", "0.0002",
                                                    "-spin", "0.0",
                                                    "-tini", "zams", 
                                                    "-tf", "end",
                                                    // "-tf", "0.000122",
                                                    "-dtout", "events",
                                                    "-xspinmode", "geneva"};
                    std::vector<char*> c_args;
                    for (auto& arg : args) {
                        c_args.push_back(&arg[0]);
                    }

                    IO* sevnio; // Eunwoo: global variable -> We can initialize Star and Binstar class anywhere.
                    sevnio = new IO;
                    sevnio->load(c_args.size(), c_args.data());

                    std::vector<std::string> init_params{std::to_string(double(p1->Mass*mass_unit)), "0.0002", "0.0", "delayed", "zams", "end", "events"};
                    size_t id = p1->PID;
                    p1->star = new Star(sevnio, init_params, id, false);
                    p1->EvolutionTime = _bin_interrupt.time_now*1e4;
                    p1->radius = p1->star->getp(Radius::ID)/(utilities::parsec_to_Rsun)/position_unit;
                    fprintf(SEVNout, "New Star class made!\n");
                    fprintf(SEVNout, "PID: %d. Mass: %e Msol, Radius: %e pc\n", p1->PID, p1->Mass*mass_unit, p1->radius*position_unit);
                }
                fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
                fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
                fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
                fprintf(mergerout, "Type: %d, \n", int(p1->star->getp(Phase::ID)));
                fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
            }
            else if (p1->star != nullptr && p2->star != nullptr) {

                Mix(p1->star, p2->star);
                UpdateEvolution(p1);
                UpdateEvolution(p2);
                fprintf(SEVNout, "Mix done!\n");

                if (p1->star->amiempty() && !p2->star->amiempty()) {
                    // p1->dm -= p1->Mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
                    p1->Mass = 0.0;
                    fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                    fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
                    fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
                    fprintf(mergerout, "Mass (Msol) - %e, \n", p2->Mass*mass_unit);
                    if (p2->star->amiremnant())
                        fprintf(mergerout, "Type: %d, \n", 8 + int(p2->star->getp(RemnantType::ID)));
                    else
                        fprintf(mergerout, "Type: %d, \n", int(p2->star->getp(Phase::ID)));
                    fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
                }
                else if (!p1->star->amiempty() && p2->star->amiempty()) {
                    // p2->dm -= p2->Mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
                    p2->Mass = 0.0;
                    fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                    fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
                    fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
                    fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
                    if (p1->star->amiremnant())
                        fprintf(mergerout, "Type: %d, \n", 8 + int(p1->star->getp(RemnantType::ID)));
                    else
                        fprintf(mergerout, "Type: %d, \n", int(p1->star->getp(Phase::ID)));
                    fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
                }
                else if (p1->star->amiempty() && p2->star->amiempty()) { // Type Ia supernova
                    p1->dm += p1->Mass; // Eunwoo: dm should be 0 after it distributes its mass to the nearby gas cells.
                    p1->Mass = 0.0;
                    p2->dm += p2->Mass;
                    p2->Mass = 0.0;
                    fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                    fprintf(mergerout, "Type Ia Supernova event! Both of the stars becomes empty!\n");
                    fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
                }
                else
                    throw std::runtime_error("None of stars are empty: Something wrong in stellar merger!");
            }
            else if (p1->star == nullptr && p2->star != nullptr) {

                p2->star->update_from_binary(Mass::ID, p1->Mass*mass_unit);
                p2->star->update_from_binary(dMcumul_binary::ID, p1->Mass*mass_unit);
                if (p2->star->aminakedhelium())
                    p2->star->jump_to_normal_tracks();
                else
                    p2->star->find_new_track_after_merger();
                UpdateEvolution(p2);
                p1->Mass = 0.0;

                fprintf(SEVNout, "Mix with no done!\n");
                fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p2->Position[0]*position_unit, p2->Position[1]*position_unit, p2->Position[2]*position_unit);
                fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p2->Velocity[0]*velocity_unit/yr*pc/1e5, p2->Velocity[1]*velocity_unit/yr*pc/1e5, p2->Velocity[2]*velocity_unit/yr*pc/1e5);
                fprintf(mergerout, "Mass (Msol) - %e, \n", p2->Mass*mass_unit);
                if (p2->star->amiremnant())
                    fprintf(mergerout, "Type: %d, \n", 8 + int(p2->star->getp(RemnantType::ID)));
                else
                    fprintf(mergerout, "Type: %d, \n", int(p2->star->getp(Phase::ID)));
                fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
            }
            else if (p1->star != nullptr && p2->star == nullptr) {

                p1->star->update_from_binary(Mass::ID, p2->Mass*mass_unit);
                p1->star->update_from_binary(dMcumul_binary::ID, p2->Mass*mass_unit);
                if (p1->star->aminakedhelium())
                    p1->star->jump_to_normal_tracks();
                else
                    p1->star->find_new_track_after_merger();
                UpdateEvolution(p1);
                p2->Mass = 0.0;

                fprintf(SEVNout, "Mix with no done!\n");
                fprintf(mergerout, "---------------Merger remnant properties---------------\n");
                fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
                fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
                fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
                if (p1->star->amiremnant())
                    fprintf(mergerout, "Type: %d, \n", 8 + int(p1->star->getp(RemnantType::ID)));
                else
                    fprintf(mergerout, "Type: %d, \n", int(p1->star->getp(Phase::ID)));
                fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
            }
        }
        fflush(mergerout);
        fflush(SEVNout);
#else
            p1->dm = mcm - p1->Mass;
            p2->dm = -p2->Mass;
            p1->Mass = mcm;
            p2->Mass = 0.0;
        }
        // fprintf(mergerout, "---------------Merger remnant properties---------------\n");
        // fprintf(mergerout, "Position (pc) - x:%e, y:%e, z:%e, \n", p1->Position[0]*position_unit, p1->Position[1]*position_unit, p1->Position[2]*position_unit);
        // fprintf(mergerout, "Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
        // fprintf(mergerout, "Mass (Msol) - %e, \n", p1->Mass*mass_unit);
        // fprintf(mergerout, "---------------------END-OF-MERGER---------------------\n\n");
        fflush(mergerout);   
#endif
    }

    // Reference for remnant spin: Hofmann et al. (2016) (https://iopscience.iop.org/article/10.3847/2041-8205/825/2/L19/pdf)
    // Using eq (2)-(6), (13)-(16)
    // n_M = 3, n_J = 4
    // Reference for remnant mass: Barausse et al. (2012) (https://iopscience.iop.org/article/10.1088/0004-637X/758/1/63/pdf)
    // Using eq (1)-(5), (12), (15)-(18)
    void remnantSpinMass(AR::BinaryTree<Particle>& _bin, Particle* p1, Particle* p2) {

        REAL k[4][5] = {{-5.9, 3.39221, 4.48865, -5.77101, -13.0459},
                        {35.1287, -72.9336, -86.0036, 93.7371, 200.975},
                        {-146.822, 387.184, 447.009, -467.383, -884.339},
                        {223.911, -648.502, -697.177, 753.738, 1166.89}};
        REAL ksi = 0.474046;

        REAL a1 = sqrt(mag(p1->a_spin));
        REAL a2 = sqrt(mag(p2->a_spin));
        REAL Mtot = p1->Mass + p2->Mass;
        REAL q = p2->Mass/p1->Mass;
        assert(q <= 1);
        REAL nu = q/(1+q)/(1+q);
        const REAL c = 299752.458 / (velocity_unit / yr * pc / 1e5);

        REAL L_ang[3]; // specific angular momentum (r_rel x v_rel)
        L_ang[0] = _bin.am.x;
        L_ang[1] = _bin.am.y;
        L_ang[2] = _bin.am.z;

        fprintf(mergerout, "L_orbit: (%e, %e, %e)\n", (p1->Mass*p2->Mass/Mtot) * L_ang[0],
                                                        (p1->Mass*p2->Mass/Mtot) * L_ang[1],
                                                        (p1->Mass*p2->Mass/Mtot) * L_ang[2]);
        fprintf(mergerout, "S_1: (%e, %e, %e)\n", (p1->Mass*p1->Mass/c) * p1->a_spin[0], 
                                                    (p1->Mass*p1->Mass/c) * p1->a_spin[1], 
                                                    (p1->Mass*p1->Mass/c) * p1->a_spin[2]);
        fprintf(mergerout, "S_2: (%e, %e, %e)\n", (p2->Mass*p2->Mass/c) * p2->a_spin[0], 
                                                    (p2->Mass*p2->Mass/c) * p2->a_spin[1], 
                                                    (p2->Mass*p2->Mass/c) * p2->a_spin[2]);
        fprintf(mergerout, "In code unit!\n");

        auto cosine = [&](REAL a[3], REAL b[3]) {
            if (mag(a) == 0 || mag(b) == 0)
                return 0.;
            else
                return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])/sqrt(mag(a))/sqrt(mag(b));
        };

        REAL cosa = cosine(p1->a_spin, p2->a_spin); // cos(alpha): This angle should be changed to the initial value. I will change this later.
        REAL cosb = cosine(L_ang, p1->a_spin);      // cos(beta)
        REAL cosg = cosine(L_ang, p2->a_spin);      // cos(gamma)

        REAL atot = (a1*cosb + a2*cosg*q*q)/(1 + q)/(1 + q); // atilde in Barausse et al. (2012)
        REAL aeff = atot + ksi*nu*(a1*cosb + a2*cosg);

        auto Z1 = [&](REAL a) {
            return 1 + pow(1 - a*a, 1./3)*(pow(1 + a, 1./3) + pow(1 - a, 1./3));
        };
        auto Z2 = [&](REAL a, REAL Z1) {
            return sqrt(3*a*a + Z1*Z1);
        };
        auto r_ISCO = [&](REAL a, REAL Z1, REAL Z2) {
            if (a > 0)
                return 3 + Z2 - sqrt(3 - Z1)*sqrt(3 + Z1 + 2*Z2);
            else if (a < 0)
                return 3 + Z2 + sqrt(3 - Z1)*sqrt(3 + Z1 + 2*Z2);
            else
                return 3 + Z2;
        };
        auto E_ISCO = [&](REAL r) {
            return sqrt(1 - 2./3/r);
        };

        REAL Z1_aeff = Z1(aeff);
        REAL Z2_aeff = Z2(aeff, Z1_aeff);
        REAL r_ISCO_aeff = r_ISCO(aeff, Z1_aeff, Z2_aeff);
        REAL E_ISCO_aeff = E_ISCO(r_ISCO_aeff);
        REAL L_ISCO_aeff = 2/3/sqrt(3)*(1 + 2*sqrt(3*r_ISCO_aeff - 2));

        REAL l = L_ISCO_aeff - 2 * atot * (E_ISCO_aeff - 1);
        for (int i = 0; i <= 3; ++i) {      // n_M = 3
            for (int j = 0; j <= 4; ++j) {  // n_J = 4
                l += k[i][j] * pow(nu, 1 + i) * pow(aeff, j);
            }
        }
        l = abs(l);

        REAL afin = 1./pow(1 + q, 2)*sqrt(a1*a1 + a2*a2*pow(q, 4) + 2*a1*a2*q*q*cosa + 2*(a1*cosb + a2*q*q*cosg)*l*q + l*l*q*q);

        REAL J_tot[3];
        for (int i=0; i<3; i++)
            J_tot[i] = (p1->Mass*p2->Mass/Mtot) * L_ang[i] + (p1->Mass*p1->Mass/c) * p1->a_spin[i] + (p2->Mass*p2->Mass/c) * p2->a_spin[i];

        REAL J_tot_norm[3];
        for (int i=0; i<3; i++)
            J_tot_norm[i] = J_tot[i]/sqrt(mag(J_tot));

        for (int i=0; i<3; i++)
            p1->a_spin[i] = afin*J_tot_norm[i];

        fprintf(mergerout, "Dimensionless spin of remnant BH: (%e, %e, %e)\n", p1->a_spin[0], p1->a_spin[1], p1->a_spin[2]);

        REAL Z1_atot = Z1(atot);
        REAL Z2_atot = Z2(atot, Z1_atot);
        REAL r_ISCO_atot = r_ISCO(atot, Z1_atot, Z2_atot);
        REAL E_ISCO_atot = E_ISCO(r_ISCO_atot);

        REAL Erad = (1 - E_ISCO_atot)*nu + 4*nu*nu*(4*0.04827 + 16*0.01707*atot*(atot+1) + E_ISCO_atot - 1);

        p1->Mass = (1 - Erad) * Mtot;

        fprintf(mergerout, "Mass of remnant BH: %e Msol\n", p1->Mass*mass_unit);
        fflush(mergerout);
    }

    // Reference: Arca Sedda et al. (2020) (https://iopscience.iop.org/article/10.3847/1538-4357/ab88b2/pdf)
    void recoilKick(AR::BinaryTree<Particle>& _bin, Particle* p1, Particle* p2) {
        
        REAL A      = 1.2e4; // km/s
        REAL B      = -0.93;
        REAL H      = 6.9e3; // km/s
        REAL ksi    = 145 * M_PI / 180; // rad
        REAL V11    = 3677.76; // km/s
        REAL VA     = 2.481e3; // km/s
        REAL VB     = 1.793e3; // km/s
        REAL VC     = 1.507e3; // km/s

        auto cosine = [&](REAL a[3], REAL b[3]) {
            if (mag(a) == 0 || mag(b) == 0)
                return 0.;
            else
                return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])/sqrt(mag(a))/sqrt(mag(b));
        };

        REAL a1 = sqrt(mag(p1->a_spin));
        REAL a2 = sqrt(mag(p2->a_spin));
        REAL q = p2->Mass/p1->Mass;
        // fprintf(mergerout, "q: %e\n", q);
        assert(q <= 1);
        REAL nu = q/(1+q)/(1+q);
        // fprintf(mergerout, "nu: %e\n", nu);

        REAL L_ang[3];
        L_ang[0] = _bin.am.x;
        L_ang[1] = _bin.am.y;
        L_ang[2] = _bin.am.z;

        REAL cosa = cosine(p1->a_spin, p2->a_spin); // cos(alpha): This angle should be changed to the initial value. I will change this later.
        REAL sina = sin(acos(cosa));
        REAL cosb = cosine(L_ang, p1->a_spin);      // cos(beta)
        REAL sinb = sin(acos(cosb));
        REAL cosg = cosine(L_ang, p2->a_spin);      // cos(gamma)
        REAL sing = sin(acos(cosg));
        // fprintf(mergerout, "cosa: %e, sina: %e, cosb: %e, sinb: %e, cosg: %e, sing: %e\n", cosa, sina, cosb, sinb, cosg, sing);

        REAL a2par  = a2 * cosg;
        // fprintf(mergerout, "a2par: %e\n", a2par);
        REAL a2per1 = a2 * sing;
        // fprintf(mergerout, "a2per1, %e\n", a2per1);
        REAL a2per2 = 0.0;
        // fprintf(mergerout, "a2per2, %e\n", a2per2);
      
        REAL a1par  = a1 * cosb;
        // fprintf(mergerout, "a1par: %e\n", a1par);
        REAL a1per1 = a1 * sinb*cosa;
        // fprintf(mergerout, "a1per1: %e\n", a1per1);
        REAL a1per2 = a1 * sinb*sina;
        // fprintf(mergerout, "a1per2: %e\n", a1per2);

        REAL KSIpar = 2 * (a2par + q*q*a1par) / (1 + q) / (1 + q);
        // fprintf(mergerout, "KSIpar: %e\n", KSIpar);
        std::random_device rd; // Obtain a random number from hardware
        std::mt19937 mt(rd()); // Seed the generator
        std::uniform_real_distribution<> distr(0.0, 1.0); // Define the range (0 to 1)
        REAL phi = 2 * M_PI * distr(mt); // phi_Delta - phi_1
        // fprintf(mergerout, "phi: %e\n", phi);

        REAL vm = A*nu*nu*sqrt(1 - 4*nu) * (1 + B*nu);
        // fprintf(mergerout, "vm: %e\n", vm);
        REAL vper = H*nu*nu / (1 + q) * (a2par - q * a1par);
        // fprintf(mergerout, "vper: %e\n", vper);
        REAL vpar = 16*nu*nu / (1 + q) * (V11 + VA*KSIpar + VB*KSIpar*KSIpar + VC*KSIpar*KSIpar*KSIpar);
        vpar *= sqrt((a2per1 - q * a1per1) * (a2per1 - q * a1per1) + (a2per2 - q * a1per2) * (a2per2 - q * a1per2)) * cos(phi);
        // fprintf(mergerout, "vpar: %e\n", vpar);

        REAL e_par[3]   = {L_ang[0]/sqrt(mag(L_ang)), L_ang[1]/sqrt(mag(L_ang)), L_ang[2]/sqrt(mag(L_ang))};
        REAL e_per1[3]  = {p2->a_spin[0] - p2->a_spin[0]*cosg, p2->a_spin[1] - p2->a_spin[1]*cosg, p2->a_spin[2] - p2->a_spin[2]*cosg};
        REAL norm1      = sqrt(mag(e_per1));
        if (norm1 != 0) {
            for (int i=0; i<3; i++)
                e_per1[i]   /= norm1; 
        }
               
        REAL e_per2[3]  =   {e_par[1] * e_per1[2] - e_par[2] * e_per1[1], 
                            e_par[2] * e_per1[0] - e_par[0] * e_per1[2], 
                            e_par[0] * e_per1[1] - e_par[1] * e_per1[0]}; // cross product: e2 = e3 x e1
        // fprintf(mergerout, "e1: %e, %e, %e\n", e_per1[0], e_per1[1], e_per1[2]);
        // fprintf(mergerout, "e2: %e, %e, %e\n", e_per2[0], e_per2[1], e_per2[2]);
        // fprintf(mergerout, "e3: %e, %e, %e\n", e_par[0], e_par[1], e_par[2]);

        REAL vkick[3];
        for (int i=0; i<3; i++)
            vkick[i] = (vm + vper*cos(ksi)) * e_per1[i] + vper*sin(ksi) * e_per2[i] + vpar * e_par[i];

        fprintf(mergerout, "GW recoil kick: (%e, %e, %e) km/s\n", vkick[0], vkick[1], vkick[2]);
        fprintf(mergerout, "\t magnitude: %e km/s\n", sqrt(mag(vkick)));

        for (int i=0; i<3; i++)
            p1->Velocity[i] += vkick[i]/(velocity_unit/yr*pc/1e5); // km/s to code unit

        // fprintf(mergerout, "Remnant velocity: (%e, %e, %e) km/s\n", p1->Velocity[0]*velocity_unit/yr*pc/1e5, p1->Velocity[1]*velocity_unit/yr*pc/1e5, p1->Velocity[2]*velocity_unit/yr*pc/1e5);
        fflush(mergerout);
    }
    
};

