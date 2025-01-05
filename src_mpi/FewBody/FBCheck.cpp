#ifdef FEWBODY
#include "../global.h"
#include "ar_interaction.hpp"

// calculate dr dv of a pair
Float calcDrDv(const double *pos1, const double *pos2, const double *vel1, const double *vel2) {
    Float dx[3],dv[3];
    dx[0] = pos1[0] - pos2[0];
    dx[1] = pos1[1] - pos2[1];
    dx[2] = pos1[2] - pos2[2];

    dv[0] = vel1[0] - vel2[0];
    dv[1] = vel1[1] - vel2[1];
    dv[2] = vel1[2] - vel2[2];

    Float drdv= dx[0]*dv[0] + dx[1]*dv[1] + dx[2]*dv[2];

    return drdv;
}

//! check the pair with distance below r_crit for ptcl in adr_dt_sorted_
/*! First check nearest neighbor distance r_min
    If r_min<r_crit, check the direction, if income, accept as group
*/
void Particle::checkNewGroup() {

    // kappa_org criterion for new group kappa_org>kappa_org_crit
    const Float kappa_org_crit = 1e-2;

    double pos1[Dim], vel1[Dim];
    // this->NewNumberOfNeighbor = 0;

    Particle* ptcl2;

    // REAL r_min = 1; // Eunwoo test;

    // check only active particles 
    // single case
    for (int i=0; i < this->NumberOfNeighbor; i++) {
		ptcl2 = &particles[this->Neighbors[i]];

        // if (ptcl2->TimeStepIrr > this->TimeStepIrr) // test_1e5_4 & 5: this must make the same result!
        if (ptcl2->TimeStepIrr*EnzoTimeStep*1e4 > tbin) // fiducial: 1e-5 but for rbin = 0.00025 pc, 1e-6 Myr seems good
            continue;

        double dt;

        dt = this->CurrentTimeIrr > ptcl2->CurrentTimeIrr ? \
									 this->CurrentTimeIrr - ptcl2->CurrentTimeIrr : ptcl2->CurrentTimeIrr - this->CurrentTimeIrr;

        double pos2[Dim], vel2[Dim];
        
		this->predictParticleSecondOrder(dt, pos1, vel1);
		ptcl2->predictParticleSecondOrder(dt, pos2, vel2);

        const Float dr = dist(pos1, pos2);
        // ASSERT(dr>0.0);

        // if (dr < r_min) // Eunwoo test
        //     r_min = dr; // Eunwoo test

        // distance criterion
        Float r_crit = rbin/position_unit;
        // Float r_crit_sq = r_crit*r_crit;
        
        if (dr < r_crit) {

            // this increase AR total step too much
            //bool add_flag=false;
            //Float semi, ecc, dr, drdv;
            //AR::BinaryTree<Tparticle>::particleToSemiEcc(semi,ecc,dr,drdv, pi, *pj, ar_manager->interaction.gravitational_constant); 
            //else {
            //    Float rp = drdv/dr*dt_limit_ + dr;
            //    if (rp < r_crit) add_flag = true;
            //}

            Float drdv = calcDrDv(pos1, pos2, vel1, vel2);
            // only inwards
            if(drdv<0.0) {

// /* // test_1e4_2
                Float fcm[3] = {this->Mass*this->a_irr[0][0] + ptcl2->Mass*ptcl2->a_irr[0][0], 
                                this->Mass*this->a_irr[1][0] + ptcl2->Mass*ptcl2->a_irr[1][0], 
                                this->Mass*this->a_irr[2][0] + ptcl2->Mass*ptcl2->a_irr[2][0]};

                AR::SlowDown sd;
                Interaction interaction;
                Float mcm = this->Mass + ptcl2->Mass;

                // sd.initialSlowDownReference(ar_manager->slowdown_pert_ratio_ref, ar_manager->slowdown_timescale_max);
                sd.initialSlowDownReference(1e-6, NUMERIC_FLOAT_MAX);

                sd.pert_in = interaction.calcPertFromMR(dr, this->Mass, ptcl2->Mass);
                sd.pert_out = interaction.calcPertFromForce(fcm, mcm, mcm);

                sd.calcSlowDownFactor();
                Float kappa_org = sd.getSlowDownFactorOrigin();

                // avoid strong perturbed case, estimate perturbation
                // if kappa_org < criterion, avoid to form new group, should be consistent as checkbreak
                if(kappa_org<kappa_org_crit) continue;
// */ // test_1e4_2

#ifdef ADJUST_GROUP_DEBUG
                if (j<index_offset_group_) {
                    std::cerr<<"Find new group: time: "<<time_
                                <<" index: "<<i<<" "<<j
                                <<" dr: "<<sqrt(dr2)
                                <<" kappa_org: "<<kappa_org<<"\n";
                }
                else {
                    auto& bin_root = groups[j-index_offset_group_].info.getBinaryTreeRoot();
                    std::cerr<<"Find new group: time: "<<time_
                                <<" dr: "<<sqrt(dr2)
                                <<" kappa_org: "<<kappa_org<<"\n"
                                <<"       index         slowdown          apo \n"
                                <<"i1 "
                                <<std::setw(8)<<i
                                <<std::setw(16)<<0
                                <<std::setw(16)<<0;
                    std::cerr<<"\ni2 "
                                <<std::setw(8)<<j
                                <<std::setw(16)<<bin_root.slowdown.getSlowDownFactorOrigin()
                                <<std::setw(16)<<bin_root.semi*(1.0+bin_root.ecc);
                    std::cerr<<std::endl;
                }
                fprintf(binout, "Find new group!\n\t");
                fprintf(binout, " separation: %e pc\n\t", dr);
                fprintf(binout, " kappa_org: %e \n\t", kappa_org);
#endif
                this->NewNeighbors[this->NewNumberOfNeighbor] = this->Neighbors[i];
                this->NewNumberOfNeighbor++;
            }
        }
    }
}

// Use when many-body (>3) group broke
// Or primordial binary search
void Particle::checkNewGroup2() {

    // kappa_org criterion for new group kappa_org>kappa_org_crit
    const Float kappa_org_crit = 1e-2;

    double pos1[Dim], vel1[Dim];
    // this->NewNumberOfNeighbor = 0;

    Particle* ptcl2;

    // check only active particles 
    // single case
    for (int i=0; i < this->NumberOfNeighbor; i++) {
		ptcl2 = &particles[this->Neighbors[i]];

        double dt;

        dt = this->CurrentTimeIrr > ptcl2->CurrentTimeIrr ? \
									 this->CurrentTimeIrr - ptcl2->CurrentTimeIrr : ptcl2->CurrentTimeIrr - this->CurrentTimeIrr;

        double pos2[Dim], vel2[Dim];
        
		this->predictParticleSecondOrder(dt, pos1, vel1);
		ptcl2->predictParticleSecondOrder(dt, pos2, vel2);

        const Float dr = dist(pos1, pos2);
        // ASSERT(dr>0.0);

        double v2 = std::pow(dist(vel1, vel2), 2);
        // determine they are bound or not
        double energy = v2/2 - (this->Mass + ptcl2->Mass)/dr;

        // distance criterion
        Float r_crit = rbin/position_unit;
        
        if (dr < r_crit && energy < 0) {
            this->NewNeighbors[this->NewNumberOfNeighbor] = this->Neighbors[i];
            this->NewNumberOfNeighbor++;
        }
    }
}


// reference: checkBreak in hermite_integrator.h (SDAR)
bool Group::CheckBreak() {

    // fprintf(binout, "CheckBreak opened!\n");

    const Float kappa_org_crit = 1e-2;
    const int n_member = sym_int.particles.getSize();

    // generate binary tree
    sym_int.info.generateBinaryTree(sym_int.particles, manager.interaction.gravitational_constant);

    auto& bin_root = sym_int.info.getBinaryTreeRoot();
    bool outgoing_flag = false; // Indicate whether it is a outgoing case or income case

    // check whether periapsis distance >  2e-3 pc
    // Periapsis is too far and binary is not close enough to use regularized technique.
    // if (bin_root.semi*(1-bin_root.ecc) > 2e-3/position_unit){ // test7
    if (bin_root.semi*(1-bin_root.ecc) > rbin/position_unit && bin_root.r > 2 * rbin/position_unit){ // test8 // fiducial
    // if (bin_root.semi*(1-bin_root.ecc) > 1.2e-3/position_unit){ // test12
        fprintf(binout, "Break group: too far periapsis!\n\t");
        fprintf(binout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
        fprintf(binout, "N_member: %d\n\t", n_member);
        fprintf(binout, "separation: %e pc\n\t", bin_root.r*position_unit);
        fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
        fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
        fprintf(binout, "ecca: %e\n\t", bin_root.ecca);
        fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
        fprintf(binout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
        fprintf(binout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
        fflush(binout);
        return true;
    }

    // check binary case 
    // ecc anomaly indicates outgoing (ecca>0) or income (ecca<0)
    if (bin_root.semi>0.0 && bin_root.ecca>0.0) {
        outgoing_flag = true;
        // check whether separation is larger than distance criterion. 
        // /*
        if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 2 * rbin/position_unit) { // test8 // fiducial
        // if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 1.2e-3/position_unit) { // test12
            fprintf(binout, "Break group: binary escape!\n\t");
            fprintf(binout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
            fprintf(binout, "N_member: %d\n\t", n_member);
            fprintf(binout, "separation: %e pc\n\t", bin_root.r*position_unit);
            fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
            fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
            fprintf(binout, "ecca: %e\n\t", bin_root.ecca);
            fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
            fprintf(binout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
            fprintf(binout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
            fflush(binout);
            return true;
        }
        // */
        /*
        if (n_member == 2) {
            if (bin_root.r > sym_int.info.r_break_crit && bin_root.r > 2e-3/position_unit) {
                fprintf(binout, "Break group: binary escape!\n\t");
                fprintf(binout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                fprintf(binout, "N_member: %d\n\t", n_member);
                fprintf(binout, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
                fprintf(binout, "ecca: %e\n\t", bin_root.ecca);
                fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(binout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(binout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                fflush(binout);
                return true;
            }
        }
        else {
            if (bin_root.r > 2e-3/position_unit) {
                fprintf(binout, "Break group: binary escape!\n\t");
                fprintf(binout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                fprintf(binout, "N_member: %d\n\t", n_member);
                fprintf(binout, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
                fprintf(binout, "ecca: %e\n\t", bin_root.ecca);
                fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(binout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(binout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                fflush(binout);
                return true;
            }
        }
        */
/*
        // in case apo is larger than distance criterion
        Float apo = bin_root.semi * (1.0 + bin_root.ecc);
        if (apo>sym_int.info.r_break_crit) {
            Float dr2, drdv;

            sym_int.info.getDrDv(dr2, drdv, *bin_root.getLeftMember(), *bin_root.getRightMember());
            ASSERT(drdv>=0.0);

            // check whether next step the separation is larger than distance criterion
            // Not sure whether it can work correctly or not:
            // rp = v_r * dt + r
            Float dr = drdv/bin_root.r*sym_int.particles.cm.TimeStepIrr*EnzoTimeStep;
            Float rp =  dr + bin_root.r;
            if (rp >sym_int.info.r_break_crit) {
                // in case r is too small, avoid too early quit of group
                Float rph = 0.5*dr + bin_root.r;
                if ( rph < sym_int.info.r_break_crit && bin_root.r <0.2*sym_int.info.r_break_crit) {
                    // std::cerr<<"Binary will escape but dr is too small, reduce cm step by half, time: "<<CurrentTime<<" N_member: "<<n_member<<" ecca: "<<bin_root.ecca<<" separation : "<<bin_root.r<<" apo: "<<apo<<" r_pred: "<<rp<<" drdv: "<<drdv<<" dt: "<<sym_int.particles.cm.TimeStepIrr<<" r_crit: "<<sym_int.info.r_break_crit<<std::endl;
                    sym_int.particles.cm.TimeLevelIrr--;
                    sym_int.particles.cm.TimeStepIrr = static_cast<REAL>(pow(2, sym_int.particles.cm.TimeLevelIrr));
                    sym_int.particles.cm.TimeBlockIrr = static_cast<ULL>(pow(2, sym_int.particles.cm.TimeLevelIrr-time_block));

                }
                else {
                    // std::cerr<<"Break group: binary will escape, time: "<<CurrentTime<<" N_member: "<<n_member<<" ecca: "<<bin_root.ecca<<" separation : "<<bin_root.r<<" apo: "<<apo<<" r_pred: "<<rp<<" drdv: "<<drdv<<" dt: "<<sym_int.particles.cm.TimeStepIrr<<" r_crit: "<<sym_int.info.r_break_crit<<std::endl;
                    fprintf(binout, "Break group: binary will escape!\n\t time: %e Myr\n\t N_member: %d\n\t ecca: %e\n\t separation: %e pc\n\t r_crit: %e pc\n", CurrentTime*EnzoTimeStep*1e4, n_member, bin_root.ecca, bin_root.r*position_unit, sym_int.info.r_break_crit*position_unit);
                    fflush(binout);
                    return true;
                }
                return false;
            }
        }
*/
    }

    // check hyperbolic case
    if (bin_root.semi<0.0) {
        // hyperbolic case, ecca is not correctly calculated
        Float dr2, drdv;
        sym_int.info.getDrDv(dr2, drdv, *bin_root.getLeftMember(), *bin_root.getRightMember());
        if (drdv>0.0) {
            outgoing_flag = true;
            // check distance criterion
            if (bin_root.r > 2 * rbin/position_unit) { // test8 // seems good! // fiducial
            // if (bin_root.r > 1.2e-3/position_unit) { // test12
                fprintf(binout, "Break group: hyperbolic escape!\n\t");
                fprintf(binout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                fprintf(binout, "N_member: %d\n\t", n_member);
                fprintf(binout, "separation: %e pc\n\t", bin_root.r*position_unit); 
                fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
                fprintf(binout, "ecca: %e\n\t", bin_root.ecca);
                fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(binout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                fflush(binout);
                return true;
            }
            /*
            if (n_member == 2) {
                if (bin_root.r > 1e-3/position_unit) { // original
                    // if (bin_root.r > 2e-3/position_unit) { // test
                    fprintf(binout, "Break group: hyperbolic escape!\n\t");
                    fprintf(binout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                    fprintf(binout, "N_member: %d\n\t", n_member);
                    fprintf(binout, "separation: %e pc\n\t", bin_root.r*position_unit); 
                    fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
                    fprintf(binout, "ecca: %e\n\t", bin_root.ecca);
                    fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                    fprintf(binout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                    fflush(binout);
                    return true;
                }
            }
            else {
                if (bin_root.r > 2e-3/position_unit) { // original
                    // if (bin_root.r > 2e-3/position_unit) { // test
                    fprintf(binout, "Break group: hyperbolic escape!\n\t");
                    fprintf(binout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                    fprintf(binout, "N_member: %d\n\t", n_member);
                    fprintf(binout, "separation: %e pc\n\t", bin_root.r*position_unit); 
                    fprintf(binout, "ecc: %e\n\t", bin_root.ecc);
                    fprintf(binout, "ecca: %e\n\t", bin_root.ecca);
                    fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                    fprintf(binout, "r_crit: %e pc\n", sym_int.info.r_break_crit*position_unit);
                    fflush(binout);
                    return true;
                }
            }
            */
/*
            // check for next step
            Float dr = drdv/bin_root.r*sym_int.particles.cm.TimeStepIrr*EnzoTimeStep;
            Float rp = dr  + bin_root.r;
            if (rp > sym_int.info.r_break_crit) {
                // in case r is too small, avoid too early quit of group
                Float rph = 0.5*dr + bin_root.r;
                if ( rph < sym_int.info.r_break_crit && bin_root.r <0.2*sym_int.info.r_break_crit) {
                    // std::cerr<<"Hyperbolic will escape but dr is too small, reduce cm step by half first, time: "<<CurrentTime<<" N_member: "<<n_member<<" drdv: "<<drdv<<" separation : "<<bin_root.r<<" r_pred: "<<rp<<" drdv: "<<drdv<<" dt: "<<sym_int.particles.cm.TimeStepIrr<<" r_crit: "<<sym_int.info.r_break_crit<<std::endl;
                    sym_int.particles.cm.TimeLevelIrr--;
                    sym_int.particles.cm.TimeStepIrr = static_cast<REAL>(pow(2, sym_int.particles.cm.TimeLevelIrr));
                    sym_int.particles.cm.TimeBlockIrr = static_cast<ULL>(pow(2, sym_int.particles.cm.TimeLevelIrr-time_block));
                
                }
                else {
                    // std::cerr<<"Break group: hyperbolic will escape, time: "<<CurrentTime<<" N_member: "<<n_member<<" drdv: "<<drdv<<" separation : "<<bin_root.r<<" r_pred: "<<rp<<" drdv: "<<drdv<<" dt: "<<sym_int.particles.cm.TimeStepIrr<<" r_crit: "<<sym_int.info.r_break_crit<<std::endl;
                    fprintf(binout, "Break group: hyperbolic will escape!\n\t time: %e Myr\n\t N_member: %d\n\t separation: %e pc\n\t r_crit: %e pc\n", CurrentTime*EnzoTimeStep*1e4, n_member, bin_root.r*position_unit, sym_int.info.r_break_crit*position_unit);
                    fflush(binout);
                    return true;
                }
                return false;
            }
*/
        }

    }
    // return false;
// /* // test_1e4_2
    // check strong perturbed binary case (only check further if it is outgoing case)
    // calculate slowdown in a consistent way like in checknewgroup to avoid switching
    // fcm may not properly represent the perturbation force (perturber mass is unknown)
    if (outgoing_flag && bin_root.semi > 0.0 && n_member == 2) {

        AR::SlowDown sd;
        auto& sd_group = sym_int.info.getBinaryTreeRoot().slowdown;
        sd.initialSlowDownReference(sd_group.getSlowDownFactorReference(),sd_group.getSlowDownFactorMax());
        sd.timescale = sd_group.timescale;
        sd.period = sd_group.period;

        // sd.pert_in = manager.interaction.calcPertFromBinary(bin_root); // original
        sd.pert_in = manager.interaction.calcPertFromMR(bin_root.r, bin_root.m1, bin_root.m2);
        Float acc_cm[3];
        for (int i = 0; i < Dim; ++i) {
            acc_cm[i] = sym_int.particles.cm.a_irr[i][0]; // a_irr or a_tot?
        }
        Float fcm[3] = {acc_cm[0]*bin_root.Mass, acc_cm[1]*bin_root.Mass, acc_cm[2]*bin_root.Mass};
        sd.pert_out= manager.interaction.calcPertFromForce(fcm, bin_root.Mass, bin_root.Mass);
        sd.calcSlowDownFactor();
        Float kappa_org = sd.getSlowDownFactorOrigin();

        if (kappa_org<kappa_org_crit) {
            // in binary case, only break when apo is larger than distance criterion
            Float apo = bin_root.semi * (1.0 + bin_root.ecc);
            // if (apo>sym_int.info.r_break_crit||bin_root.semi<0.0) {
            if (apo>sym_int.info.r_break_crit && bin_root.r > 2 * rbin/position_unit) { // test8 // fiducial
            // if ((apo>sym_int.info.r_break_crit && bin_root.r > 1.2e-3/position_unit)||bin_root.semi<0.0) { // test12
                auto& sd_root = sym_int.info.getBinaryTreeRoot().slowdown;

                fprintf(binout, "Break group: strong perturbed!\n\t");
                fprintf(binout, "time: %e Myr\n\t", CurrentTime*EnzoTimeStep*1e4);
                fprintf(binout, "N_member: %d\n\t", n_member);
                fprintf(binout, "pert_in: %e \n\t", sd_root.pert_in);
                fprintf(binout, "pert_out: %e \n\t", sd_root.pert_out);
                fprintf(binout, "kappa_org: %e \n\t", kappa_org);
                fprintf(binout, "separation: %e pc\n\t", bin_root.r*position_unit);
                fprintf(binout, "semi: %e pc\n\t", bin_root.semi*position_unit);
                fprintf(binout, "ecc: %e \n\t", bin_root.ecc);
                fprintf(binout, "ecca: %e \n\t", bin_root.ecca);
                fprintf(binout, "peri: %e pc\n\t", bin_root.semi*(1-bin_root.ecc)*position_unit);
                fprintf(binout, "apo: %e pc\n\t", bin_root.semi*(1+bin_root.ecc)*position_unit);
                fprintf(binout, "r_break: %e pc\n", sym_int.info.r_break_crit*position_unit);
                fflush(binout);
                return true;
            }
        }
    }
// */ // test_1e4_2
    return false;

}
#endif