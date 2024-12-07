#pragma once

#include "Common/Float.h"

namespace AR{
    
    //! Slow-down parameter control class
    /*! Determine the slow-down factor due to the perturbation and internal force
      \f$ \kappa = k_0 / [F_{pert,max}/F_{inner}] \f$
    */
    class SlowDown{
    private:
        Float time_update_;    // update time (int) for new slowdown factor
        Float kappa_;          // slow-down factor
        Float kappa_org_;      // original slow-down factor without kappa limit (1.0, kappa_max)
        Float kappa_max_;      // maximum kappa factor
        Float kappa_ref_;      ///< reference kappa factor; slow-down factor kappa = max(1,kref/perturbation_factor)
        Float timescale_max_;  ///< maximum timescale 
    public:
        Float pert_in;           // inner strength
        Float pert_out;          // perturbation strength
        Float timescale;         // perturbation strength change timescale
        Float period;            // orbital period

        //! defaulted constructor
        SlowDown(): time_update_(Float(0.0)), kappa_(Float(1.0)), kappa_org_(Float(1.0)), kappa_max_(Float(1.0)), kappa_ref_(Float(1.0e-6)), timescale_max_(NUMERIC_FLOAT_MAX), pert_in(0.0), pert_out(0.0), timescale(NUMERIC_FLOAT_MAX), period(NUMERIC_FLOAT_MAX) {}
    
        //! clear function
        void clear(){
            time_update_ = Float(0.0);
            kappa_ = kappa_org_ = kappa_max_ = Float(1.0);
            kappa_ref_ = Float(1.0e-6);
            pert_in = pert_out = Float(0.0);
            timescale_max_ = timescale = period = NUMERIC_FLOAT_MAX;
        }

        //! initialize slow-down parameters
        /*! Set slow-down parameters, slow-down method will be switched on
          @param[in] _kappa_ref: reference kappa factor; slow-down factor kappa = max(1,kref/perturbation_factor)
          @param[in] _timescale_max: maximum timescale 
        */
        void initialSlowDownReference(const Float _kappa_ref, const Float _timescale_max) {
            ASSERT(_kappa_ref>0.0);
            ASSERT(_timescale_max>0.0);
            kappa_ref_ = _kappa_ref;
            timescale_max_ = _timescale_max;
        }

        //! set update time for new slowdown factor
        void setUpdateTime (const Float _time) {
            time_update_ = _time;
        }

        //! get update time for new slowdown factor
        Float getUpdateTime() const {
            return time_update_;
        }

        //! advance update time by one period
        void increaseUpdateTimeOnePeriod() {
            time_update_ += period;
        }

        //! manually set kappa
        void setSlowDownFactor(const Float _kappa) {
            kappa_ = _kappa;
            kappa_org_ = _kappa;
        }

        //! calculate slowdown factor based on perturbation and inner acceleration
        /* if it is a hyperbolic encounter, (ebin_>0), set slowdown factor to 1.0
          \return slowdown factor
         */
        Float calcSlowDownFactor() {
            // hyberbolic case or no perturbation
            if(pert_in<=0.0) {
                kappa_org_ = kappa_max_ = kappa_ = Float(1.0);
            }
            else { 
                ASSERT(period>0.0);
                kappa_max_ = std::max(Float(1.0), timescale/period);
                if (pert_out==0.0) kappa_org_ = kappa_max_;
                else kappa_org_ = kappa_ref_*pert_in/pert_out;
                kappa_ = std::min(kappa_org_, kappa_max_);
                kappa_ = std::max(Float(1.0), kappa_);
            }
            return kappa_;
        }

        //! Get slow-down factor
        /*!
          \return get adjusted kappa by keeping phase corrected
        */
        Float getSlowDownFactor() const {
            return kappa_;
        }

        //! Get original slow-down factor
        /*!
          \return kappa_origin
        */
        Float getSlowDownFactorOrigin() const {
            return kappa_org_;
        }

        //! Get sd reference factor
        Float getSlowDownFactorReference() const {
            return kappa_ref_;
        }

        //! Get slow-down fact maximum
        Float getSlowDownFactorMax() const {
            return kappa_max_;
        }

        Float getPertIn() const {
            return pert_in;
        }

        Float getPertOut() const {
            return pert_out;
        }

        Float getTimescaleMax() const {
            return timescale_max_;
        }

        //! write class data with BINARY format
        /*! @param[in] _fout: file IO for write
         */
        void writeBinary(FILE *_fout) {
            fwrite(this, sizeof(*this),1,_fout);
        }

        //! read class data with BINARY format
        /*! @param[in] _fin: file IO for read
         */
        void readBinary(FILE *_fin) {
            size_t rcount = fread(this, sizeof(*this), 1, _fin);
            if (rcount<1) {
                std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
                abort();
            }
        }
        
        //! get backup data size 
        /*! \return the data array size for backupSlowDownFactorAndTimeReal()
         */
        static int getBackupDataSize() {
            return 2;
        }

        //! backup real time and force ratio
        /*! @param[in] _bk: backup data array, should be size of getBackupDataSize() (2)
          \return backup array size
         */
        int backup(Float* _bk) {
            _bk[0] = kappa_;
            _bk[1] = time_update_;
            return 2;
        }

        //! restore real time and force ratio
        /*! @param[in] _bk: restore data array[3]
          \return backup array size
         */
        int restore(Float* _bk) {
            kappa_       =   _bk[0];
            time_update_ =   _bk[1];
            return 2;
        }
    
        //! print slowdown data
        /*! Print slowdown data 
          @param[in] fout: ofstream for printing
          @param[in] precision: printed precision for one variable
          @param[in] width: printing width for one variable
        */
        void print(std::ostream & fout, const int precision=15, const int width=23) {
            ASSERT(width>0);
            ASSERT(precision>0);
            fout<<"kappa= "<<std::setw(width)<<kappa_
                <<"kappa_org= "<<std::setw(width)<<kappa_org_
                <<"kappa_max= "<<std::setw(width)<<kappa_max_
                <<"kappa_ref= "<<std::setw(width)<<kappa_ref_;
        }
    
        //! print titles of class members using column style
        /*! print titles of class members in one line for column style
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        static void printColumnTitle(std::ostream & _fout, const int _width=20) {
            _fout<<std::setw(_width)<<"SD_factor"
                 <<std::setw(_width)<<"SD_factor_org"
                 <<std::setw(_width)<<"SD_factor_max";
        }

        //! print data of class members using column style
        /*! print data of class members in one line for column style. Notice no newline is printed at the end
          @param[out] _fout: std::ostream output object
          @param[in] _width: print width (defaulted 20)
        */
        void printColumn(std::ostream & _fout, const int _width=20){
            _fout<<std::setw(_width)<<kappa_
                 <<std::setw(_width)<<kappa_org_
                 <<std::setw(_width)<<kappa_max_;
        }

    };
}
