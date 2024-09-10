#pragma once

#include "Common/binary_tree.h"
#include "tidal_tensor.hpp"
#include "pseudoparticle_multipole.hpp"
#include <cassert>


//! class to store necessary information for using artificial particles
/*!
                  single    c.m.             members          initial     artificial
      mass_backup 0         mass      (+)    mass     (+)     -LARGE       default is 0.0 / data stored (-/0)
      status      0         n_members (+)    c.m. adr (-)     -LARGE       position in artificial particle array / data stored (+)
 */
class ArtificialParticleInformation{
private:
    REAL mass_backup;
    REAL status;

public:
    //! initializer
    ArtificialParticleInformation(): mass_backup(-std::numeric_limits<float>::max()), status(-std::numeric_limits<float>::max()) {}

    //! set particle type to member
    /*! @param[in] _mass: mass to backup
     */
    void setParticleTypeToMember(const REAL _mass = std::numeric_limits<float>::max(), const REAL _status = -std::numeric_limits<float>::max()) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_status<0.0);
#endif
        mass_backup = _mass;
        status = _status;
    };

    //! return whether the particle type is member
    bool isMember() const {
        return (status<0.0);
    }

    //! set particle type to artificial
    /*! @param[in] _status: status to save
     */
    void setParticleTypeToCM(const REAL _mass_backup = std::numeric_limits<float>::max(), const REAL _status=std::numeric_limits<float>::max()) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_status>0.0&&_mass_backup>0.0);
#endif
        status = _status;
        mass_backup = _mass_backup;
    };

    //! return whether the particle type is c.m.
    bool isCM() const {
        return (status>0.0 && mass_backup>0.0);
    }

    //! set backup mass
    void setMassBackup(const REAL _mass) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert((isMember()||isCM())&&_mass>0.0);
#endif
        mass_backup = _mass;
    }

    //! get backup mass
    REAL getMassBackup() const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(isMember()||isCM());
#endif
        return mass_backup;
    }

    //! set particle type to single
    void setParticleTypeToSingle() {
        mass_backup = 0.0;
        status = 0.0;
    };

    //! return whether the particle type is single
    bool isSingle() const {
        return (status==0.0 && mass_backup==0.0);
    }

    //! set particle type to artificial
    /*! @param[in] _status: status to save
     */
    void setParticleTypeToArtificial(const REAL _status=std::numeric_limits<float>::max(), const REAL _mass_backup = 0.0) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(_status>0.0);
#endif
        status = _status;
        mass_backup = _mass_backup;
    };

    //! return whether the particle type is artificial
    bool isArtificial() const {
        return (status>0.0);
    }

    //! store one data (only in artificial particles)
    /*! positive data stored in status, otherwise in mass_backup;
     */
    void storeData(const REAL _data) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(isArtificial());
#endif
        if (_data>0.0) status = _data;
        else      mass_backup = _data;
    }

    //! get stored data (only in artificial particles)
    /*! @param[in] is_positive: indicate whether stored data is positive (true) or negative (false)
     */
    REAL getData(const bool is_positive) const {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert(isArtificial());
#endif
        if (is_positive) return status;
        else             return mass_backup;
    }

    //! set status
    void setStatus(const REAL _status) {
#ifdef ARTIFICIAL_PARTICLE_DEBUG
        assert((isMember()&&_status<0.0)||(isArtificial()&&_status>0.0));
#endif        
        status = _status;
    }

    //! get status
    REAL getStatus() const {
        return status;
    }

    //! return whether the particle is unused
    bool isUnused() const {
        return (status<0.0 && mass_backup<0.0);
    }

    //! set particle type to unused
    void setParticleTypeToUnused() {
        status = - NUMERIC_FLOAT_MAX;
        mass_backup = - NUMERIC_FLOAT_MAX;
    }

    //! print titles of class members using column style
    /*! print titles of class members in one line for column style
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    static void printColumnTitle(std::ostream & _fout, const int _width=20) {
        _fout<<std::setw(_width)<<"mass_bk"
             <<std::setw(_width)<<"status";
    }

    //! print column title with meaning (each line for one column)
    /*! @param[out] _fout: std::ostream output object
      @param[in] _counter: offset of the number counter for each line to indicate the column index (defaulted 0)
      @param[in] _offset: the printing whitespace offset for each line (defaulted 0)
      \return: the total counter of columns
     */
    static int printTitleWithMeaning(std::ostream & _fout, const int _counter=0, const int _offset=0) {
        int counter = _counter;
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". mass_bk: artificial particle parameter 1 [formatted] (0.0)\n";
        counter++;
        _fout<<std::setw(_offset)<<" "<<counter<<". status: artificial particle parameter 2 [formatted] (0.0)\n";
        return counter;
    }

    //! print data of class members using column style
    /*! print data of class members in one line for column style. Notice no newline is printed at the end
      @param[out] _fout: std::ostream output object
      @param[in] _width: print width (defaulted 20)
     */
    void printColumn(std::ostream & _fout, const int _width=20){
        _fout<<std::setw(_width)<<mass_backup
             <<std::setw(_width)<<status;
    }

    //! write class data with ASCII format
    /*! @param[in] _fout: file IO for write
     */
    void writeAscii(FILE* _fout) const{
        fprintf(_fout, "%26.17e %26.17e ", 
                this->mass_backup, this->status);
    }

    //! write class data with BINARY format
    /*! @param[in] _fout: file IO for write
     */
    void writeBinary(FILE* _fin) const{
        fwrite(this, sizeof(ArtificialParticleInformation), 1, _fin);
    }

    //! read class data with ASCII format
    /*! @param[in] _fin: file IO for read
     */
    void readAscii(FILE* _fin) {
        long long int rcount=fscanf(_fin, "%lf %lf ",
                              &this->mass_backup, &this->status);
        if (rcount<2) {
            std::cerr<<"Error: Data reading fails! requiring data number is 2, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! read class data with BINARY format
    /*! @param[in] _fin: file IO for read
     */
    void readBinary(FILE* _fin) {
        size_t rcount = fread(this, sizeof(ArtificialParticleInformation), 1, _fin);
        if (rcount<1) {
            std::cerr<<"Error: Data reading fails! requiring data number is 1, only obtain "<<rcount<<".\n";
            abort();
        }
    }

    //! print parameters
    void print(std::ostream & _fout) const{
        _fout<<" mass_bk="<<mass_backup
             <<" status="<<status;
    }    
    
};