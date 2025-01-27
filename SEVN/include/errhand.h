//
// Created by mario on 07/11/18.
//

#ifndef SEVN_REVISED_ERRHAND_H
#define SEVN_REVISED_ERRHAND_H


#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <string>
#include <exception>
#include <sevnlog.h>


//NB, THis class is now deprecated since all the logging and errors are managed by SevnLogging.
// we mantain it  for now for compatibility but all the all the old call to sevenerr have to be repalced
//TODO the namespade should be sevnerr
namespace sevnstd {

class sevnerr : public std::exception {

public:

    explicit sevnerr(std::string s="") {
        set_mess(s);
    }

    inline  sevnerr istance(std::string& s){
        return sevnerr(s);
    }

    const char * what () const throw(){
        return &mess[0];
    }

    virtual inline  std::string default_mess() const {return "SEVN error: ";}


protected:

    void inline set_mess(std::string& s){
        mess=default_mess()+s;
    }
    std::string mess;

};


class sevnio_error : public sevnerr {

    public:
        explicit sevnio_error(std::string s="") {
            set_mess(s);
        }
        inline sevnio_error istance(std::string& s){
            return sevnio_error(s);
        }
        inline  std::string default_mess() const override {return "SEVN IO error: ";}
};


class sse_error : public sevnerr {

    public:
        explicit sse_error(std::string s="") {
            set_mess(s);
        }
        inline sse_error istance(std::string& s){
            return sse_error(s);
        }
        inline  std::string default_mess() const override {return "SSE error: ";}
    };


class bse_error : public sevnerr {

    public:
        explicit bse_error(std::string s="") {
            set_mess(s);
        }
        inline bse_error istance(std::string& s){
            return bse_error(s);
        }
        inline  std::string default_mess() const override {return "BSE error: ";}
};


class rl_error : public bse_error {

    public:
        explicit rl_error(std::string s="") {
            set_mess(s);
        }
        inline rl_error istance(std::string& s){
            return rl_error(s);
        }
        inline  std::string default_mess() const override {return "RocheLobe error: ";}
};

class ce_error : public bse_error {

public:
    explicit ce_error(std::string s="") {
        set_mess(s);
    }

    inline ce_error istance(std::string& s){
        return ce_error(s);
    }
    inline  std::string default_mess() const override {return "CommonEnvelope error: ";}
};

class sn_error : public sse_error {

public:
    explicit sn_error(std::string s=""){
        set_mess(s);
    }

    inline sn_error istance(std::string& s){
        return sn_error(s);
    }

    inline  std::string default_mess() const override {return "SN error: ";}

};

class jtrack_error : public sevnerr {

    public:
        explicit jtrack_error(std::string s="") {
        set_mess(s);
    }
        inline jtrack_error istance(std::string& s){
            return jtrack_error(s);
        }
        inline std::string default_mess() const override { return "Change tracks error: ";}
};

class params_error : public sevnerr {

public:
    explicit params_error(std::string s=""){
        set_mess(s);
    }
    inline params_error istance(std::string& s){
        return params_error(s);
    }
    inline std::string default_mess() const override { return "SEVN parameters error: ";}
};

class notimplemented_error : public  sevnerr {
public:
    explicit notimplemented_error(std::string s=""){
        set_mess(s);
    }
    inline notimplemented_error istance(std::string& s){
        return notimplemented_error(s);
    }
    inline std::string default_mess() const override { return "SEVN not implemented error: ";}
};

class sanity_error : public  sevnerr {
public:
    explicit sanity_error(std::string s=""){
        set_mess(s);
    }
    inline sanity_error istance(std::string& s){
        return sanity_error(s);
    }
    inline std::string default_mess() const override { return "SEVN sanity check error: ";}
};


}


//This is an example of how override throw of exception class, maybe at a certain point sevnerr can inerhit from exception
/*
class SevnException : public std::exception
{
public:
    const char * what() const throw () override  {
        return "SEVN EXCEPTION";
    }
};
 */

#endif //SEVN_REVISED_ERRHAND_H
