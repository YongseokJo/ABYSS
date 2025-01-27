//
// Created by Giuliano Iorio on 2020-01-25.
//


#include <sevnlog.h>
#include <iostream>
#include <string>
#include <stdexcept>
#include <errhand.h>


//Initialise the logging counter
unsigned int sevnstd::SevnLogging::count_debug=0;
unsigned int sevnstd::SevnLogging::count_info=0;
unsigned int sevnstd::SevnLogging::count_warning=0;
unsigned int sevnstd::SevnLogging::count_error=0;
unsigned int sevnstd::SevnLogging::count_custom_log=0;


//Set the default to debug level if debug is defined or INFO in the other case.
#ifdef DEBUG
    int sevnstd::SevnLogging::log_level = 10;
#else
    int sevnstd::SevnLogging::log_level = 20;
#endif


void sevnstd::SevnLogging::set_level(std::string level)   {

    unsigned int lvl;

    if(level=="debug")
        lvl=_LOG_LEVEL::_debug;
    else if(level=="info")
        lvl=_LOG_LEVEL::_info;
    else if(level=="warning")
        lvl=_LOG_LEVEL::_warning;
    else if(level=="error")
        lvl=_LOG_LEVEL::_error;
    else if(level=="critical")
        lvl=_LOG_LEVEL::_critical;
    else
        critical("Loglevel "+level+"not allowed. Possible options are: debug, info, warning, error.",
                __FILE__,__LINE__,sevnstd::params_error());

    set_level(lvl);

}

void sevnstd::SevnLogging::log(int level, std::string errstate, const char *file_input, int line_input, int stop)  const{



    if (level>=log_level) {

        #ifdef _OPENMP
                int num_thread=omp_get_thread_num();
        #else
                int num_thread=0;
        #endif


        std::ostringstream oss;
        oss << " LOG::LOG (lvl "<<level <<", Thread " << num_thread << "): " << std::endl;
        oss << "   Message : " << errstate << std::endl;
        if (file_input)
            oss << " From file: " << std::string(file_input) << std::endl;
        if (line_input >= 0)
            oss << " From line: " << line_input << std::endl;
        if (stop)
            throw std::runtime_error(oss.str());
        else
            std::cout << oss.str();
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        count_custom_log++;
    }

}



void sevnstd::SevnLogging::critical(std::string errstate, const char *file_input, int line_input)  const{

    #ifdef _OPENMP
        int num_thread=omp_get_thread_num();
    #else
        int num_thread=0;
    #endif

    std::ostringstream oss;
    oss << " LOG::CRITICAL (Thread " << num_thread<< "): " << std::endl;
    oss << "   Message : " << errstate << std::endl;
    if (file_input)
        oss << " From file: " << std::string(file_input) << std::endl;
    if (line_input >= 0)
        oss << " From line: " << line_input << std::endl;

    throw sevnstd::sevnerr(oss.str());

}

void sevnstd::SevnLogging::error(std::string errstate, const char *file_input, int line_input , bool stop)  const{


    if (_LOG_LEVEL::_error>=log_level) {

        #ifdef _OPENMP
                int num_thread=omp_get_thread_num();
        #else
                int num_thread=0;
        #endif


        std::ostringstream oss;
        oss << " LOG::ERROR (Thread " << num_thread << "): " << std::endl;
        oss << "   Message: " << errstate << std::endl;
        if (file_input)
            oss << " From file: " << std::string(file_input) << std::endl;
        if (line_input >= 0)
            oss << " From line: " << line_input << std::endl;
        if (stop)
            throw sevnstd::sevnerr(oss.str());
        else
            std::cerr << oss.str();
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        count_error++;
    }

}


void sevnstd::SevnLogging::warning(std::string errstate, const char *file_input, int line_input)  const{


    //if (_LOG_LEVEL::_warning>=log_level and count_warning<=MAX_N_WARNING) {
    if (_LOG_LEVEL::_warning>=log_level) {

        #ifdef _OPENMP
                int num_thread=omp_get_thread_num();
        #else
                int num_thread=0;
        #endif

        std::ostringstream oss;
        oss << " LOG::WARNING (Thread " << num_thread << "): " << std::endl;
        oss << "   Message: " << errstate << std::endl;
        if (file_input)
            oss << " From file: " << std::string(file_input) << std::endl;
        if (line_input >= 0)
            oss << " From line: " << line_input << std::endl;
        std::cerr << oss.str();
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        count_warning++;
    }

}


void sevnstd::SevnLogging::info(std::string errstate, const char *file_input, int line_input)  const{


    if (_LOG_LEVEL::_info>=log_level) {

        #ifdef _OPENMP
                int num_thread=omp_get_thread_num();
        #else
                int num_thread=0;
        #endif

        std::ostringstream oss;
        oss << " LOG::INFO (Thread " << num_thread << "): " << std::endl;
        oss << "   " << errstate << std::endl;
        if (file_input)
            oss << " From file: " << std::string(file_input) << std::endl;
        if (line_input >= 0)
            oss << " From line: " << line_input << std::endl;
        std::cout << oss.str();
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        count_info++;
    }

}


void sevnstd::SevnLogging::debug(std::string errstate, const char *file_input, int line_input)  const{


    if (_LOG_LEVEL::_debug>=log_level) {

        #ifdef _OPENMP
                int num_thread=omp_get_thread_num();
        #else
                int num_thread=0;
        #endif

        std::ostringstream oss;
        oss << " LOG::DEBUG (Thread " << num_thread<< "): " << std::endl;
        oss << "   Message: " << errstate << std::endl;
        if (file_input)
            oss << " From file: " << std::string(file_input) << std::endl;
        if (line_input >= 0)
            oss << " From line: " << line_input << std::endl;
        std::cout << oss.str();
        #ifdef _OPENMP
        #pragma omp atomic
        #endif
        count_debug++;
    }

}
