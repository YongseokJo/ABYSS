//
// Created by Giuliano on 29/11/19.
// Header to store logging functions and definitions (mostly for debug)
//

#ifndef SEVN_SEVNLOG_H
#define SEVN_SEVNLOG_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <string>
#include <sstream>
#include <iostream>
#include "errhand.h"


//TODO Remove this and correctly use SevnLogging
//GI291119: Define the DEBUG_LOG functions to print the Debug message only if enable in the compilation
#ifdef DEBUG
#define DEBUG_LOG(str) do { std::cout <<  "DEBUG: FILE::" << __FILE__ << " LINE::"  <<__LINE__ << std::endl << " -> " << str << " <- " << std::endl; } while (false)
#else
#define DEBUG_LOG(str) do {  } while (false)
#endif

//TODO Add the possibility to directly flush the output in some log file(s).
//TODO It is really thread safe?
namespace sevnstd{

    class sevnerr;

    /*!  A thread safe(really?) Logging class to handle the message output and the error. It is
     * based on a log level scheme. There is a static attribute log_level that has
     * some  integer value by default (20 or 10 id DEBUG has been enabled).
     * Then a given message is sent in output only if its level is larger than log_level.
     * A new log_level can be set with the method set_level.
     * The possible logging message are debug(lvl 10), info(lvl 20), warning(lvl 30),
     * error(lvl 40), critical (no level always printed). Critical raises automatically an exception, while
     * in error is optional. A general log method can be used to print  output with a custom level.
     * */

    class SevnLogging {

    public:

        /**
         * Default class constructor.
         */
        SevnLogging() {}

        //TODO Is ok for a single instance to modify the log_level. Is maybe better to modify only a local log_level and the use this in the various log function?
        /**
         * Class constructor that set the static log level attribute.
         * @param level log level to set.
         */
        SevnLogging(int level) {
            set_level(level);
        }
        /**
        * Default class destructor.
        */
        ~SevnLogging() {}

        /**
         * Logs a message to std::cout with integer level  on this logger.
         * @param level   message level.
         * @param errstate  message to log.
         * @param file_input Null or __FILE__. if __FILE__ is used the message logs the name of the file where the log is called.
         * @param line_input Null or __LINE__. if __LINE__ is used the message logs the row number in the file where the log is called.
         * @param stop if != 0, throw a runtime exception.
         */
        void log(int level, std::string errstate, const char *file_input = nullptr, int line_input = -1, int stop = 0)  const;

        /**
         * Logs a message to std::cout with level DEBUG (lvl 10) on this logger.
         * @param errstate message to log.
         * @param file_input Null or __FILE__. if __FILE__ is used the message logs the name of the file where the log is called.
         * @param line_input Null or __LINE__. if __LINE__ is used the message logs the row number in the file where the log is called.
         */
        void debug(std::string errstate, const char *file_input = nullptr, int line_input = -1)  const;


        /**
         * Logs a message to std::cout with level INFO (lvl 20) on this logger.
         * @param errstate message to log.
         * @param file_input Null or __FILE__. if __FILE__ is used the message logs the name of the file where the log is called.
         * @param line_input Null or __LINE__. if __LINE__ is used the message logs the row number in the file where the log is called.
         */
        void info(std::string errstate, const char *file_input = nullptr, int line_input = -1)  const;


        /**
         * Logs a message with std::cerr level WARNING (lvl 30) on this logger.
         * @param errstate message to log.
         * @param file_input Null or __FILE__. if __FILE__ is used the message logs the name of the file where the log is called.
         * @param line_input Null or __LINE__. if __LINE__ is used the message logs the row number in the file where the log is called.
         */
        void warning(std::string errstate, const char *file_input = nullptr, int line_input = -1)  const;

        /**
         * Logs a message with std::cerr level ERROR (lvl 40) on this logger and throw an exception.
         * @param errstate message to log.
         * @param file_input Null or __FILE__. if __FILE__ is used the message logs the name of the file where the log is called.
         * @param line_input Null or __LINE__. if __LINE__ is used the message logs the row number in the file where the log is called.
         * @param stop if true throw a sevnerr exception.
         */
        void error(std::string errstate, const char *file_input = nullptr, int line_input = -1, bool stop = true)  const;

        /**
         *
         * @tparam E exception derived from the class sevnerr
         * @param errstate message to log.
         * @param file_input file_input Null or __FILE__. if __FILE__ is used the message logs the name of the file where the log is called.
         * @param line_input line_input Null or __LINE__. if __LINE__ is used the message logs the row number in the file where the log is called.
         * @param stop if if true throw a err exception (see below).
         * @param err exception derived from the class sevnerr
         */
        template<class E>
        void error(std::string errstate, const char *file_input = nullptr, int line_input = -1, bool stop = true, E&& err= nullptr)  const{

            #ifdef _OPENMP
                        int num_thread=omp_get_thread_num();
            #else
                        int num_thread=0;
            #endif

            std::ostringstream oss;
            oss << " LOG::ERROR (Thread " << num_thread << "): " << std::endl;
            oss << "   Message : " << errstate << std::endl;
            if (file_input)
                oss << " From file: " << std::string(file_input) << std::endl;
            if (line_input >= 0)
                oss << " From line: " << line_input << std::endl;
            std::string err_mess=oss.str();
            if (stop)
                throw err.istance(err_mess);
            else
                std::cerr << oss.str();

#ifdef _OPENMP
#pragma omp atomic
#endif
            count_error++;
        }


        /**
         * Logs a message with std::cerr level CRITICAL (no level) on this logger and throw an sevnerr exception.
         * This message will never be filtered out.
         * @param errstate message to log.
         * @param file_input Null or __FILE__. if __FILE__ is used the message logs the name of the file where the log is called.
         * @param line_input Null or __LINE__. if __LINE__ is used the message logs the row number in the file where the log is called.
         */
        void critical(std::string errstate, const char *file_input = nullptr, int line_input = -1) const;


        /**
         * Logs a message with std::cerr level CRITICAL (no level) on this logger and throw an  exception E.
         * @tparam E exception derived from the class sevnerr
         * @param errstate message to log.
         * @param file_input file_input Null or __FILE__. if __FILE__ is used the message logs the name of the file where the log is called.
         * @param line_input line_input Null or __LINE__. if __LINE__ is used the message logs the row number in the file where the log is called.
         * @param err exception derived from the class sevnerr
         */
        template<class E>
        void critical(std::string errstate, const char *file_input = nullptr, int line_input = -1, E&& err= nullptr) const{

            #ifdef _OPENMP
                        int num_thread=omp_get_thread_num();
            #else
                        int num_thread=0;
            #endif

            std::ostringstream oss;
            oss << " LOG::CRITICAL (Thread " << num_thread << "): " << std::endl;
            oss << "   Message : " << errstate << std::endl;
            if (file_input)
                oss << " From file: " << std::string(file_input) << std::endl;
            if (line_input >= 0)
                oss << " From line: " << line_input << std::endl;
            std::string err_mess=oss.str();
            throw err.istance(err_mess);
        }




        ///Variadic prints
        //debug
        void inline pdebug()  const {

            #ifdef _OPENMP
                        int num_thread=omp_get_thread_num();
            #else
                        int num_thread=0;
            #endif

            std::cout<<"\nLOG::DEBUG (Thread " << num_thread << ")"<< std::endl;
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            count_debug++;}

        template<typename T, typename... Tail>
        void pdebug(T head, Tail... tail)  const{
            if (_LOG_LEVEL::_debug>=log_level) {
                std::cout << head << " ";
                pdebug(tail...);
            }
        }


        //info
        void inline pinfo()  const {

            #ifdef _OPENMP
                        int num_thread=omp_get_thread_num();
            #else
                        int num_thread=0;
            #endif

            std::cout<<"\nLOG::INFO (Thread " << num_thread << ")"<< std::endl;
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            count_info++;}

        template<typename T, typename... Tail>
        void pinfo(T head, Tail... tail)  const {
            if (_LOG_LEVEL::_info>=log_level) {
                std::cout << head << " ";
                pinfo(tail...);
            }
        }

        //warning
        void inline pwarning()  const {


            #ifdef _OPENMP
                        int num_thread=omp_get_thread_num();
            #else
                        int num_thread=0;
            #endif

            std::cerr<<"\nLOG::WARNING (Thread " << num_thread << ")"<< std::endl;
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            count_warning++;
        }

        template<typename T, typename... Tail>
        void pwarning(T head, Tail... tail)  const {
            //if (_LOG_LEVEL::_warning>=log_level and count_warning<=MAX_N_WARNING) {
            if (_LOG_LEVEL::_warning>=log_level) {
                std::cerr << head << " ";
                pwarning(tail...);
            }
        }


        ////WARNING C++17 feature
/*        *//**
         * Logs a message to std::cout with level DEBUG (lvl 10) on this logger.
         * @tparam Args  pack of Variadic arguments
         * @param args args to be printed
         *//*
        template<typename... Args>
        void pdebug(Args... args){
            if (_LOG_LEVEL::_debug>=log_level){
                std::cout << " LOG::DEBUG (Thread " << omp_get_thread_num() << "): " << std::endl;
                std::cout << "  Message:";
                ((std::cout << " "<<args), ...);
                #pragma omp atomic
                count_debug++;
            }
        }

        *//**
         * Logs a message to std::cout with level DEBUG (lvl 10) on this logger.
         * @tparam Args  pack of Variadic arguments
         * @param args args to be printed
         *//*
        template<typename... Args>
        void pinfo(Args... args){
            if (_LOG_LEVEL::_info>=log_level){
                std::cout << " LOG::INFO (Thread " << omp_get_thread_num() << "): " << std::endl;
                std::cout << "  Message:";
                ((std::cout << " "<<args), ...);
                #pragma omp atomic
                count_info++;
            }
        }*/






        /**
         * Get the current level of this logger
         * @return log level.
         */
        inline int get_level() { return log_level;};

        /**
        * Get the current counter of debug calls
        * @return current counter of debug calls
        */
        inline unsigned  int get_Ndebug() { return count_debug;};
        /**
        * Get the current counter of info calls
        * @return current counter of info calls
        */
        inline unsigned  int get_Ninfo() { return count_info;};
        /**
        * Get the current counter of warning calls
        * @return current counter of warning calls
        */
        inline unsigned  int get_Nwarning() { return count_warning;};
        /**
        * Get the current counter of error calls
        * @return current counter of error calls
        */
        inline unsigned  int get_Nerror() { return count_error;};
        /**
        * Get the current counter of custom log calls
        * @return current counter of custom log calls
        */
        inline unsigned  int get_Ncustom() { return count_custom_log;};

        /**
         * Public interface to change log level
         * @param level  string, can be: dubug, info, warning,error.
         */
        void set_level(std::string level);

    protected:


        /** An enum storing the various log level.
         */
        enum _LOG_LEVEL {
            _notset = 0,       /**< lvl 0, Notset general value */
            _debug = 10,       /**< lvl 10, Debug level */
            _info = 20,        /**< lvl 20, Info level */
            _warning = 30,     /**< lvl 30, Warning level*/
            _error = 40,       /**< lvl 40, Error level */
            _critical = 100,   /**< lvl 100, Only critical level */
        };

        //const unsigned int MAX_N_WARNING=10;

        /**
         * Set the static  log_level.
         * @param level
         */
        inline void set_level(int level) {log_level=level;};



        static int log_level;    /*!< Current log level  */

        //GI This counter should be thread safe because each update is proteceted by the openmp atomic directive.
        static unsigned int count_debug;  /*!< Counter storing how many  times a debug log has been called*/
        static unsigned int count_info;   /*!< Counter storing how many  times a info log has been called*/
        static unsigned int count_warning; /*!< Counter storing how many  times a warning log has been called*/
        static unsigned int count_error;     /*!< Counter storing how many  times an error log has been called*/
        static unsigned int count_custom_log;   /*!< Counter storing how many  times a custom log has been called*/



        //NB the critical has not a counter since it always throws an exception, so we cannot have more than one call at runtime.
    };


}



#endif //SEVN_SEVNLOG_H