#ifndef base_init_h
#define base_init_h 1

#include "error.h"

#if DBG_ON
class TraceFileInit {
private:
        const char *file_;
        int line_;
public:
        TraceFileInit(const char *file, int line = -1):
                file_(file),
                line_(line)
        {
                if(!getenv("BASE_INIT_LOG")) {
                        file_ = 0;
                        return;
                }
                Rcpp::Rcerr << "Construct objects for " << file_;
                if(line_ > 0)
                        Rcpp::Rcerr << " at " << line_;
                Rcpp::Rcerr << endl;
        }
        ~TraceFileInit() {
                if(!file_)
                        return;
                Rcpp::Rcerr << "Destruct objects for " << file_;
                if(line_ > 0)
                        Rcpp::Rcerr << " at " << line_;
                Rcpp::Rcerr << endl;
        }
};
#endif // DBG_OM
extern void  at_finish(void (*func)());
extern const char * version();
extern bool  arg_version(int argc, char *argv[]);

#if DBG_ON
# define TRACE_INIT(ID) static TraceFileInit ID(__FILE__, __LINE__);
#else
# define TRACE_INIT(ID)
#endif // DBG_ON
#if TRACE
# if DBG_ON
#  define BASE_CC_FILE \
        static TraceFileInit trace_file_init_(__FILE__); \
        static ActiveMode &error_mode = get_error_mode(); 
//      static ActiveMode &memory_mode = get_memory_mode(); 
//	static ActiveMode &trace_mode = get_trace_mode(); 
//	static ActiveMode &stack_mode = get_stack_mode(); 
//	static ActiveMode &dbg_mode = get_dbg_mode();
# else // DBG_ON
#  define BASE_CC_FILE \
        static ActiveMode &error_mode = get_error_mode(); 
//        static ActiveMode &memory_mode = get_memory_mode(); 
//        static ActiveMode &trace_mode = get_trace_mode(); 
//        static InactiveMode &dbg_mode = get_dbg_mode();
# endif // DBG_ON
#else // TRACE
# if DBG_ON
#  define BASE_CC_FILE \
        static TraceFileInit trace_file_init_(__FILE__); \
        static ActiveMode &error_mode = get_error_mode();
//      static ActiveMode &memory_mode = get_memory_mode(); 
//      static InactiveMode &trace_mode = get_trace_mode();
//      static ActiveMode &stack_mode = get_stack_mode();
//	static ActiveMode &dbg_mode = get_dbg_mode();
# else // DBG_ON
#  define BASE_CC_FILE \
        static ActiveMode &error_mode = get_error_mode();
//        static ActiveMode &memory_mode = get_memory_mode();
//        static InactiveMode &trace_mode = get_trace_mode();
//        static InactiveMode &dbg_mode = get_dbg_mode();
# endif // DBG_ON
#endif // TRACE

#ifdef ENV_FILE
extern ActiveMode base_error_mode;
#define error_mode base_error_mode
//extern ActiveMode base_memory_mode;
//#define memory_mode base_memory_mode
#if DBG_ON
extern ActiveMode base_stack_mode;
#define stack_mode base_stack_mode
#endif // DBG_ON
#if TRACE
extern ActiveMode base_trace_mode;
#define trace_mode base_trace_mode
#endif // TRACE
#if DBG_ON
extern ActiveMode base_dbg_mode;
#define dbg_mode base_dbg_mode
#endif // DBG_ON
#endif // ENV_FILE

#endif // base_init_h
