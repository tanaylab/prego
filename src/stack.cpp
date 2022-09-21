#define ENV_FILE
#include "port.h"
#include "stktrace.h"
#include <stdlib.h>

#if DBG_ON

void call_tracers(ActionTrace::FuncWhere at_start);

struct CallTrace {
        const char *file;
        const char *func;
        int line;
};

static const int max_stack_depth = 128;
static CallTrace stack[max_stack_depth];
static int stack_top = -1;
#define max_tracers 100
int last_tracer = -1;
ActionTrace *tracers[max_tracers];

void stack_push(const char *file, int line, const char *func) {
        ENV_FUNC("stack_push")
        ASSERT(stack_top < max_stack_depth - 1,
                "function call stack overflow");
        stack_top++;
        stack[stack_top].file = file;
        stack[stack_top].line = line;
        stack[stack_top].func = func;
        Stk << file << '(' << line << "): " << func << " {" << endl; // }
#if HAS_HEAP_FUNCS
        memory_mode.deact();
        memory_mode.deact();
        if(memory_mode.to_act() && !mem_verify()) {
                cerr << "Detected at entry to:" << endl;
                stack_print(error_mode.out());
        }
        memory_mode.act();
        memory_mode.act();
#endif // HAS_HEAP_FUNCS
        call_tracers(ActionTrace::FuncStart);
}
void stack_pop() {
        ENV_FUNC("stack_pop")
        call_tracers(ActionTrace::FuncEnd);
#if HAS_HEAP_FUNCS
        memory_mode.deact();
        memory_mode.deact();
        if(memory_mode.to_act() && !mem_verify()) {
                cerr << "Detected at entry to:" << endl;
                stack_print(error_mode.out());
        }
        memory_mode.act();
        memory_mode.act();
#endif // HAS_HEAP_FUNCS
        Stk << stack[stack_top].file << '('
                << stack[stack_top].line << "): " // {
                << stack[stack_top].func << " }" << endl;
        ASSERT(stack_top-- >= 0, "function call stack underflow");
}
void stack_print(ostream &out) {
        out << "Function call stack: {" << endl; // }
        for(int i = 0; i <= stack_top; i++) {
                out << stack[i].file << '('
                    << stack[i].line << "): "
                    << stack[i].func << endl;
        } // {
        out << '}' << endl;
}
ActionTrace::ActionTrace() {
        ENV_FUNC("ActionTrace::ActionTrace")
        ASSERT(last_tracer + 1 < max_tracers,
                "Request for " << (max_tracers + 1) << " tracers");
        tracers[++last_tracer] = this;
}
ActionTrace::~ActionTrace() {
        ENV_FUNC("ActionTrace::~StackTrace")
        DBG_ASSERT(last_tracer > -1, "Deleting a non exsistent tracer");
        int tracer = 0;
        for(; tracers[tracer] != this && tracer <= last_tracer; tracer++)
        DBG_ASSERT(tracer <= last_tracer,
                "Deleted tracer wasn't found in the tracers array");
                
        for(; tracer < last_tracer; tracer++) 
                tracers[tracer] = tracers[tracer + 1];
        tracers[last_tracer] = 0;
        last_tracer--;
}
void call_tracers(ActionTrace::FuncWhere at_start) {
        ENV_FUNC("call_tracers")
        for (int t = 0; t <= last_tracer; t++)
                (*tracers[t])(at_start);
}
#else // DBG_ON 
ActionTrace::ActionTrace() {
        ENV_FUNC("ActionTrace::ActionTrace")
}
ActionTrace::~ActionTrace() {
        ENV_FUNC("ActionTrace::~StackTrace")
}
#endif // DBG_ON

const char *const sea_msg_text = "this function should not be called";
const char *const sea_msg() { return(sea_msg_text); }
