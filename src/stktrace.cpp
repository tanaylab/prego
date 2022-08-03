#define ENV_FILE
#include "port.h"
#include "stktrace.h"

#if DBG_ON 
TraceMemWrite::TraceMemWrite(void *addr, int size) {
        ENV_FUNC("TraceMemWrite::TraceMemWrite")
        addr_ = addr;
        size_ = size;
        DBG_ASSERT(size < MAX_SIZE, "Maximal size is " << int(MAX_SIZE));
        memcpy(mem_, addr_, size);
}
void TraceMemWrite::operator()(FuncWhere at_start) {
        ENV_FUNC("TraceMemWrite::operator()")
        if(memcmp(mem_, addr_, size_) == 0)
                return;

        stack_print(error_mode.out());
        if(at_start == FuncStart)
                error_mode.out() << "Start of function: ";
        else
                error_mode.out() << "End of function: ";

        error_mode.out() << "Adress " << int(addr_) <<
                ", Range " << size_ << " changed. " << endl;
        memcpy(mem_, addr_, size_);
}
#else // DBG_ON
TraceMemWrite::TraceMemWrite(void *, int) {
        ENV_FUNC("TraceMemWrite::TraceMemWrite")
        addr_ = 0;
        size_ = 0;
}
void TraceMemWrite::operator()(FuncWhere ) {
        ENV_FUNC("TraceMemWrite::operator()")
}
#endif // DBG_ON
TraceMemWrite::~TraceMemWrite() {
        ENV_FUNC("TraceMemWrite::~TraceMemWrite")
}
