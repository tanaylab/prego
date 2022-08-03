#ifndef base_stktrace_h
#define base_stktrace_h 1

class ActionTrace {
public:
        enum FuncWhere { FuncStart, FuncEnd };
public:
         ActionTrace();
        virtual  ~ActionTrace();
        virtual void  operator()(FuncWhere at_start) = 0;
};

#define MAX_SIZE 1000
class TraceMemWrite : public ActionTrace {
private:
        void *addr_;
        int size_;
        char mem_[MAX_SIZE];
public:
         TraceMemWrite(void *addr, int size);
        virtual  ~TraceMemWrite();
        virtual void  operator()(FuncWhere at_start);
};

#endif //base_stktrace_h
