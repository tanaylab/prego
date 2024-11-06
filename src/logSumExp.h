inline void log_sum_log(double &l1, double l2) {
    if(l1 > l2) {
        if(!std::isinf(l2)) {
            l1 += log(1 + exp(l2-l1));
        }
    } else {
        if(std::isinf(l1)) {
            l1 = l2;
        } else {
            l1 = l2 + log(1 + exp(l1-l2));
        }
    }
}

inline double log_sum_exp(const double* x, int n) {
    if(n == 0) return -std::numeric_limits<double>::infinity();
    if(n == 1) return x[0];
    
    // Find maximum using indices for better cache usage
    double max_val = x[0];
    for(int i = 1; i < n; ++i) {
        if(x[i] > max_val) max_val = x[i];
    }

    if(std::isinf(max_val)) return max_val;

    // Compute sum of exp(x[i] - max_val)
    double sum = 0.0;
    constexpr int BLOCK_SIZE = 4;  // Process 4 elements at a time
    
    // Process blocks of 4
    int i = 0;
    for(; i <= n-BLOCK_SIZE; i += BLOCK_SIZE) {
        sum += exp(x[i] - max_val);
        sum += exp(x[i+1] - max_val);
        sum += exp(x[i+2] - max_val);
        sum += exp(x[i+3] - max_val);
    }
    
    // Handle remaining elements
    for(; i < n; ++i) {
        sum += exp(x[i] - max_val);
    }
    
    return max_val + log(sum);
}