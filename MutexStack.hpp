template<class data, typename my_complex>
struct MutexRectangles {
    std::stack<data> S;
    mutable std::mutex m;

    void push(data val) {
        std::lock_guard<std::mutex> lock(m);
        S.emplace(val);
    }

    data pop(my_complex NIL) {
        std::lock_guard<std::mutex> lock(m);
        if(S.empty()) {
            return data(NIL, NIL, 0, 0, 0, 0);
        }
        data res = S.top();
        S.pop();
        return res;
    }

    int size() {
        std::lock_guard<std::mutex> lock(m);
        return S.size();
    }
};

template<class answer, typename my_complex>
struct MutexAnswer {
    std::stack<answer> S;
    mutable std::mutex m;

    void push(answer val, int &cnt, int total_zeros) {
        std::lock_guard<std::mutex> lock(m);
        S.emplace(val);
        cnt++;
    }

    answer pop(my_complex NIL) {
        std::lock_guard<std::mutex> lock(m);
        if(S.empty()) {
            return answer(NIL, 0);
        }
        answer res = S.top();
        S.pop();
        return res;
    }

    int size() {
        std::lock_guard<std::mutex> lock(m);
        return S.size();
    }
};
