#pragma once

#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <map>
#include <unordered_map>

struct TERM{
    long double coef;
    int pow;
    TERM(long double c = 0, int p = 0) : coef(c), pow(p) {}
};

class POLY;

class DENSE_POLY{
    public:
    
        std::vector<long double> coefficients;
        
        int min_deg = 0;
        
        int degree() const;
        
        auto begin() const;
        
        auto end() const;
        
        auto begin();
        
        auto end();
        
        POLY dry() const;
        
        long double operator() (const long double point) const;
        
        DENSE_POLY negate();
        
        DENSE_POLY operator-() const;
        
        DENSE_POLY operator+ (const DENSE_POLY& other) const;
        
        DENSE_POLY operator- (const DENSE_POLY& other) const;
        
        DENSE_POLY& operator+= (const DENSE_POLY& other);
        
        DENSE_POLY& operator-= (const DENSE_POLY& other);
        
        DENSE_POLY clean();
        
        DENSE_POLY operator* (const DENSE_POLY& other) const;
        
        DENSE_POLY operator*= (const DENSE_POLY& other);
        
        DENSE_POLY operator+ (const long double val) const;
        
        DENSE_POLY operator- (const long double val) const;
        
        DENSE_POLY operator* (const long double val) const;
        
        DENSE_POLY operator/ (const long double val) const;
        
        DENSE_POLY& operator+= (const long double val);
        
        DENSE_POLY& operator-= (const long double val);
        
        DENSE_POLY& operator*= (const long double val);
        
        DENSE_POLY& operator/= (const long double val);
};

class POLY{
    public:
        
        std::vector<TERM> terms;
        
        long double operator() (long double x) const;
        
        POLY& negate();
        
        auto begin() const;
        
        auto end() const;
        
        auto begin();
        
        auto end();
        
        POLY operator-() const;
        
        POLY operator+ (const POLY& v2) const;
        
        POLY operator- (const POLY& v2) const;
        
        POLY operator* (const POLY& v2) const;
        
        // POLY operator/ (const POLY& v2) const;
        //     ....not yet made the function
        
        POLY& operator+= (POLY v);
        
        POLY& operator-= (const POLY& v2);
        
        POLY& operator*= (POLY v);
        
        // POLY& operator/= (POLY v);
        //     ...not yet made the logic
        
        POLY operator+ (long double val) const;
        
        POLY operator- (long double val) const;
        
        POLY operator* (long double val) const;
        
        POLY operator/ (long double val) const;
        
        POLY& operator+= (long double val);
        
        POLY& operator-= (long double val);
        
        POLY& operator*= (long double val);
        
        POLY& operator/= (long double val);
        
        POLY operator^ (int power) const;
        
        POLY operator^= (int power);
        
        DENSE_POLY wet() const;
        
        /*
        long double operator[] (int power) const;
        not doing it rn.
        */
        
        void print(std::ostream& os) const;
        
        int degree() const;
        
        POLY derivative(POLY v) const;
        
        POLY integrate(POLY v) const;
        
        long double definite_integrate(long double a, long double b) const;
        
        POLY compose(POLY v) const;
};

extern std::vector<std::vector<long double>> nCr;
