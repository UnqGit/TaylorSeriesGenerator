#include "holder.h"

std::vector<long double> exp_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1);
    der[0] = (std::exp(polycoef(x)));
    if(polycoef.isConstant()){
        if(d > 0) std::fill(der.begin()+1, der.end(), 0);
        return der;
    }
    POLY dr = polycoef.derivative();
    int minP = d < polycoef.degree() ? d : polycoef.degree();
    std::vector<long double> derivs(minP);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    for(int i = 1; i < d + 1; i++){
        int min = i < derivs.size() ? i : derivs.size();
        std::vector<long double> NCR = nCr[i-1];
        for(int j = 0; j < min; j++)
            der[i]+=NCR[j]*der[i-j-1]*derivs[j];
    }
    return der;
}

std::vector<long double> sinh_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> first = exp_derivatives(polycoef, x, d);
    std::vector<long double> second = exp_derivatives(-polycoef, x, d);
    std::vector<long double> final(first.size(), 0);
    for(int i = 0; i < final.size(); i++)
        final[i] = (first[i] - second[i])/2.0;
    return final;
}

std::vector<long double> cosh_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> first = exp_derivatives(polycoef, x, d);
    std::vector<long double> second = exp_derivatives(-polycoef, x, d);
    std::vector<long double> final(first.size(), 0);
    for(int i = 0; i < final.size(); i++)
        final[i] = (first[i] + second[i])/2.0;
    return final;
}

std::vector<long double> sqrt_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1);
    der[0] = std::sqrt(polycoef(x));
    if(polycoef.isConstant()){
        if(d>0)std::fill(der.begin()+1, der.end(), 0.0);
        return der;
    }
    POLY dr = polycoef.derivative();
    int minP = d < polycoef.degree() ? d : polycoef.degree();
    std::vector<long double> derivs(minP);
    for(int i = 1; i <= minP; i++){
        der[i] = dr(x)/2.0;
        dr = dr.derivative();
    }
    for(int i = 1; i < d + 1; i++){
        std::vector<long double> NCR = nCr[i-1];
        for(int j = 1; j < i; j++)
            der[i]-=NCR[j-1]*der[j]*der[i-j];
        der[i]/=der[0];
    }
    return der;
}

std::vector<long double> sin_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1);
    long double sinHELP = std::sin(polycoef(x)), cosHELP = std::cos(polycoef(x));
    der[0] = (sinHELP);
    if(polycoef.isConstant()){
        if(d > 0) std::fill(der.begin()+1, der.end(), 0);
        return der;
    }
    POLY dr = polycoef.derivative();
    if(d > 0) der[1] = cosHELP*(dr(x));
    POLY old_A = dr.derivative();
    POLY old_B = dr*dr;
    POLY A = old_A;
    POLY B = old_B;
    for(int i = 2; i < d + 1; i++){
        der[i] = cosHELP*A(x) - sinHELP*B(x);
        old_A = A;
        old_B = B;
        A = old_A.derivative() - old_B*dr;
        B = old_A*dr + old_B.derivative();
    }
    return der;
}

std::vector<long double> cos_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1);
    long double sinHELP = std::sin(polycoef(x)), cosHELP = std::cos(polycoef(x));
    der[0] = (cosHELP);
    if(polycoef.isConstant()){
        if(d > 0) std::fill(der.begin()+1, der.end(), 0);
        return der;
    }
    POLY dr = polycoef.derivative();
    if(d > 0) der[1] = -sinHELP*(dr(x));
    POLY old_A = dr.derivative();
    POLY old_B = dr*dr;
    POLY A = old_A;
    POLY B = old_B;
    for(int i = 2; i < d + 1; i++){
        der[i] = -sinHELP*A(x)-cosHELP*B(x);
        old_A = A;
        old_B = B;
        A = old_A.derivative() - old_B*dr;
        B = old_A*dr + old_B.derivative();
    }
    return der;
}

std::vector<long double> tan_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1, 0.0);
    der[0] = (std::tan(polycoef(x)));
    long double helper = der[0]*der[0] + 1.0;
    if(polycoef.isConstant()){
        if(d > 0) std::fill(der.begin()+1, der.end(), 0);
        return der;
    }
    POLY dr = polycoef.derivative();
    int minP = d < polycoef.degree() ? d : polycoef.degree();
    std::vector<long double> derivs(minP);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    std::vector<long double> factor(minP, 0.0);
    if(d > 0){
        der[1] = derivs[0]*helper;
        factor[0] = 2;
    }
    auto p = [&](){
        std::vector<long double> newfactor(factor.begin(), factor.end());
        for(int i = 1; i < factor.size(); i++){
            if(factor[i]==0){
                newfactor[i] = 2 + factor[i - 1];
                break;
            }
            else newfactor[i] += factor[i - 1];
        }
        factor = std::move(newfactor);
    };
    for(int i = 2; i < d + 1; i++){
        int min = i < derivs.size() + 1 ? i : derivs.size() + 1;
        for(int j = 0; j < min - 1; j++){
            std::vector<long double> NCR = nCr[i - 2 - j];
            long double result = 0.0;
            for(int k = 0; k < NCR.size(); k++){
                result += NCR[k] * der[k] * der[NCR.size() - k];
            }
            der[i] += result * factor[j] * derivs[j];
        }
        if (i <= derivs.size()) der[i] += helper*derivs[i - 1];
        p();
    }
    return der;
}

std::vector<long double> tanh_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1, 0.0);
    der[0] = (std::tanh(polycoef(x)));
    long double helper = 1.0 - der[0]*der[0];
    if(polycoef.isConstant()){
        if(d > 0) std::fill(der.begin()+1, der.end(), 0);
        return der;
    }
    POLY dr = polycoef.derivative();
    int minP = d < polycoef.degree() ? d : polycoef.degree();
    std::vector<long double> derivs(minP);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    std::vector<long double> factor(minP, 0.0);
    if(d > 0){
        der[1] = derivs[0]*helper;
        factor[0] = 2;
    }
    auto p = [&](){
        std::vector<long double> newfactor(factor.begin(), factor.end());
        for(int i = 1; i < factor.size(); i++){
            if(factor[i]==0){
                newfactor[i] = 2 + factor[i - 1];
                break;
            }
            else newfactor[i] += factor[i - 1];
        }
        factor = std::move(newfactor);
    };
    for(int i = 2; i < d + 1; i++){
        int min = i < derivs.size() + 1 ? i : derivs.size() + 1;
        for(int j = 0; j < min - 1; j++){
            std::vector<long double> NCR = nCr[i - 2 - j];
            long double result = 0.0;
            for(int k = 0; k < NCR.size(); k++){
                result += NCR[k] * der[k] * der[NCR.size() - k];
            }
            der[i] += result * factor[j] * derivs[j];
        }
        der[i] = ((i<=derivs.size())?helper*derivs[i - 1] : 0) - der[i];
        p();
    }
    return der;
}

void derivatives(const std::string& func, int degree, POLY& p1, POLY& p2, long double point){
    std::vector<long double> derivatives;
    if(func == "exp"){
        derivatives = exp_derivatives(p1, point, degree);
    }
    else if(func == "sinh"){
        derivatives = sinh_derivatives(p1, point, degree);
    }
    else if(func == "cosh"){
        derivatives = cosh_derivatives(p1, point, degree);
    }
    else if(func == "sqrt"){
        derivatives = sqrt_derivatives(p1, point, degree);
    }
    else if(func == "sin"){
        derivatives = sin_derivatives(p1, point, degree);
    }
    else if(func == "cos"){
        derivatives = cos_derivatives(p1, point, degree);
    }
    else if(func == "tan"){
        derivatives = tan_derivatives(p1, point, degree);
    }
    else if(func == "tanh"){
        derivatives = tanh_derivatives(p1, point, degree);
    }
}
