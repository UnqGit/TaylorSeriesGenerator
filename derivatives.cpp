#include "holder.h"

std::vector<long double> exp_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1, 0);
    der[0] = (std::exp(polycoef(x)));
    if(polycoef.isConstant()) return der;
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

std::vector<long double> pow_const_base_derivatives(const POLY& polycoef, long double const_base, long double x, int d){
    std::vector<long double> der(d+1, 0);
    der[0] = std::pow(const_base, polycoef(x));
    if(polycoef.isConstant()) return der;
    POLY dr = polycoef.derivative();
    int minP = std::min(d, polycoef.degree());
    std::vector<long double> derivs(minP, 0);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    long double helper = std::log(const_base);
    for(int i = 1; i < d + 1; i++){
        int min = std::min(i, minP);
        std::vector<long double> NCR = nCr[i-1];
        for(int j = 0; j < min; j++){
            der[i] += NCR[j]*der[i-1-j]*derivs[j];
        }
        der[i]*=helper;
    }
    return der;
}

std::vector<long double> sqrt_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1, 0);
    der[0] = std::sqrt(polycoef(x));
    if(polycoef.isConstant()) return der;
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

// std::vector<long double> sin_derivatives(const POLY& polycoef, long double x, int d){
//     std::vector<long double> der(d+1, 0);
//     long double sinHELP = std::sin(polycoef(x)), cosHELP = std::cos(polycoef(x));
//     der[0] = (sinHELP);
//     if(polycoef.isConstant()) return der;
//     POLY dr = polycoef.derivative();
//     if(d > 0) der[1] = cosHELP*(dr(x));
//     POLY old_A = dr.derivative();
//     POLY old_B = dr*dr;
//     POLY A = old_A;
//     POLY B = old_B;
//     for(int i = 2; i < d + 1; i++){
//         der[i] = cosHELP*A(x) - sinHELP*B(x);
//         old_A = A;
//         old_B = B;
//         A = old_A.derivative() - old_B*dr;
//         B = old_A*dr + old_B.derivative();
//     }
//     return der;
// }

// std::vector<long double> cos_derivatives(const POLY& polycoef, long double x, int d){
//     std::vector<long double> der(d+1, 0);
//     long double sinHELP = std::sin(polycoef(x)), cosHELP = std::cos(polycoef(x));
//     der[0] = (cosHELP);
//     if(polycoef.isConstant()) return der;
//     POLY dr = polycoef.derivative();
//     if(d > 0) der[1] = -sinHELP*(dr(x));
//     POLY old_A = dr.derivative();
//     POLY old_B = dr*dr;
//     POLY A = old_A;
//     POLY B = old_B;
//     for(int i = 2; i < d + 1; i++){
//         der[i] = -sinHELP*A(x)-cosHELP*B(x);
//         old_A = A;
//         old_B = B;
//         A = old_A.derivative() - old_B*dr;
//         B = old_A*dr + old_B.derivative();
//     }
//     return der;
// }

std::vector<long double> sin_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> result(d + 1, 0);
    int minP = std::min(d, polycoef.degree()) + 1;
    POLY dr = polycoef;
    std::vector<long double> derivs(minP, 0);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    long double helper = std::cos(derivs[0]);
    for(int i = 0; i < d + 1; i++){
        
    }
}

std::vector<long double> tan_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1, 0);
    der[0] = (std::tan(polycoef(x)));
    long double helper = der[0]*der[0] + 1.0;
    if(polycoef.isConstant()) return der;
    POLY dr = polycoef.derivative();
    int minP = d < polycoef.degree() ? d : polycoef.degree();
    std::vector<long double> derivs(minP);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }    
    for(int i = 1; i < d + 1; i++){
        int min = std::min(i, minP + 1);
        std::vector<long double> factor = nCr[i - 1];
        for(int j = 0; j < min - 1; j++){
            std::vector<long double> NCR = nCr[i - 2 - j];
            long double result = 0.0;
            for(int k = 0; k < NCR.size(); k++){
                result += NCR[k] * der[k] * der[NCR.size() - k];
            }
            der[i] += result * 2.0 * factor[j] * derivs[j];
        }
        if (i <= derivs.size()) der[i] += helper*derivs[i - 1];
    }
    return der;
}

std::vector<long double> tanh_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1, 0);
    der[0] = (std::tanh(polycoef(x)));
    long double helper = 1.0 - der[0]*der[0];
    if(polycoef.isConstant()) return der;
    POLY dr = polycoef.derivative();
    int minP = d < polycoef.degree() ? d : polycoef.degree();
    std::vector<long double> derivs(minP);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    for(int i = 1; i < d + 1; i++){
        int min = std::min(i, minP + 1);
        std::vector<long double> factor = nCr[i - 1];
        for(int j = 0; j < min - 1; j++){
            std::vector<long double> NCR = nCr[i - 2 - j];
            long double result = 0.0;
            for(int k = 0; k < NCR.size(); k++){
                result += NCR[k] * der[k] * der[NCR.size() - k];
            }
            der[i] -= result * 2.0 * factor[j] * derivs[j];
        }
        if(i<=derivs.size()) der[i] += helper*derivs[i - 1];
    }
    return der;
}

std::vector<long double> ln_const_base_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d+1, 0);
    long double help = polycoef(x);
    der[0] = std::log(help);
    if(polycoef.isConstant()) return der;
    POLY dr = polycoef.derivative();
    int minP = std::min(d, polycoef.degree());
    std::vector<long double> derivs(minP, 0);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    if(d > 0) der[1] = derivs[0]/help;
    for(int i = 2; i < d + 1; i++){
        std::vector<long double> NCR = nCr[i-1];
        int min = std::min(i , minP + 1) - 1;
        for(int j = 0; j < min; j++){
            der[i] -= NCR[j+1]*der[i - 1 - j]*derivs[min - 1 - j];
        }
        if(i <= minP) der[i] += derivs[i - 1];
        der[i]/=help;
    }
    return der;
}

std::vector<long double> atan_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d + 1, 0);
    long double pValue = polycoef(x);
    der[0] = std::atan(pValue);
    if(polycoef.isConstant()) return der;
    int minP = std::min(polycoef.degree(), d) + 1;
    std::vector<long double> derivs(d + 1, 0);
    derivs[0] = pValue;
    POLY dr = polycoef.derivative();
    for(int i = 1; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    long double pFactor = 1.0 + pValue * pValue;
    for(int i = 1; i < d + 1; i++){
        long double result = 0;
        std::vector<long double> factor = nCr[i - 1];
        for(int j = 1; j < i; j++){
            std::vector<long double> NCR = nCr[j - 1];
            long double res = 0;
            for(int k = 1; k < j + 1; k++){
                res += NCR[k - 1] * derivs[k - 1] * derivs[j + 1 - k];
            }
            result += factor[j] * der[i - j] * res;
        }
        der[i] = derivs[i] - 2.0 * result;
        der[i] /= pFactor;
    }
    return der;
}

std::vector<long double> pow_const_pow_derivatives(const POLY& polycoef, long double const_pow, long double x, int d){
    std::vector<long double> der(d+1, 0);
    long double helper = polycoef(x);
    der[0] = std::pow(helper, const_pow);
    if(polycoef.isConstant()) return der;
    int minP = std::min(polycoef.degree(), d);
    POLY dr = polycoef.derivative();
    std::vector<long double> derivs(minP, 0);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    for(int i = 1; i < d + 1; i++){
        std::vector<long double> NCR = nCr[i - 1];
        long double result = 0.0;
        int min = std::min(i, minP);
        for(int j = 1; j < min; j++){
            result += derivs[j - 1]*der[i - j]*(NCR[j - 1]*const_pow - NCR[j]);
        }
        if(i <= minP) result += const_pow*derivs[i - 1]*der[0];
        der[i] = result/helper;
    }
    return der;
}

std::vector<long double> division(std::vector<long double> numerator, std::vector<long double> denominator){
    std::vector<long double> result(numerator.size(), 0);
    long double helper = denominator[0];
    for(int i = 0; i < result.size(); i++){
        result[i] = numerator[i];
        std::vector<long double> NCR = nCr[i];
        for(int j = 0; j < i; j++){
            result[i] -= NCR[j]*denominator[i - j]*result[j];
        }
        result[i] /= helper;
    }
    return result;
}

std::vector<long double> asin_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d + 1, 0);
    int minP = std::min(polycoef.degree(), d) + 1;
    POLY dr = polycoef;
    std::vector<long double> derivs(minP, 0);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    der[0] = std::asin(derivs[0]);
    long double helper = std::sqrt(1 - derivs[0]*derivs[0]);
    if(polycoef.isConstant()) return der;
    for(int n = 1; n < d + 1; n++){
        std::vector<long double> factor = nCr[n - 1];
        long double outer_result = 0.0;
        for(int i = 1; i < n; i++){
            std::vector<long double> NCR = nCr[n-i-1];
            long double result = 0.0;
            int min = std::min(n - i, minP);
            for(int j = 0; j < min; j++){
                result += NCR[j]*derivs[j]*der[n-i-j];
            }
            outer_result += factor[i-1]*der[i]*result;
        }
        if(n <= minP) outer_result += derivs[n];
        der[n] = outer_result/helper;
    }
    return result;
}

std::vector<long double> asin_derivatives(const POLY& polycoef, long double x, int d){
    std::vector<long double> der(d + 1, 0);
    int minP = std::min(polycoef.degree(), d) + 1;
    POLY dr = polycoef;
    std::vector<long double> derivs(minP, 0);
    for(int i = 0; i < minP; i++){
        derivs[i] = dr(x);
        dr = dr.derivative();
    }
    der[0] = std::acos(derivs[0]);
    long double helper = std::sqrt(1 - derivs[0]*derivs[0]);
    if(polycoef.isConstant()) return der;
    for(int n = 1; n < d + 1; n++){
        std::vector<long double> factor = nCr[n - 1];
        long double outer_result = 0.0;
        for(int i = 1; i < n; i++){
            std::vector<long double> NCR = nCr[n-i-1];
            long double result = 0.0;
            int min = std::min(n - i, minP);
            for(int j = 0; j < min; j++){
                result += NCR[j]*derivs[j]*der[n-i-j];
            }
            outer_result += factor[i-1]*der[i]*result;
        }
        if(n <= minP) outer_result += derivs[n];
        der[n] = -outer_result/helper;
    }
    return result;
}

void derivatives(const std::string& func, int degree, POLY& p1, POLY& p2, long double point){
    std::vector<long double> derivatives;
    degree = degree > 100 ? 100 : degree;
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
    else if(func == "asin"){
        derivatives = sin_derivatives(p1, point, degree);
    }
    else if(func == "acos"){
        derivatives = cos_derivatives(p1, point, degree);
    }
    else if(func == "tan"){
        derivatives = tan_derivatives(p1, point, degree);
    }
    else if(func == "tanh"){
        derivatives = tanh_derivatives(p1, point, degree);
    }
    else if(func == "atan"){
        derivatives = atan_derivatives(p1, point, degree);
    }
    else if(func == "log"){
        if(p2.isConstant()){
            derivatives = ln_const_base_derivatives(p1, point, degree);
            if(!p2.isZero()){
                long double val = std::log(p2(point));
                std::transform(derivatives.begin(), derivatives.end(), derivatives.begin(), [a = val](auto p){return p/a;});
            }
        }
        else{
            derivatives = division(ln_const_base_derivatives(p1, point, degree), ln_const_base_derivatives(p2, point, degree));
        }
    }
    else if(func == "pow"){
        if(p2.isConstant()){
            long double val = p2(point);
            if(val >= 0 && val==std::floor(val)){
                int VAL = static_cast<int>(val);
                POLY result = p1^VAL;
                std::cout << result;
                std::cout << std::endl;
                return;
            }
            else derivatives = pow_const_pow_derivatives(p1, val, point, degree);
        }
        else if(p1.isConstant()){
            derivatives = pow_const_base_derivatives(p2, p1(point), point, degree);
        }
    }
    for(const auto& c : derivatives) std::cout << c << " ";
    std::cout << std::endl;
}
