struct TERM{
    long double coef;
    int pow;
};

class DENSE_POLY;

struct complex {
    long double real, imag;
    complex(long double r = 0, long double i = 0) : real(r), imag(i) {}
    complex operator+(const complex& other) const { return {real + other.real, imag + other.imag}; }
    complex operator-(const complex& other) const { return {real - other.real, imag - other.imag}; }
    complex operator*(const complex& other) const {
        return {real * other.real - imag * other.imag, real * other.imag + imag * other.real};
    }
    complex& operator*=(const complex& other) {
        long double r = real * other.real - imag * other.imag;
        long double i = real * other.imag + imag * other.real;
        real = r; imag = i; return *this;
    }
    complex& operator/=(long double val) { real /= val; imag /= val; return *this; }
};

long double SQaM(long double a, int b){
    if(b == 0)return 1.0;
    long double result = 1.0;
    while(b>0){
        if(b%2)result*=a;
        a*=a;
        b>>=1;
    }
    return result;
}

long double horner_poly(std::vector<long double> v, int x, int min){
    long double result = 0.0;
    for(size_t i = v.size() - 1; i >= 0; i--){
        result = result*x+v[i];
    }
    result*=SQaM(x, min);
    return result;
}

long double poly(std::vector<TERM> v1, long double x){
    if(x == 0) return polycoef.back().pow == 0 ? polycoef.back().coef : 0;
    if(x == 1) return std::accumulate(polycoef.begin(), polycoef.end(), 0.0L, [](long double total, const auto& p){return total+p.coef;});
    if(x == -1) return std::accumulate(polycoef.begin(), polycoef.end(), 0.0L, [](long double total, const auto& p){return (p.pow%2)?total-p.coef:total+p.coef;});
    long double result = 0.0;
    if(polycoef.front().pow==polycoef.size()-1){
        for(size_t i = 0; i < polycoef.size(); i++)
            result=result*x+polycoef[i].coef;
    }
    else{
        for(size_t i = 0; i < polycoef.size(); i++){
            result += polycoef[i].coef;
            int j = polycoef[i].pow;
            if(i + 1 < polycoef.size()) j -= polycoef[i+1].pow;
            result*=SQaM(x, j);
        }
    }
    return result;
}

std::vector<TERM> poly_sum(std::vector<TERM> v1, std::vector<TERM> v2){
    std::vector<TERM> result;
    auto poly_sum_direct = [&](){
        int max = v1.front().pow > v2.front().pow ? v1.front().pow : v2.front().pow;
        std::vector<long double> temp(max+1, 0);
        for(const auto& [coefs, pows] : v1)temp[pows]+=coefs;
        for(const auto& [coefs, pows] : v2)temp[pows]+=coefs;
        std::vector<TERM> result, default_vec{{0,0}};
        for(int i = temp.size()-1; i >= 0; i--) if(temp[i]!=0) result.push_back({temp[i], i});
        return result.empty()?default_vec:result;
    }
    auto poly_sum_unomap = [&](){
        std::unordered_map<int, long double> res;
        for(const auto& [coefs, pows] : v1) res[pows]+=coefs;
        for(const auto& [coefs, pows] : v2) res[pows]+=coefs;
        std::vector<TERM> result, default_vec{{0,0}};
        for(const auto& [pows , coefs] : res)if(coef!=0)result.push_back({coefs, pows});
        std::sort(result.begin(), result.end(), [](const auto& a, const auto& b){return a.pow>b.pow;});
        return result.empty()?default_vec:result;
    }
    return (std::max(v1.front().pow, v2.front().pow)-v1.size()-v2.size() < 5.5*(v1.size()+v2.size())*(32 - __builtin_clz(static_cast<unsigned int>(v1.size()+v2.size())))) ? poly_sum_direct(v1, v2) : poly_sum_unomap(v1, v2);
}

std::vector<TERM> poly_diff(std::vector<TERM> v1, std::vector<TERM> v2){
    std::vector<TERM> result;
    auto poly_sum_direct = [&](){
        int max = v1.front().pow > v2.front().pow ? v1.front().pow : v2.front().pow;
        std::vector<long double> temp(max+1, 0);
        for(const auto& [coefs, pows] : v1)temp[pows]=coefs;
        for(const auto& [coefs, pows] : v2)temp[pows]-=coefs;
        std::vector<TERM> result, default_vec{{0,0}};
        for(int i = temp.size()-1; i >= 0; i--) if(temp[i]!=0) result.push_back({temp[i], i});
        return result.empty()?default_vec:result;
    }
    auto poly_sum_unomap = [&](){
        std::unordered_map<int, long double> res;
        for(const auto& [coefs, pows] : v1) res[pows]=coefs;
        for(const auto& [coefs, pows] : v2) res[pows]-=coefs;
        std::vector<TERM> result, default_vec{{0,0}};
        for(const auto& [pows , coefs] : res)if(coef!=0)result.push_back({coefs, pows});
        std::sort(result.begin(), result.end(), [](const auto& a, const auto& b){return a.pow>b.pow;});
        return result.empty()?default_vec:result;
    }
    return (std::max(v1.front().pow, v2.front().pow)-v1.size()-v2.size() < 5.5*(v1.size()+v2.size())*(32 - __builtin_clz(static_cast<unsigned int>(v1.size()+v2.size())))) ? poly_sum_direct(v1, v2) : poly_sum_unomap(v1, v2);
}

std::vector<TERM> unomap_poly_mult(const std::vector<TERM>& p1, const std::vector<TERM>& p2){
    std::unordered_map<int, long double> res;
    for(size_t i = 0; i < p1.size(); i++)
        for(size_t j = 0; j < p2.size(); j++)
            res[p1[i].pow+p2[j].pow]+=p1[i].coef*p2[j].coef;
    std::vector<TERM> result, default_vec;
    for(auto& [k, v] : res)
        if(v!=0)result.push_back({v, k});
    std::sort(result.begin(), result.end(), [](const auto& p1, const auto& p2){return p1.pow>p2.pow;});
    return result.empty()?default_vec:result;
}

void fft(std::vector<complex>& a, bool invert){
    int n = a.size();
    for(int i = 1, j = 0; i < n; ++i){
        int bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    static constexpr long double TAU = 2*M_PIl;
    for(int len = 2; len <= n; len <<= 1){
        long double ang = TAU / len * (invert ? -1 : 1);
        complex wlen(std::cos(ang), std::sin(ang));
        for(int i = 0; i < n; i += len){
            complex w(1);
            for(int j = 0; j < len / 2; j++){
                complex u = a[i + j];
                complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
    if(invert){
        for (auto& x : a) x /= n;
    }
}

std::vector<long double> fft_multiply(const std::vector<long double>& a, const std::vector<long double>& b){
    std::vector<complex> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 2LL << (31 - __builtin_clz(static_cast<unsigned int>(a.size() + b.size() - 1)));
    fa.resize(n);
    fb.resize(n);
    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++) fa[i] *= fb[i];
    fft(fa, true);
    std::vector<long double> result(n);
    for (int i = 0; i < n; i++) result[i] = fa[i].real;
    return result;
}

std::vector<long double> to_dense(const std::vector<TERM>& poly) {
    int max_deg = poly.front().second;
    int offset = poly.back().second;
    std::vector<long double> dense(max_deg - offset + 1, 0);
    for (const auto& [coef, pow] : poly)
        dense[pow - offset] = coef;
    return dense;
}

std::vector<TERM> to_sparse(const std::vector<long double>& poly, int offset) {
    std::vector<std::pair<long double, int>> result;
    for(int i = poly.size()-1; i >= 0; i--) {
        if(std::abs(poly[i])>3e-14) result.push_back({poly[i], i+offset});
    }
    return result;
}

std::vector<TERM> poly_fft_mult(const std::vector<TERM>& v1, const std::vector<TERM>& v2) {
    std::vector<long double> d1 = to_dense(v1);
    std::vector<long double> d2 = to_dense(v2);
    std::vector<long double> result_dense = fft_multiply(d1, d2);
    return to_sparse(result_dense, v1.back().pow+v2.back().pow);
}

std::vector<TERM> poly_mult(std::vector<TERM>& v1, std::vector<TERM>& v2){
    int n = (32 - __builtin_clz(static_cast<unsigned int>(v1.front().pow + v2.front().pow - 1)));
    int k = 1LL<<n;
    int m = v1.size();
    int l = v2.size();
    return (4*m*l+(m+l)*(31-__builtin_clz(static_cast<unsigned int>(m+l))) > 3*(n*(5+4*k)+3*n) - v1.back().pow - v2.back().pow)?poly_fft_mult(v1, v2):unomap_poly_mult(v1, v2);
}

std::vector<TERM> poly_raise(std::vector<TERM> base_poly, int power){
    if(power==0){
        std::vector<TERM> res{{1,0}};
        return res;
    }
    if(power==1)return base_poly;
    if(base_poly.size()==1){
        std::vector<TERM> res{{SQaM(base_poly[0].coef), base_poly[0].pow*power}};
        return res;
    }
    if(base_poly.size()==2){
        std::vector<long double> NCR = nCr(power);
        long double max_coef = SQaM(base_poly[0].coef);
        int max_pow = power*base_poly[0].pow;
        long double diff = base_poly[0].coef/base_poly[1].coef;
        int pow_diff = base_poly[1].pow-base_poly[0].pow;
        std::vector<TERM> res(power+1);
        for(int i = 0; i < power + 1; i++){
            res[i]={NCR[i]*max_coef, max_pow};
            max_coef*=diff;
            max_pow+=pow_diff;
        }
        return res;
    }
    std::vector<TERM> temp = base_poly, res = {{1,0}};
    while(power>0){
        if(power&1)res=poly_mult(res, temp);
        temp=poly_mult(temp, temp);
        power>>=1;
    }
    return res;
}

std::vector<TERM> poly_compose(std::vector<TERM> v1, std::vector<TERM> v2){
    std::unordered_map<int, long double> res;
    std::vector<TERM> med{{1, 0}};
    long double outer_coef;
    int exp, prev_exp = 0;
    if(v1.back().pow==0) res[0]=v1.back().coef;
    else{
        outer_coef = v1.back().coef;
        prev_exp = v1.back().pow;
        med = poly_raise(v2, prev_exp);
        for(const auto& [coefs, pows] : med)
            res[pows]+=outer_coef*coefs;
    }
    for(int i = v1.size()-2; i >= 0; i--){
        outer_coef = v1[i].coef;
        exp = v1[i].pow;
        std::vector<TERM> temp = poly_raise(v2, exp-prev_exp);
        med = poly_mult(med, temp);
        for(const auto& [coefs, pows] : med)
            res[pows]+=outer_coef*coefs;
        prev_exp = exp;
    }
    std::vector<TERM> result, default_vec{{0,0}};
    for(auto& [pows, coefs] : res)
        if(coefs!=0)result.push_back({coefs, pows});
    std::sort(result.begin(), result.end(), [](const auto& a, const auto& b){ return a.pow > b.pow;});
    return result.empty()?default_vec:result;
}

std::vector<TERM> poly_ind_int(std::vector<TERM> v){
    for(auto& [coefs, pows] : v){
        coefs/=(pows+1);
        pow++;
    }
    return v;
}

long double poly_def_int(std::vector<TERM> v, long double a, long double b){
    std::vector<TERM> V = poly_ind_int(v);
    return poly(V, b)-poly(V, a);
}

std::string power_display(int i){
    std::string result = "";
    if(i == 0) return result;
    if(i == 1) return "x";
    std::string super[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
    while(i > 0){
        result = super[i%10] + result;
        i/=10;
    }
    result = "x"+result;
    return result;
}

class POLY;

class DENSE_POLY{
    public:
    
        std::vector<long double> coefficients;
        
        int min_deg = 0;
        
        int degree(){
            return this->coefficients.size()+min_deg-1;
        }
        
        POLY wet(){
            POLY result;
            std::vector<TERM>& alias = result.terms;
            for(size_t i = this->coefficients.size()-1; i >= 0; i--)
                if(coefficients[i]!=0) alias.push_back({coefficients[i], i+min_deg});
            std::vector<TERM> default_vec{{0,0}};
            if(alias.empty()) alias = default_vec;
            return result;
        }
        
        long double operator() (const long double point){
            return horner_poly(this->coefficients, point, min_deg);
        }
}

class POLY{
    public:
        
        std::vector<TERM> terms;
        
        long double operator() (long double x){
            return poly(this->terms, x);
        }
        
        void negate(){
            std::vector<TERM>& alias = this->terms;
            std::transform(alias.begin(), alias.end(), alias.begin(), [](auto p){return p.coef=-p.coef;});
            return *this;
        }
        
        POLY operator-() const{
            POLY result;
            result.terms = this->terms;
            std::transform(result.terms.begin(), result.terms.end(), result.terms.begin(), [](auto p){return p.coef=-p.coef;});
            return result;
        }
        
        POLY operator+ (const POLY& v2){
            POLY result;
            result.terms = poly_sum(this->terms, v2.terms);
            return result;
        }
        
        POLY operator- (const POLY& v2){
            POLY result;
            result.terms = poly_diff(this->terms, v2.terms);
            return result;
        }
        
        POLY operator* (const POLY& v2){
            POLY result;
            result.terms = poly_mult(this->terms, v2.terms);
            return result;
        }
        
        POLY operator/ (const POLY& v2){
            //....not yet made the function
        }
        
        POLY& operator+= (POLY v){
            this->terms = poly_sum(this->terms, v.terms);
        }
        
        POLY& operator-= (const POLY& v2){
            this->terms = poly_diff(this->terms, v2.terms);
        }
        
        POLY& operator*= (POLY v){
            this->terms = poly_mult(this->terms, v.terms);
        }
        
        POLY& operator/= (POLY v){
            //...not yet made the logic
        }
        
        POLY operator+ (long double val){
            POLY result;
            result.terms = this->terms;
            if(result.terms.back().pow!=0)result.terms.push_back({val, 0});
            else result.terms.back().coef += val;
            return result;
        }
        
        POLY operator- (long double val){
            POLY result;
            result.terms = this->terms;
            if(result.terms.back().pow!=0)result.terms.push_back({val, 0});
            else result.terms.back().coef -= val;
            return result;
        }
        
        POLY operator* (long double val){
            POLY result;
            result.terms = this->terms;
            std::transform(result.terms.begin(), result.terms.end(), result.terms.begin(), [a = val](auto p){return p.coef*a;});
            return result;
        }
        
        POLY operator/ (long double val){
            POLY result;
            result.terms = this->terms;
            std::transform(result.terms.begin(), result.terms.end(), result.terms.begin(), [a = val](auto p){return p.coef/a;});
            return result;
        }
        
        POLY& operator+= (long double val){
            std::vector<TERM> alias& = this->terms;
            if(alias.back().pow!=0) alias.push_back({val, 0});
            else alias.back().coef += val;
            return *this;
        }
        
        POLY& operator-= (long double val){
            std::vector<TERM> alias& = this->terms;
            if(alias.back().pow!=0) alias.push_back({val, 0});
            else alias.back().coef -= val;
            return *this;
        }
        
        POLY& operator*= (long double val){
            std::vector<TERM> alias& = this->terms;
            std::transform(alias.begin(), alias.end(), alias.begin(), [a = val](auto p){return p.coef*a;});
            return *this;
        }
        
        POLY& operator/= (long double val){
            std::vector<TERM> alias& = this->terms;
            std::transform(alias.begin(), alias.end(), alias.begin(), [a = val](auto p){return p.coef/a;});
            return *this;
        }
        
        POLY operator^ (int power){
            POLY result;
            result.terms = poly_raise(this->terms, power);
            return result;
        }
        
        POLY operator^= (int power){
            this->terms = poly_raise(this->terms, power);
            return *this;
        }
        
        DENSE_POLY dry(){
            std::vector<TERM>& alias = this->terms;
            std::vector<long double> temp(alias.front().pow-alias.back().pow+1, 0);
            for(auto& [coefs, pows] : alias){
                temp[pows] = coefs;
            }
            DENSE_POLY dense;
            dense.coefficients = temp;
            dense.min_deg = alias.back().pow;
            return dense;
        }
        
        long double operator[] (int power){
            if(power > this->terms.front().pow) return 0;
            if(power < 0) return 0;
            int idx = binary_search(this->terms, power);
            if(idx < 0) return 0;
            else return this->terms[idx].coef;
        }
        
        void print(std::ostream& os) const{
            std::vector<TERM>& alias = this->terms;
            if(alias.empty()){
                os << 0;
                return;
            }
            if(alias.front().pow==0) os << alias[0].coef;
            else{
                if(std::abs(alias[0].coef)==1) os << (alias[0].coef>0?"":"-");
                else os << (alias[0].coef>0?"":"-") << std::abs(alias[0].coef);
                os << power_display(alias[0].pow);
                for(size_t i = 1; i < alias.size(); i++){
                    os << (alias[i].coef>0?" + ":" - ");
                    if((alias[i].coef)==1&&alias[i].pow!=0){}
                    else os << std::abs(alias[i].coef);
os << power_display(alias[i].pow);
                    if(i % 10 == 0) os << "\n";
                }
            }
        }
        
        int degree(){
            return this->terms.front().pow;
        }
        
        POLY derivative(POLY v){
            POLY result;
            result.terms = poly_der(v.terms);
            return result;
        }
        
        POLY integrate(POLY v){
            POLY result;
            result.terms = poly_indef_int(v.terms);
            return result;
        }
        
        long double definite_integrate(long double a, long double b){
            return poly_def_int(this->terms, a, b);
        }
        
        POLY compose(POLY v){
            POLY result;
            result.terms = poly_compose(this->terms, v);
            return result;
        }
}

std::ostream& operator<<(std::ostream& os, const POLY& p) {
    p.print(os);
    return os;
}
