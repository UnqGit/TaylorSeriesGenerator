#include "holder.h"

void start_text(){
    std::cout << "In the program you can enter coefficients for in multiple ways:\n";
    std::cout << "\t1.) INPUT: 3 2 0 1; the function would be 3x³ + 2x² + 0x¹ + 1.\n";
    std::cout << "\t2.) INPUT: 3,2,0,1; function would be same as above(i.e, highest degree first).\n";
    std::cout << "\t3.) INPUT: 3x^3 + 2x^2 + 1; these terms can be in any order(power denoted via '^').\n";
    std::cout << "\t4.) INPUT: 3x3+2x2+1; this is also valid(number following x is regarded as power).\n";
    std::cout << "\t4.) INPUT: 3x³+2x²+1; this is also valid(superscipt as power).\n";
    std::cout << "\t5.) INPUT: 3x^3 2x^2 1; this is also valid.\n";
    std::cout << "\t6.) INPUT: 3X³ 2X^2 + x1 + 1; this is also valid.\n\n";
    std::cout << "NOTE: things like:\n";
    std::cout << "\t1.) 2² are not valid and would be considered as 22.\n";
    std::cout << "\t2.) Parenthesis(\"()\") are not allowed.\n";
    std::cout << "\t3.) x²x is invalid.\n";
    std::cout << "\t4.) x-+x² is also invalid.\n\n";
    std::cout << "The available functions will be\npow, sqrt, exp, log, tan, sin, cos, atan, acos, asin, sinh, cosh, tanh\n\n";
}

void input_corrector(std::string& input){
    input = std::regex_replace(input, std::regex("[,\t\n]+"), " ");
    input.erase(0, input.find_first_not_of(" 0"));
    if(input.empty()||input==" ") input = "0";
    else if(input.back()==' ')input.pop_back();
    for(char& ch : input) ch = std::tolower(ch);
}

bool invalid_input(std::string& input){
    if(std::regex_search(input, std::regex("[^-x+ 0123456789\\.\\^⁰¹²³⁴⁵⁶⁷⁸⁹]"))) return true;
    std::string super[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
    for(int i = 0; i <= 9; i++){
        input = std::regex_replace(input, std::regex(super[i]), std::to_string(i));
        std::string temp = "x"+std::to_string(i);
        std::string todo = "x^"+std::to_string(i);
        input = std::regex_replace(input, std::regex(temp), todo);
    }
    if(std::regex_search(input, std::regex("x\\^[- ]"))) return true;
    input = std::regex_replace(input, std::regex("x\\^\\+?"), "#");
    if(input.find('^')!=std::string::npos) return true;
    input = std::regex_replace(input, std::regex("x"), "#1");
    input = std::regex_replace(input, std::regex("([+-])"), " $1");
    input = std::regex_replace(input, std::regex(" #"), " 1#");
    input = std::regex_replace(input, std::regex("\\-#"), "-1#");
    input = std::regex_replace(input, std::regex("\\+#"), "+1#");
    input = std::regex_replace(input, std::regex("\\s+"), " ");
    input = std::regex_replace(input, std::regex("([+-]|#)\\s"), "$1");
    if(input[0]=='#') input.insert(input.begin(), '1');
    for(int i = 0; i < input.size(); i++){
        if(input[i]=='#'||input[i]=='.'){
          if(i+1==input.length()) return true;
          else if(!(std::isdigit(input[i+1]))) return true;
            while(i+1<input.length()){
                if(input[i+1]==' ')break;
                if(!std::isdigit(input[i+1])) return true;
                i++;
            }
        }
     }
    if(std::regex_search(input, std::regex("[+-][^0-9#\\.]"))) return true;
    if(input.back()=='-'||input.back()=='+') return true;
    return false;
}

bool cond(const std::string& part,const std::string& str){
    return str.find(part)!=std::string::npos;
}
bool function_checker(std::string& func){
    if(cond("ex", func)){func = "exp"; return 0;}
    if(cond("lo", func)){func = "log"; return 0;}
    if(cond("at", func)){func = "atan"; return 0;}
    if(cond("as", func)){func = "asin"; return 0;}
    if(cond("ac", func)){func = "acos"; return 0;}
    if(cond("sq", func)){func = "sqrt"; return 0;}
    if(cond("inh", func)){func = "sinh"; return 0;}
    if(cond("nh", func)){func = "tanh"; return 0;}
    if(cond("sh", func)){func = "cosh"; return 0;}
    if(cond("ta", func)){func = "tan"; return 0;}
    if(cond("si", func)){func = "sin"; return 0;}
    if(cond("co", func)){func = "cos"; return 0;}
    if(cond("po", func)){func = "pow"; return 0;}
    return 1;
}

std::string get_function(){
    std::cout << "Enter the function that you want: ";
    std::string func;
    do{
        std::getline(std::cin, func);
        input_corrector(func);
        if(!function_checker(func)) break;
        std::cout << "Incorrect function name, enter again: ";
        std::cin.clear();
    }while(true);
    std::cin.clear();
    return func;
}

std::string get_poly(){
    std::string input;
    do{
        std::getline(std::cin, input);
        input_corrector(input);
        if(!invalid_input(input)) break;
        std::cout << "Wrong entry, please enter in the correct format as stated at the starting:\n";
        std::cin.clear();
    }while(true);
    std::cin.clear();
    return input;
}

POLY get_vec(const std::string& input){
    std::vector<TERM> poly_terms;
    std::vector<std::string> separation;
    std::stringstream ss(input);
    std::string temp, temp2;
    int power;
    while(ss >> temp) separation.push_back(temp);
    if(input.find('#')!=std::string::npos){
        for(int i = 0; i < separation.size(); i++){
            temp = separation[i];
            if(temp.find('#')!=std::string::npos){
                int z = temp.find('#');
                TERM t = TERM(std::stold(temp.substr(0, z)), std::stoi(temp.substr(z + 1, temp.size() - z)));
                poly_terms.push_back(t);
            }
            else{
                TERM t = TERM(std::stold(temp), 0);
                poly_terms.push_back(t);
            }
        }
    }
    else{
        for(int i = 0; i < separation.size(); i++){
            TERM t = TERM(std::stold(separation[i]), separation.size()-i - 1);
            poly_terms.push_back(t);
        }
    }
    POLY p = POLY(poly_terms);
    return p;
}

std::tuple<POLY, POLY> get_poly(std::string& func){
    POLY polycoef, polycoef2;
    std::string prompt, input;
    prompt = "Enter the";
    prompt+=(func=="pow"?" base":func=="log"?" function":"");
    prompt+=" polynomial/coefficients:\n";
    std::cout << prompt;
    input = get_poly();
    polycoef = get_vec(input);
    if(func=="pow"||func=="log"){
        prompt = "Enter the ";
        prompt+=(func=="pow"?"exponent":"base(leave empty or enter 0 to get natural log)");
        prompt+=" polynomial/coefficients:\n";
        std::cout << prompt;
        input = get_poly();
        polycoef2 = get_vec(input);
    }
    return std::tuple<POLY, POLY>(polycoef, polycoef2);
}

long double get_point(const std::string& prompt,const std::string& func,const std::function<long double(long double)>& poly){
    std::cout << prompt;
    long double a;
    while(!(std::cin >> a) || ((func=="sqrt"||func=="pow")&&(poly(a)<0)) || ((func=="log")&&(poly(a)<=0)) || ((func=="acos"||func=="asin")&&(std::abs(poly(a))>1))){
        if(std::cin.fail()) std::cout << "Invalid number";
        else std::cout << "Number is out of domain of function";
        std::cout << " please enter again:\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return a;
}

int get_degree(const std::string& prompt){
    int degree;
    std::cout << prompt;
    while(!(std::cin >> degree) || degree < 0 || degree > 100){
        std::cout << "Invalid number please enter again(valid range: 0-100):\n";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    return degree;
}
