#include "holder.h"

int main(){
    start_text();
    std::string function = get_function();
    auto [poly1, poly2] = get_poly(function);
    long double point = get_point(function, poly1);
    int degree = get_degree(function);
}
