#include <iostream>
#include <cassert>
#include <vector>
#include <utility> 
#include <cmath>
#include "problib.h"

std::vector<std::vector<double>>v={{0.8,0.2},{0.6,0.4}};
std::vector<std::vector<double>>w={{0,1,0,0},{0.5,0,0.5,0},{0,1,0,0},{0,0,0,1}};

void testIsIrreducible(){
    try{
        ProbLib::MarkovChain mc(v);
        bool value = mc.is_irreducible(v);
        assert(value == true);
        std::cout<<"MC is_irreducible test passed: " << value << std::endl;;
    }catch(const std::exception& e){
        std::cerr << "MC is_irreducible test failed: "<< std::endl;
    }

    try{
        ProbLib::MarkovChain mc(w);
        bool value = mc.is_irreducible(w);
        assert(value == false);
        std::cout<<"MC is_irreducible test passed: " << value << std::endl;;
    }catch(const std::exception& e){
        std::cerr << "MC is_irreducible test failed: "<< e.what() << std::endl;
    }
}

void testIsAperiodic(){
    try{
        ProbLib::MarkovChain mc(v);
        bool value = mc.is_aperiodic(v);
        assert(value == true);
        std::cout<<"MC is_aperiodic test passed: " << value << std::endl;;
    }catch(const std::exception& e){
        std::cerr << "MC is_aperiodic test failed: "<< e.what() << std::endl;
    }

    try{
        ProbLib::MarkovChain mc(w);
        bool value = mc.is_aperiodic(w);
        assert(value == false);
        std::cout<<"MC is_aperiodic test passed: " << value << std::endl;;
    }catch(const std::exception& e){
        std::cerr << "MC is_aperiodic test failed: "<< e.what() << std::endl;
    }
}

void testIsErgodic(){
    try{
        ProbLib::MarkovChain mc(v);
        bool value = mc.is_ergodic(v);
        assert(value == true);
        std::cout<<"MC is_ergodic test passed: " << value << std::endl;;
    }catch(const std::exception& e){
        std::cerr << "MC is_ergodic test failed: "<< e.what() << std::endl;
    }

    try{
        ProbLib::MarkovChain mc(w);
        bool value = mc.is_ergodic(w);
        assert(value == false);
        std::cout<<"MC is_ergodic test passed: " << value << std::endl;;
    }catch(const std::exception& e){
        std::cerr << "MC is_ergodic test failed: "<< e.what() << std::endl;
    }
}

void testSteadyState(){
    std::vector<std::vector<double>>ergSS={{0.75,0.25}};
    try{
        ProbLib::MarkovChain mc(v);
        std::vector<std::vector<double>> ss=mc.steady_state(v).second;
        bool test;
        int i=0;
        int j=0;
        for(const auto& row : ss) {
            for(double val : row) {
                test=(fabs(val-ergSS[i][j])<1e-10);
                if(!test) break;
                j++;
            }
            i++;
        }
        assert(test);
        std::cout<<"MC Steady State test passed"<<std::endl;
    }catch(const std::exception& e){
        std::cerr << "MC Steady State test failed: "<< e.what() << std::endl;
    }
}

void testnStep(){
    std::vector<std::vector<double>>nS={{0.76,0.24},{0.72,0.28}};
    try{
        ProbLib::MarkovChain mc(v);
        std::vector<std::vector<double>>nStep=mc.n_step_transition(2);
        bool test;
        int i=0;
        int j=0;
        for(const auto& row : nStep) {
            for(double val : row) {
                test=(fabs(val-nS[i][j])<1e-10);
                if(!test) break;
                j++;
            }
            i++;
            j=0;
        }
        assert(test);
        std::cout<<"MC n-Step test passed"<<std::endl;
    }catch(const std::exception& e){
        std::cerr << "MC n-Step test failed: "<< e.what() << std::endl;
    }
}


int main(){

    testIsIrreducible();
    testIsAperiodic();
    testIsErgodic();
    testSteadyState();
    testnStep();

    return 0;
}