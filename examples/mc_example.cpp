#include "problib.h"
#include <iostream>

int main() {
    // Example transition matrix
    try{
        std::vector<std::vector<double>> transition_matrix = {
            {0, 0.5, 0, 0.5},
            {0.5, 0, 0.5, 0},
            {0, 0.5, 0, 0.5},
            {0.5, 0, 0.5, 0}
        };

        // Create a MarkovChain object
        ProbLib::MarkovChain mc(transition_matrix);

        //Check if the Markov Chain is irreducible
        if (mc.is_irreducible(transition_matrix)) {
            std::cout << "The Markov chain is irreducible." << std::endl;
        } else {
            std::cout << "The Markov chain is not irreducible." << std::endl;
        }

        //Check if the MarkovChain is aperiodic
        if (mc.is_aperiodic(transition_matrix)) {
            std::cout << "The Markov chain is aperiodic." << std::endl;
        } else {
            std::cout << "The Markov chain is not aperiodic." << std::endl;
        }

        // Check if the Markov chain is ergodic
        if (mc.is_ergodic(transition_matrix)) {
            std::cout << "The Markov chain is ergodic." << std::endl;
        } else {
            std::cout << "The Markov chain is not ergodic." << std::endl;
        }

        // Calculate the 2-step transition matrix
        std::cout << "2-step transition matrix:" << std::endl;
        auto result_matrix = mc.n_step_transition(2);
        for (const auto& row : result_matrix) {
            for (double val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }

        // Compute the steady-state distribution
        auto result = mc.steady_state(transition_matrix);
        std::cout << "Message: " << result.first << std::endl;
        for (const auto& row : result.second) {
            for (double val : row) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    } catch (const std::exception& e){
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
