#ifndef _OPTIMIZER_
#define _OPTIMIZER_

#include <iostream>

#include "soga.h"
#include "Chromosome.h"

using std::cout;
using std::endl;

template <typename GA1, typename GA2>
class Optimizer
{
private:
    GA1 ga1;
    GA2 ga2;
    Individual optimum1, optimum2;
    Statistics stats1, stats2;


public:
    Optimizer(Vec2D const minPos, Vec2D const maxPos, string const& filePath1, string const& filePath2, size_t popSize, vector<string> const & buildList1, vector<string> const & buildList2)
        : ga1(minPos, maxPos, filePath1, filePath2, popSize, buildList1, buildList2),
          ga2(minPos, maxPos, filePath2, filePath1, popSize, buildList2, buildList1)
    {}

    void optimize(size_t const sel, size_t const co, size_t const mut, size_t const iterations, size_t const genPerIt)
    {
        setSelection(sel);
        setCrossover(co);
        setMutation(mut);

        cout << "Performing Optimization with the following parameters" << endl << endl;
        cout << "Selection Operator: " << getSelectionOperatorName() << endl;
        cout << "Crossover Operator: " << getCrossoverOperatorName() << endl;
        cout << "Mutation Operator: " << getMutationOperatorName() << endl;
//        cout << "Number of Iterations: " << iterations << endl;
//        cout << "Generations per Iteration: " << genPerIt << endl;

        for(size_t i = 0; i < iterations; ++i)
        {
//            std::cout << "Progress: " << static_cast<double>(i)/iterations*100 << "%" << "\r" << std::flush;
//            printf("%c[2K", 27);
            ga1.optimize(optimum2.chromosome, genPerIt);
            ga2.optimize(optimum1.chromosome, genPerIt);
            stats1 = ga1.getStatistics();
            stats2 = ga2.getStatistics();
            optimum1 = stats1.optimum;
            optimum2 = stats2.optimum;
        }
//        ga1.optimize(optimum2.chromosome, 0);
//        ga2.optimize(optimum1.chromosome, 0);
//        stats1 = ga1.getStatistics();
//        stats2 = ga2.getStatistics();
//        optimum1 = stats1.optimum;
//        optimum2 = stats2.optimum;

        printStatistics();
        cout << "\n" << endl;
    }

    void printStatistics()
    {
        string separator;
        for(int i = 0; i < 50; ++i) separator += "-";
        separator += '\n';
        cout << separator << "Statistics Player 1:" << endl;
        stats1.print();
        cout << separator << endl;
        cout << separator << "Statistics Player 2:" << endl;
        stats2.print();
        cout << separator << endl;
    }

    GA1 const& getGA1() const
    {
        return ga1;
    }
    GA2 const& getGA2() const
    {
        return ga2;
    }

    size_t getNumberOfSelectionOperators()
    {
        ga1.getNumberOfSelectionOperators();
    }

    size_t getNumberOfCrossoverOperators()
    {
        return ga1.getNumberOfCrossoverOperators();
    }

    size_t getNumberOfMutationOperators()
    {
        return ga1.getNumberOfMutationOperators();
    }



    void setSelection(size_t const value)
    {
        ga1.setSelection(value);
        ga2.setSelection(value);
    }

    void setCrossover(size_t const value)
    {
        ga1.setCrossover(value);
        ga2.setCrossover(value);
    }

    void setMutation(size_t const value)
    {
        ga1.setMutation(value);
        ga2.setMutation(value);
    }

    string getSelectionOperatorName()
    {
        return ga1.getSelectionOperatorName();
    }
    string getCrossoverOperatorName()
    {
        return ga1.getCrossoverOperatorName();
    }
    string getMutationOperatorName()
    {
        return ga1.getMutationOperatorName();
    }
};

#endif
