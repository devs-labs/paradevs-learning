#include <algorithm>
#include <cstdlib>
#include <map>
#include <fstream>
#include <iostream>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include <HMM.hpp>

#include <graph_manager.hpp>
#include <models.hpp>

#include <paradevs/common/RootCoordinator.hpp>

using namespace paradevs::tests::pdevs;
using namespace paradevs::common;

enum IO { IN = 0, OUT = 1, INF = 2, END = 3 };

struct Observation : std::pair < double, IO >
{
    Observation(double value, IO io) : std::pair < double, IO >(value, io)
    { }

    bool operator==(const Observation& o) const
    { return o.first == first and o.second == second; }
};

typedef std::vector < Observation > Sequence;
typedef std::vector < Sequence > Sequences;
typedef std::map < Observation, unsigned int > Table;
typedef std::map < double, unsigned int > DoubleTable;
typedef std::vector < double > Doubles;

Table get_unique(const Sequence& sequence)
{
    Table unique;
    Sequence::const_iterator it = sequence.begin();

    while (it != sequence.end()) {
        if (unique.find(*it) == unique.end()) {
            unique[*it] = 1;
        } else {
            unique[*it]++;
        }
        ++it;
    }
    return unique;
}

DoubleTable get_unique(const Doubles& sequence)
{
    DoubleTable unique;
    Doubles::const_iterator it = sequence.begin();

    while (it != sequence.end()) {
        if (unique.find(*it) == unique.end()) {
            unique[*it] = 1;
        } else {
            unique[*it]++;
        }
        ++it;
    }
    return unique;
}

Sequence read_initial_sequence(const std::string& name, bool ignore = false)
{
    Sequence sequence;
    std::ifstream f("seq");
    char line[100000];
    std::vector < std::string > numbers;

    f.getline(line, 100000);

    std::string l = line;

    boost::split(numbers, l, boost::is_any_of(" \t\n"));

    for (std::vector < std::string >::const_iterator it = numbers.begin();
         it != numbers.end(); ++it) {
        if (*it != "") {
            std::vector < std::string > value;

            boost::split(value, *it, boost::is_any_of("|"));
            if (value[0] == name) {
                if (value[2] != "INF") {
                    if (not ignore) {
                        if (value[1] == "O") {
                            sequence.push_back(
                                Observation(
                                    boost::lexical_cast < double >(
                                        value[2]),
                                    value[1] == "I" ? IN : OUT));
                        }
                    } else {
                        if (value[1] == "O") {
                            sequence.push_back(
                                Observation(
                                    boost::lexical_cast < double >(
                                        value[2]), OUT));
                        }
                    }
                } else {
                    if (not ignore) {
                        sequence.push_back(Observation(0, INF));
                    }
                }
            }
        }
    }
    return sequence;
}

Sequences build_initial_sequences(const std::string& name)
{
    Sequence sequence = read_initial_sequence(name);
    Sequences sequences;
    unsigned int length = sequence.size() / 20;
    Sequence::const_iterator it = sequence.begin();
    unsigned int added = 0;
    unsigned int index = 0;

    sequences.push_back(Sequence());
    while (it != sequence.end()) {
        if (added == length) {
            sequences[index].push_back(Observation(0, END));
            ++index;
            added = 0;
            sequences.push_back(Sequence());
        }
        sequences[index].push_back(*it);
        ++added;
        ++it;
    }
    if (sequence.size() % 20 != 0) {
        sequences[index].push_back(Observation(0, END));
    }
    return sequences;
}

void save_sequences(const std::string& name, const Sequences& sequences)
{
    std::ofstream f((boost::format("seq-%1%.input") %
                         name).str().c_str());

    for (Sequences::const_iterator it = sequences.begin();
         it != sequences.end(); ++it) {
        for (Sequence::const_iterator it2 = it->begin();
             it2 != it->end(); ++it2) {
            f << it2->first << "|" << it2->second << " ";
        }
        f << std::endl;
    }
    f.close();
}

Doubles build_probs(unsigned int size)
{
    Doubles sequence;
    double sum = 0;

    for (unsigned int i = 0; i < size; ++i) {
        double p = rand() % 1000;

        sequence.push_back(p);
        sum += p;
    }
    for (unsigned int i = 0; i < size; ++i) {
        sequence[i] /= sum;
    }
    return sequence;
}

void build_init(const std::string& name, unsigned int state_number)
{
    Sequence sequence = read_initial_sequence(name);

    {
        std::ofstream f((boost::format("seq-%1%-init.emit") %
                         name).str().c_str());
        Table unique = get_unique(sequence);

        for (unsigned int i = 0; i < state_number; ++i) {
            unsigned int j = 0;
            Doubles p = build_probs(unique.size());

            for (Table::const_iterator it = unique.begin();
                 it != unique.end(); ++it) {
                f << "S" << i << "\t"
                  << it->first.first << "|" << it->first.second
                  << "\t" << p[j] << std::endl;
                ++j;
            }
        }
        f << "S" << state_number << "\t0|3\t1" << std::endl;
        f.close();
    }

    {
        std::ofstream f((boost::format("seq-%1%-init.trans") %
                         name).str().c_str());

        f << "INIT" << std::endl;
        f << "INIT\tS0\t1" << std::endl;
        for (unsigned int i = 0; i < state_number + 1; ++i) {
            Doubles p = build_probs(state_number + 1);

            for (unsigned int j = 0; j < state_number + 1; ++j) {
                f << "S" << i << "\t" << "S" << j << "\t" << p[j] << std::endl;
            }
        }
        f << "S" << state_number << "\t" << "FINAL\t1" << std::endl;
        f.close();
    }

    save_sequences(name, build_initial_sequences(name));
}

void learn(const std::string& name, unsigned int maxIterations)
{
    Hmm hmm;
    std::string output = (boost::format("seq-%1%-result") %
                          name).str();
    ifstream istrm((boost::format("seq-%1%.input") % name).str().c_str());
    std:: vector < std::vector < unsigned long >* > trainingSequences;

    hmm.loadProbs((boost::format("seq-%1%-init") % name).str().c_str());
    hmm.readSeqs(istrm, trainingSequences);
    hmm.baumWelch(trainingSequences, maxIterations);
    hmm.saveProbs(output.c_str());
}

void compute_stats(const Doubles& sequence)
{
    unsigned int n = sequence.size();
    DoubleTable D = get_unique(sequence);
    DoubleTable::iterator it = D.begin();
    double m = 0;

    while (it != D.end()) {
        m += it->first * (it->second / (double)(n));

        // std::cout << it->first << "\t"
        //           << (it->second / (double)(n))
        //           << std::endl;
        ++it;
    }

    std::cout << "average = " << m << std::endl;
}

void compute_stats(const Sequence& sequence)
{
    unsigned int n = sequence.size();
    Table D = get_unique(sequence);
    Table::iterator it = D.begin();
    double m = 0;

    while (it != D.end()) {
        m += it->first.first * (it->second / (double)(n));

        // std::cout << it->first.first << "\t"
        //           << (it->second / (double)(n))
        //           << std::endl;
        ++it;
    }

    std::cout << "average = " << m << std::endl;
}

Sequence generate_devs()
{
    Sequence sequence;

    paradevs::common::RootCoordinator <
        DoubleTime, paradevs::pdevs::Coordinator <
            DoubleTime,
                        // GeneratorGraphManager >
        TwoModelsGraphManager >
        > rc(0, 10000, "root", paradevs::common::NoParameters(),
             paradevs::common::NoParameters());

    rc.run();
    return sequence;
}

typedef std::vector < double > Dates;
typedef std::vector < double > Durations;

Dates generate_dates(double max, const Dates& events, Hmm& hmm)
{
    Dates dates;
    double t = 0;
    unsigned long state = hmm.getInitState();
    Dates::const_iterator ite = events.begin();

    while (t < max) {
        unsigned long obs;

        if (hmm.genSeq(state, obs)) {
            std::string obsStr = hmm.getStr(obs);
            int type;
            std::vector < std::string > obsV;

            boost::split(obsV, obsStr, boost::is_any_of("|"));
            type = atoi(obsV[1].c_str());
            if (type != IN) {
                if (type != END) {
                    if (type != INF) {
                        double value = atof(obsV[0].c_str());

                        if (dates.empty()) {
                            dates.push_back(value);
                        } else {
                            dates.push_back(dates.back() + value);
                        }
                        t = dates.back();
                        while (ite != events.end() and t >= *ite) {
                            ++ite;
                        }
                    } else {
                        // double duration = *ite - dates.back();
                        // std::string key =
                        //     (boost::format("%1%|0") % duration).str();

                        // // std::cout << "INF => " << duration << std::endl;
                        // // std::cout << "    => "
                        // //           << key << " "
                        // //           << hmm.getId(key) << " "
                        // //           << std::endl;

                        // OneDTable::const_iterator it =
                        //     hmm._transition[state]->begin();
                        // double max = 0;
                        // unsigned long new_state;

                        // while (it != hmm._transition[state]->end()) {
                        //     OneDTable::const_iterator it2 =
                        //         hmm._emission[it->first]->find(hmm.getId(key));
                        //     double p = std::exp(it2->second) *
                        //         std::exp(it->second);

                        //     if (max < p) {
                        //         max = p;
                        //         new_state = it->first;
                        //     }
                        //     ++it;
                        // }
                        // state = new_state;
                        // dates.push_back(*ite);
                        t = *ite;
                        ++ite;
                    }
                } else {
                    state = hmm.getInitState();
                }
            } else {
                state = hmm.getInitState();
            }
        }
    }
    return dates;
}

Durations compute_durations(const Dates& dates)
{
    Durations durations;
    double previous = 0;
    Dates::const_iterator it = dates.begin();

    while (it != dates.end()) {
        durations.push_back(*it - previous);
        previous = *it;
        ++it;
    }
    return durations;
}

// merge date lists
Dates union_dates(const std::vector < Dates >& D)
{
    Dates U;
    std::vector < Dates >::const_iterator it = D.begin();

    while (it != D.end()) {
        Dates::const_iterator it2 = it->begin();

        while (it2 != it->end()) {
            if (std::find(U.begin(), U.end(), *it2) == U.end()) {
                U.push_back(*it2);
            }
            ++it2;
        }
        ++it;
    }
    std::sort(U.begin(), U.end());
    return U;
}

void generate_hmm(double max, const std::vector < std::string >& names)
{
    std::map < std::string, Hmm > hmms;
    std::vector < std::string >::const_iterator it = names.begin();

    while (it != names.end()) {
        hmms[*it] = Hmm();
        hmms[*it].loadProbs((boost::format("seq-%1%-result") %
                       *it).str().c_str());
        ++it;
    }

    Dates A_dates, B_dates, C_dates, D_dates, E_dates, F_dates, G_dates,
        H_dates;
    Durations A_durations, B_durations, C_durations, D_durations, E_durations,
        F_durations, G_durations, H_durations;

    A_dates = generate_dates(max, {}, hmms["a"]);
    A_durations = compute_durations(A_dates);
    compute_stats(A_durations);

    B_dates = generate_dates(max, {}, hmms["b"]);
    B_durations = compute_durations(B_dates);
    compute_stats(B_durations);

    C_dates = generate_dates(max, {}, hmms["c"]);
    C_durations = compute_durations(C_dates);
    compute_stats(C_durations);

    Dates AB_dates = union_dates({A_dates, B_dates});

    D_dates = generate_dates(max, AB_dates, hmms["d"]);
    D_durations = compute_durations(D_dates);
    compute_stats(D_durations);

    Dates BC_dates = union_dates({B_dates, C_dates});

    E_dates = generate_dates(max, BC_dates, hmms["e"]);
    E_durations = compute_durations(E_dates);
    compute_stats(E_durations);

    Dates DE_dates = union_dates({D_dates, E_dates});

    F_dates = generate_dates(max, DE_dates, hmms["f"]);
    F_durations = compute_durations(F_dates);
    compute_stats(F_durations);

    G_dates = generate_dates(max, E_dates, hmms["g"]);
    G_durations = compute_durations(G_dates);
    compute_stats(G_durations);

    Dates FG_dates = union_dates({F_dates, G_dates});

    H_dates = generate_dates(max, FG_dates, hmms["h"]);
    H_durations = compute_durations(H_dates);
    compute_stats(H_durations);
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout <<
            "USAGE: paradevs-learning [-D | -L x n i | -G n ... | -M x | -A]"
                  << std::endl;
        return -1;
    }
    if (strcmp(argv[1], "-D") == 0) {
        generate_devs();
    } else if (strcmp(argv[1], "-L") == 0) {
        if (argc != 5) {
            std::cout <<
                "USAGE: paradevs-learning [-D | -L x n i | -G n ... | -M x | -A n i]"
                      << std::endl;
            return -1;
        }
        build_init(argv[2], atoi(argv[3]));
        learn(argv[2], atoi(argv[4]));
    } else if (strcmp(argv[1], "-G") == 0) {
        if (argc < 3) {
            std::cout <<
                "USAGE: paradevs-learning [-D | -L x n i | -G n ... | -M x | -A n i]"
                      << std::endl;
            return -1;
        }
        std::vector < std::string > names;

        for (int i = 0; i < argc - 3; ++i) {
            names.push_back(argv[i + 3]);
        }
        generate_hmm(atoi(argv[2]), names);
    } else if (strcmp(argv[1], "-M") == 0) {
        if (argc != 3) {
            std::cout <<
                "USAGE: paradevs-learning [-D | -L x n i | -G n ... | -M x | -A n i]"
                      << std::endl;
            return -1;
        }
        compute_stats(read_initial_sequence(argv[2], true));
    }
    return 0;
}
