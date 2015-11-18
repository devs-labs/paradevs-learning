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

// construit une table avec les différentes valeurs d'une séquence et
// compte le nombre d'itérations
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

// construit une table avec les différentes valeurs d'une liste de
// valeurs réelles et compte le nombre d'itérations
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

// relecture du fichier de séquence initiale (seq)
// chaque élément de la séquence est séparé par un espace
// chaque élément est de la forme x|y|z où x est le nom du modèle, y
// est le type d'observation (O ou I) et z la durée depuis le dernier
// événement (ou INF si la durée est l'état est infinie)
// un objet Sequence est retourné avec la liste des observations
// si le paramètre ignore est true, les INF sont ignorés ; par défaut,
// ce n'est pas le cas
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

// la séquence initiale est lue et découpée en sous-séquences de
// longueur 20 ;
// chaque sous-séquence est terminée par une observation END pour que
// le HMM généré possède un état terminal
Sequences build_initial_sequences(const std::string& name)
{
    Sequence sequence = read_initial_sequence(name, false);
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

// fonction de sauvegarde de la liste des sous-séquences dans les
// fichiers .input
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

// construction un vecteur d'une taille donnée avec des probabilités aléatoires
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

// construction des matrices initiales emit et trans avec les
// probabilités aléatoires et en fonction du nombre d'états cachés
// les matrices sont aussi en fonction des symboles des séquences
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

// construit les matrices emit et trans à l'aide de l'algorithme Baum-Welch
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

// calcule une moyenne à partir d'une liste de réels
// la liste est parcourue pour identifier les différentes valeurs réelles
// puis calcul de la moyenne pondérée
void compute_stats(const Doubles& sequence)
{
    unsigned int n = sequence.size();
    DoubleTable D = get_unique(sequence);
    DoubleTable::iterator it = D.begin();
    double m = 0;

    while (it != D.end()) {
        m += it->first * (it->second / (double)(n));
        ++it;
    }

    std::cout << m << std::endl;
}

// calcule une moyenne à partir d'une séquence
// la séquence est parcourue pour identifier les différentes durées
// puis calcul de la moyenne pondérée
void compute_stats(const Sequence& sequence)
{
    unsigned int n = sequence.size();
    Table D = get_unique(sequence);
    Table::iterator it = D.begin();
    double m = 0;

    while (it != D.end()) {
        m += it->first.first * (it->second / (double)(n));
        ++it;
    }

    std::cout << "average = " << m << std::endl;
}

// simule le graphe de modèles
// ATTENTION !!! TOTALEMENT SPECIFIQUE AU GRAPHE DE MODELES DE
// L'EXEMPLE !!!
Sequence generate_devs(double max)
{
    Sequence sequence;

    paradevs::common::RootCoordinator <
        DoubleTime, paradevs::pdevs::Coordinator <
            DoubleTime,
                        // GeneratorGraphManager >
                        TwoModelsGraphManager >
        > rc(0, max, "root", paradevs::common::NoParameters(),
             paradevs::common::NoParameters());

    rc.run();
    return sequence;
}

typedef std::vector < double > Dates;
typedef std::vector < double > Durations;

// génére les dates des événements de sortie en fonction du HMM sur
// une période [0, max]
// le paramètre events est la liste des dates d'arrivée des événements d'entrée
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
                    // ce n'est pas un état infini
                    if (type != INF) {
                        double value = atof(obsV[0].c_str());

                        if (dates.empty()) {
                            dates.push_back(value);
                        } else {
                            dates.push_back(dates.back() + value);
                        }
                        t = dates.back();
                        // détermine la prochaine date d'arrivée d'un
                        // événement d'entrée
                        while (ite != events.end() and t >= *ite) {
                            ++ite;
                        }
                    } else {
                        // c'est un état infini
                        if (ite != events.end()) {
                            t = *ite;
                            ++ite;
                        }
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

// à partir d'une liste de dates, on calcule une liste de durées
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

// fusionne une liste de listes de dates en une unique liste ordonée
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

// simule les n modèles et calcule les moyennes
// ATTENTION !!! TOTALEMENT SPECIFIQUE AU GRAPHE DE MODELES DE
// L'EXEMPLE !!!
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

    {
        Dates A_dates, B_dates, C_dates, D_dates, E_dates, F_dates, G_dates,
            H_dates;
        Durations A_durations, B_durations, C_durations, D_durations,
            E_durations, F_durations, G_durations, H_durations;

        A_dates = generate_dates(max, {}, hmms["a"]);
        A_durations = compute_durations(A_dates);
        std::cout << "a: ";
        compute_stats(A_durations);

        B_dates = generate_dates(max, {}, hmms["b"]);
        B_durations = compute_durations(B_dates);
        std::cout << "b: ";
        compute_stats(B_durations);

        C_dates = generate_dates(max, {}, hmms["c"]);
        C_durations = compute_durations(C_dates);
        std::cout << "c: ";
        compute_stats(C_durations);

        Dates AB_dates = union_dates({A_dates, B_dates});

        D_dates = generate_dates(max, AB_dates, hmms["d"]);
        D_durations = compute_durations(D_dates);
        std::cout << "d: ";
        compute_stats(D_durations);

        Dates BC_dates = union_dates({B_dates, C_dates});

        E_dates = generate_dates(max, BC_dates, hmms["e"]);
        E_durations = compute_durations(E_dates);
        std::cout << "e: ";
        compute_stats(E_durations);

        Dates DE_dates = union_dates({D_dates, E_dates});

        F_dates = generate_dates(max, DE_dates, hmms["f"]);
        F_durations = compute_durations(F_dates);
        std::cout << "f: ";
        compute_stats(F_durations);

        G_dates = generate_dates(max, E_dates, hmms["g"]);
        G_durations = compute_durations(G_dates);
        std::cout << "g: ";
        compute_stats(G_durations);

        Dates FG_dates = union_dates({F_dates, G_dates});

        H_dates = generate_dates(max, FG_dates, hmms["h"]);
        H_durations = compute_durations(H_dates);
        std::cout << "h: ";
        compute_stats(H_durations);
    }

    // {
    //     Dates A_dates, B_dates, C_dates, D_dates, E_dates, F_dates, G_dates,
    //         H_dates;

    //     A_dates = generate_dates(max, {}, hmms["a"]);
    //     std::cout << "a: ";
    //     compute_stats(compute_durations(A_dates));

    //     B_dates = generate_dates(max, {}, hmms["a"]);
    //     std::cout << "b: ";
    //     compute_stats(compute_durations(B_dates));

    //     C_dates = generate_dates(max, {}, hmms["a"]);
    //     std::cout << "c: ";
    //     compute_stats(compute_durations(C_dates));

    //     Dates AB_dates = union_dates({A_dates, B_dates});

    //     D_dates = generate_dates(max, AB_dates, hmms["d"]);
    //     std::cout << "d: ";
    //     compute_stats(compute_durations(D_dates));

    //     Dates BC_dates = union_dates({B_dates, C_dates});

    //     E_dates = generate_dates(max, BC_dates, hmms["d"]);
    //     std::cout << "e: ";
    //     compute_stats(compute_durations(E_dates));

    //     Dates DE_dates = union_dates({D_dates, E_dates});

    //     F_dates = generate_dates(max, DE_dates, hmms["d"]);
    //     std::cout << "f: ";
    //     compute_stats(compute_durations(F_dates));

    //     G_dates = generate_dates(max, E_dates, hmms["g"]);
    //     std::cout << "g: ";
    //     compute_stats(compute_durations(G_dates));

    //     Dates FG_dates = union_dates({F_dates, G_dates});

    //     H_dates = generate_dates(max, FG_dates, hmms["d"]);
    //     std::cout << "h: ";
    //     compute_stats(compute_durations(H_dates));
    // }
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
        if (argc != 3) {
            std::cout <<
                "USAGE: paradevs-learning [-D | -L x n i | -G n ... | -M x | -A n i]"
                      << std::endl;
            return -1;
        }
        generate_devs(atof(argv[2]));
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
