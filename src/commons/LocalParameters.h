#ifndef LOCALPARAMETERS_H
#define LOCALPARAMETERS_H

#include <Parameters.h>

class LocalParameters : public Parameters {
public:
    static void initInstance() {
        new LocalParameters;
    }

    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initInstance();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }
    PARAMETER(PARAM_KINGDOMS)
    std::string kingdoms;

    std::vector<MMseqsParameter*> conterminatordna;
    std::vector<MMseqsParameter*> conterminatorprotein;
    std::vector<MMseqsParameter*> extractalignments;
    std::vector<MMseqsParameter*> createstats;
private:
    LocalParameters() :
            Parameters(),
            PARAM_KINGDOMS(PARAM_KINGDOMS_ID,"--kingdoms", "Compare across kingdoms", "",typeid(std::string), (void *) &kingdoms, "[,]"){

        // extractalignments
        extractalignments.push_back(&PARAM_BLACKLIST);
        extractalignments.push_back(&PARAM_KINGDOMS);
        extractalignments.push_back(&PARAM_THREADS);
        extractalignments.push_back(&PARAM_V);
        // createstats
        createstats.push_back(&PARAM_BLACKLIST);
        createstats.push_back(&PARAM_KINGDOMS);
        createstats.push_back(&PARAM_THREADS);
        createstats.push_back(&PARAM_V);
        // conterminatordna
        conterminatordna = removeParameter(searchworkflow, PARAM_MAX_SEQS);
        conterminatordna = combineList(conterminatordna, createdb);
        conterminatordna = combineList(conterminatordna, createtaxdb);
        conterminatordna = combineList(conterminatordna, createstats);
        conterminatordna = combineList(conterminatordna, extractalignments);
        // conterminatorprotein
        conterminatorprotein = removeParameter(linclustworkflow, PARAM_MAX_SEQS);
        conterminatorprotein = combineList(conterminatordna, createdb);
        conterminatorprotein = combineList(conterminatordna, createtaxdb);
        conterminatorprotein = combineList(conterminatordna, createstats);
        conterminatorprotein = combineList(conterminatordna, extractalignments);
    }
    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
