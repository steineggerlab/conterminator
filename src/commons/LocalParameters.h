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

    std::vector<MMseqsParameter*> conterminatordna;
    std::vector<MMseqsParameter*> conterminatorprotein;
    std::vector<MMseqsParameter*> extractalignments;
    std::vector<MMseqsParameter*> createstats;
private:
    LocalParameters() :
            Parameters(){

        // extractalignments
        extractalignments.push_back(&PARAM_BLACKLIST);
        extractalignments.push_back(&PARAM_TAXON_LIST);
        extractalignments.push_back(&PARAM_THREADS);
        extractalignments.push_back(&PARAM_V);
        // createstats
        createstats.push_back(&PARAM_BLACKLIST);
        createstats.push_back(&PARAM_TAXON_LIST);
        createstats.push_back(&PARAM_THREADS);
        createstats.push_back(&PARAM_V);
        // conterminatordna
        conterminatordna = removeParameter(searchworkflow, PARAM_MAX_SEQS);
        conterminatordna = combineList(conterminatordna, createstats);
        conterminatordna = combineList(conterminatordna, extractalignments);
        // conterminatorprotein
        conterminatorprotein = removeParameter(linclustworkflow, PARAM_MAX_SEQS);
        conterminatorprotein = combineList(conterminatordna, createstats);
        conterminatorprotein = combineList(conterminatordna, extractalignments);
    }

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
