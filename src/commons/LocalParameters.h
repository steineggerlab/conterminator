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

    std::vector<MMseqsParameter*> conterminatorworkflow;
    std::vector<MMseqsParameter*> conterminatorSearch;
    std::vector<MMseqsParameter*> extractalignments;
    std::vector<MMseqsParameter*> createstats;
    std::vector<MMseqsParameter*> distanceton;
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
        // distanceton
        distanceton.push_back(&PARAM_EXTRACT_MODE);
        distanceton.push_back(&PARAM_THREADS);
        distanceton.push_back(&PARAM_V);
        // conterminator
        conterminatorSearch = removeParameter(searchworkflow, PARAM_MAX_SEQS);
        conterminatorworkflow = combineList(conterminatorworkflow, conterminatorSearch);
        conterminatorworkflow = combineList(conterminatorworkflow, createstats);
        conterminatorworkflow = combineList(conterminatorworkflow, extractalignments);
    }

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
