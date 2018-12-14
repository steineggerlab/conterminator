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
    std::vector<MMseqsParameter*> distanceton;
private:
    LocalParameters() :
            Parameters(){



        // strucclust
        conterminatorSearch = removeParameter(searchworkflow, PARAM_MAX_SEQS);
        conterminatorworkflow = combineList(conterminatorworkflow, conterminatorSearch);
        conterminatorworkflow.push_back(&PARAM_REMOVE_TMP_FILES);
        conterminatorworkflow.push_back(&PARAM_RUNNER);
        // distanceton
        distanceton.push_back(&PARAM_EXTRACT_MODE);
        distanceton.push_back(&PARAM_THREADS);
        distanceton.push_back(&PARAM_V);
    }

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};


#endif
