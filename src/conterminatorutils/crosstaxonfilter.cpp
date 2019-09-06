#include "TaxonUtils.h"
#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "filterdb.h"
#include <algorithm>
#include "Matcher.h"
#include "IntervalArray.h"
#include <set>
#include <omptl/omptl_algorithm>
#include <mmseqs/src/taxonomy/TaxonomyExpression.h>


#ifdef OPENMP
#include <omp.h>
#endif


int crosstaxonfilter(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);


    std::vector<std::pair<unsigned int, unsigned int>> mapping;
    if (FileUtil::fileExists(std::string(par.db1 + "_mapping").c_str()) == false) {
        Debug(Debug::ERROR) << par.db1 + "_mapping" << " does not exist. Please create the taxonomy mapping!\n";
        EXIT(EXIT_FAILURE);
    }
    bool isSorted = Util::readMapping(par.db1 + "_mapping", mapping);
    if (isSorted == false) {
        std::stable_sort(mapping.begin(), mapping.end(), ffindexFilter::compareFirstInt());
    }
    std::vector<std::string> ranks = Util::split(par.lcaRanks, ":");

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    NcbiTaxonomy * t = NcbiTaxonomy::openTaxonomy(par.db1);

    Debug(Debug::INFO) << "Add taxonomy information ...\n";

    TaxonomyExpression taxonomyExpression(par.taxonList);
    size_t taxTermCount = taxonomyExpression.getTaxTerms().size();

    std::vector<std::string> blacklistStr = Util::split(par.blacklist, ",");
    std::vector<int> blackList;
    for (size_t i = 0; i < blacklistStr.size(); ++i) {
        int currTaxa = Util::fast_atoi<int>(blacklistStr[i].c_str());
        blackList.push_back(currTaxa);
    }
    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        size_t *taxaCounter = new size_t[taxTermCount];
        std::vector<TaxonUtils::TaxonInformation> elements;
        char buffer[4096];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            elements.clear();
            memset(taxaCounter, 0, taxTermCount * sizeof(size_t));
            unsigned int queryKey = reader.getDbKey(i);

            unsigned int queryTaxon = TaxonUtils::getTaxon(queryKey, mapping);
            if(queryTaxon == 0 || queryTaxon == UINT_MAX ){
                continue;
            }
            unsigned int queryAncestorTermId = UINT_MAX;

            for (size_t j = 0; j < taxTermCount; ++j) {
                int taxIndex = taxonomyExpression.isAncestorOf(*t, queryTaxon);
                if (taxIndex != -1) {
                    queryAncestorTermId = taxIndex;
                    break;
                }
            }

            if (queryAncestorTermId == UINT_MAX) {
                continue;
            }
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getEntryLen(i);

            if (length == 1) {
                continue;
            }
            // find taxonomical information
            TaxonUtils::assignTaxonomy(elements, data, mapping, *t, taxonomyExpression, blackList, taxaCounter);
            std::sort(elements.begin(), elements.end(), TaxonUtils::TaxonInformation::compareByTaxAndStart);
            int distinctTaxaCnt = 0;

            // find max. taxa
            for (size_t taxTermId = 0; taxTermId < taxTermCount; taxTermId++) {
                bool hasTaxa = (taxaCounter[taxTermId] > 0);
                distinctTaxaCnt += hasTaxa;
            }
            if (distinctTaxaCnt > 1) {
                writer.writeStart(thread_idx);
                for(size_t i = 0; i < elements.size(); i++){
                    Matcher::result_t res = Matcher::parseAlignmentRecord(elements[i].data, true);
                    size_t len =Matcher::resultToBuffer(buffer, res, true);
                    writer.writeAdd(buffer, len, thread_idx);
                }
                //res.dbKey=overlappingAlnRes.
                writer.writeEnd(queryKey, thread_idx);
            }
        }
        delete [] taxaCounter;
    }

    delete t;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}

