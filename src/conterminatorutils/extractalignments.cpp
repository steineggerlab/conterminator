#include "TaxonUtils.h"
#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include <algorithm>
#include "Matcher.h"
#include "IntervalArray.h"
#include <set>
#include <omptl/omptl_algorithm>
#include "LocalParameters.h"


#ifdef OPENMP
#include <omp.h>
#endif


int extractalignments(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);


    std::vector<std::pair<unsigned int, unsigned int>> mapping;
    if (FileUtil::fileExists(std::string(par.db1 + "_mapping").c_str()) == false) {
        Debug(Debug::ERROR) << par.db1 + "_mapping" << " does not exist. Please create the taxonomy mapping!\n";
        EXIT(EXIT_FAILURE);
    }
    bool isSorted = Util::readMapping(par.db1 + "_mapping", mapping);
    if (isSorted == false) {
        std::stable_sort(mapping.begin(), mapping.end(), TaxonUtils::compareToFirstInt);
    }
    std::vector<std::string> ranks = Util::split(par.lcaRanks, ":");

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, reader.getDbtype());
    writer.open();

    NcbiTaxonomy * t = NcbiTaxonomy::openTaxonomy(par.db1);

    Debug(Debug::INFO) << "Add taxonomy information ...\n";


    std::vector<std::string> blacklistStr = Util::split(par.blacklist, ",");
    std::vector<int> blackList;
    for (size_t i = 0; i < blacklistStr.size(); ++i) {
        int currTaxa = Util::fast_atoi<int>(blacklistStr[i].c_str());
        blackList.push_back(currTaxa);
    }

    struct Contamination{
        Contamination(unsigned int key, int start, int end, unsigned int len)
                : key(key), start(start), end(end), len(len) {}
        Contamination(){};
        unsigned int key;
        int start;
        int end;
        unsigned int len;

        // need for sorting the results
        static bool compareContaminationByKeyStartEnd(const Contamination &first, const Contamination &second) {
            //return (first.eval < second.eval);
            if(first.key < second.key )
                return true;
            if(second.key < first.key )
                return false;
            if(first.start < second.start )
                return true;
            if(second.start < first.start )
                return false;
            if(first.end < second.end )
                return true;
            if(second.end < first.end )
                return false;
            return false;
        }
    };

    std::vector<Contamination> allContaminations;
    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        KingdomExpression kingdomExpression(par.kingdoms, *t);
        size_t taxTermCount = kingdomExpression.getTaxTerms().size();
        std::vector<Contamination> privateContaminations;
        size_t *taxaCounter = new size_t[taxTermCount];
        IntervalArray ** speciesRanges = new IntervalArray*[taxTermCount];
        for(size_t i = 0; i < taxTermCount; i++){
            speciesRanges[i] = new IntervalArray();
        }
        std::vector<TaxonUtils::TaxonInformation> elements;

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
            unsigned int queryLen = reader.getSeqLen(i);

            unsigned int queryTaxon = TaxonUtils::getTaxon(queryKey, mapping);
            if(queryTaxon == 0 || queryTaxon == UINT_MAX ){
                continue;
            }
            unsigned int queryAncestorTermId = UINT_MAX;

            for (size_t j = 0; j < taxTermCount; ++j) {
                int taxIndex = kingdomExpression.isAncestorOf(queryTaxon);
                if (taxIndex != -1) {
                    queryAncestorTermId = taxIndex;
                    break;
                }
            }

            if (queryAncestorTermId == UINT_MAX) {
                continue;
            }
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getSeqLen(i);

            if (length == 1) {
                continue;
            }
            // find taxonomical information
            TaxonUtils::assignTaxonomy(elements, data, mapping, *t, kingdomExpression, blackList, taxaCounter);
            std::sort(elements.begin(), elements.end(), TaxonUtils::TaxonInformation::compareByTaxAndStart);
            int distinctTaxaCnt = 0;

            // find max. taxa
            for (size_t taxTermId = 0; taxTermId < taxTermCount; taxTermId++) {
                bool hasTaxa = (taxaCounter[taxTermId] > 0);
                distinctTaxaCnt += hasTaxa;
            }
            if (distinctTaxaCnt > 1) {
                for (size_t i = 0; i < taxTermCount; i++) {
                    speciesRanges[i]->reset();
                }

                // fill up interval tree with elements
                for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                    if (static_cast<unsigned int>(elements[elementIdx].termId) != queryAncestorTermId) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                        speciesRanges[elements[elementIdx].termId]->insert(res.qStartPos, res.qEndPos);
                    }
                }
                for (size_t i = 0; i < taxTermCount; i++) {
                    speciesRanges[i]->buildRanges();
                }

                for (size_t i = 0; i < taxTermCount; i++) {
                    for (size_t j = 0; j < speciesRanges[i]->getRangesSize(); j++) {
                        IntervalArray::Range range = speciesRanges[i]->getRange(j);
                        privateContaminations.push_back(Contamination(queryKey, range.start, range.end, queryLen));
                    }
                }
            }
        }

#pragma omp critical
        {
            allContaminations.insert(allContaminations.end(), privateContaminations.begin(), privateContaminations.end());
        };
        delete[] taxaCounter;

        for(size_t i = 0; i < taxTermCount; i++){
            delete speciesRanges[i];
        }
        delete [] speciesRanges;
    }

    omptl::sort(allContaminations.begin(), allContaminations.end(), Contamination::compareContaminationByKeyStartEnd);
    size_t writePos = -1;
    unsigned int prevKey = UINT_MAX;

    for(size_t i = 0; i < allContaminations.size(); i++){
        if(prevKey != allContaminations[i].key){
            writePos++;
            allContaminations[writePos]=allContaminations[i];
        }else{
            if(allContaminations[i].start <= (allContaminations[writePos].end+1)){
                allContaminations[writePos].end = std::max(allContaminations[i].end, allContaminations[writePos].end);
            }else{
                writePos++;
                allContaminations[writePos]=allContaminations[i];
            }
        }
        prevKey = allContaminations[i].key;
    }
    writePos++;
    allContaminations.resize (writePos);

    {
        char buffer[4096];

        for(size_t i = 0; i < allContaminations.size(); i++) {
            Matcher::result_t res(allContaminations[i].key, 255,
                                  0.0, 0.0,
                                  1.0, 0.0,
                                  0,
                                  allContaminations[i].start,
                                  allContaminations[i].end,
                                  allContaminations[i].len,
                                  allContaminations[i].start,
                                  allContaminations[i].end,
                                  allContaminations[i].len, "");
            size_t len = Matcher::resultToBuffer(buffer, res, false, false);
            //res.dbKey=overlappingAlnRes.
            writer.writeData(buffer, len, allContaminations[i].key, 0);
        }
    }
    delete t;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}

