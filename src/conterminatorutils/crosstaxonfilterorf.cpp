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
#include <mmseqs/src/commons/Orf.h>
#include "LocalParameters.h"


#ifdef OPENMP
#include <omp.h>
#endif


int crosstaxonfilterorf(int argc, const char **argv, const Command &command) {
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
    DBReader<unsigned int> orfHeader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    orfHeader.open(DBReader<unsigned int>::NOSORT);
    orfHeader.readMmapedDataInMemory();
    DBReader<unsigned int> reader(par.db3.c_str(), par.db3Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);
    DBWriter writer(par.db4.c_str(), par.db4Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    NcbiTaxonomy * t = NcbiTaxonomy::openTaxonomy(par.db1);

    std::vector<std::string> blacklistStr = Util::split(par.blacklist, ",");
    std::vector<int> blackList;
    for (size_t i = 0; i < blacklistStr.size(); ++i) {
        int currTaxa = Util::fast_atoi<int>(blacklistStr[i].c_str());
        blackList.push_back(currTaxa);
    }
    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        KingdomExpression kingdomExpression(par.kingdoms, *t);
        size_t taxTermCount = kingdomExpression.getTaxTerms().size();
        size_t *taxaCounter = new size_t[taxTermCount];
        char buffer[4096];
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            memset(taxaCounter, 0, taxTermCount * sizeof(size_t));
            unsigned int queryKey = reader.getDbKey(i);
            char *queryHeader = orfHeader.getDataByDBKey(queryKey, thread_idx);
            Orf::SequenceLocation qloc = Orf::parseOrfHeader(queryHeader);
            unsigned int queryTaxon = TaxonUtils::getTaxon(qloc.id, mapping);
            if(queryTaxon == 0 || queryTaxon == UINT_MAX ){
                continue;
            }

            int queryAncestorTermId = kingdomExpression.isAncestorOf(queryTaxon);

            if (queryAncestorTermId == -1) {
                continue;
            }
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getEntryLen(i);

            if (length == 1) {
                continue;
            }
            // find taxonomical information
            char * dataToRead = data;
            int distinctTaxaCnt = 1;
            while(*dataToRead != '\0' && distinctTaxaCnt == 1){
                Util::parseKey(dataToRead, buffer);
                unsigned int targetKey = Util::fast_atoi<unsigned int>(buffer);
                char *targetHeader = orfHeader.getDataByDBKey(targetKey, thread_idx);
                Orf::SequenceLocation tlog = Orf::parseOrfHeader(targetHeader);
                unsigned int targetTaxon = TaxonUtils::getTaxon(tlog.id, mapping);

                int targetAncestorTermId = kingdomExpression.isAncestorOf(targetTaxon);

                if(targetAncestorTermId == -1){
                    dataToRead = Util::skipLine(dataToRead);
                    continue;
                }
                distinctTaxaCnt += (queryAncestorTermId != targetAncestorTermId);
                dataToRead = Util::skipLine(dataToRead);
            }
            if (distinctTaxaCnt > 1) {
                writer.writeData(data, length, queryKey, thread_idx);
            }
        }
        delete [] taxaCounter;
    }

    delete t;
    writer.close();
    reader.close();
    orfHeader.close();
    return EXIT_SUCCESS;
}

