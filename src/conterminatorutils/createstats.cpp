#include "TaxonUtils.h"
#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include <omptl/omptl_algorithm>
#include <set>
#include <limits>

#ifdef OPENMP
#include <omp.h>
#endif


int createstats(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    // bacteria, archaea, eukaryotic, virus
    par.parseParameters(argc, argv, command, true, 0, 0);

    NcbiTaxonomy * t = NcbiTaxonomy::openTaxonomy(par.db1);
    TaxonomyExpression taxonomyExpression(par.taxonList);
    size_t taxTermCount = taxonomyExpression.getTaxTerms().size();

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

//    DBReader<unsigned int> seqdb(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
//    seqdb.open(DBReader<unsigned int>::NOSORT);
//    seqdb.close();

    DBReader<unsigned int> header(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    header.open(DBReader<unsigned int>::NOSORT);
    header.readMmapedDataInMemory();


    DBReader<unsigned int> sequences(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX);
    sequences.open(DBReader<unsigned int>::NOSORT);

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    std::vector<std::string> blacklistStr = Util::split(par.blacklist, ",");
    std::vector<int> blackList;
    for (size_t i = 0; i < blacklistStr.size(); ++i) {
        int currTaxa = Util::fast_atoi<int>(blacklistStr[i].c_str());
        blackList.push_back(currTaxa);
    }

    Debug::Progress progress(reader.getSize());
    size_t *totalContermCounter = new size_t[taxTermCount];
    memset(totalContermCounter, 0, taxTermCount * sizeof(size_t));
#pragma omp parallel
    {
        std::string resultData;
        resultData.reserve(4096);
        size_t *taxaCounter = new size_t[taxTermCount];
        unsigned int thread_idx = 0;
        std::vector<TaxonUtils::TaxonInformation> elements;
        std::set<unsigned int> idDetected;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            resultData.clear();
            elements.clear();
            idDetected.clear();
            unsigned int queryKey = reader.getDbKey(i);

            unsigned int queryTaxon = TaxonUtils::getTaxon(queryKey, mapping);
            if (queryTaxon == 0 || queryTaxon == UINT_MAX) {
                continue;
            }
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getSeqLens(i);

            if (length == 1) {
                continue;
            }
            // find taxonomical information
            memset(taxaCounter, 0, taxTermCount * sizeof(size_t));
            TaxonUtils::assignTaxonomy(elements, data, mapping, *t, taxonomyExpression, blackList, taxaCounter);
            // recount
            memset(taxaCounter, 0, taxTermCount * sizeof(size_t));
            for(size_t i = 0; i < elements.size(); i++){
                char dbKeyBuffer[255 + 1];
                Util::parseKey(elements[i].data, dbKeyBuffer);
                const unsigned int key = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const bool existsAlready = idDetected.find(key) != idDetected.end();
                if(existsAlready == false) {
                    taxaCounter[elements[i].termId]++;
                    idDetected.insert(key);
                }
            }
            std::sort(elements.begin(), elements.end(), TaxonUtils::TaxonInformation::compareByTaxAndStart);
            // find max. taxa
            size_t minTaxCnt = std::numeric_limits<size_t>::max();
            int minTaxTerm = 0;
            int maxTaxId = 0;
            size_t maxTaxCnt = 0;
            for (size_t taxTermId = 0; taxTermId < taxTermCount; taxTermId++) {
                resultData.append(SSTR(taxaCounter[taxTermId]));
                if(taxaCounter[taxTermId] < minTaxCnt && taxaCounter[taxTermId] != 0){
                    minTaxCnt = taxaCounter[taxTermId];
                    minTaxTerm = taxTermId;
                }

                if(taxaCounter[taxTermId] >= maxTaxCnt){
                    maxTaxCnt = taxaCounter[taxTermId];
                    maxTaxId = taxTermId;
                }
                resultData.push_back('\t');
            }
            if(minTaxTerm == 0 || maxTaxId == 0){
                continue;
            }
            __sync_fetch_and_add(&(totalContermCounter[minTaxTerm]), 1);

            std::set<std::pair<unsigned int, int> > minDbKeys;
            std::map<int, size_t> maxTaxon;
            unsigned int maxDbKey=UINT_MAX;
            for(size_t i = 0; i < elements.size(); i++){
                if(elements[i].termId == minTaxTerm){
                    std::pair<unsigned int, int> pair = std::make_pair(elements[i].dbKey, elements[i].currTaxa);
                    minDbKeys.insert(pair);
                }
                if(elements[i].termId == maxTaxId) {
                    maxDbKey = std::min(elements[i].dbKey, maxDbKey);
                    maxTaxon[elements[i].currTaxa]++;
                }
            }
            std::map<int, size_t>::iterator maxTaxonIt;
            maxTaxCnt = 0;
            for ( maxTaxonIt = maxTaxon.begin(); maxTaxonIt != maxTaxon.end(); maxTaxonIt++ )
            {
                if(maxTaxonIt->second >= maxTaxCnt){
                    maxTaxId = maxTaxonIt->first;
                    maxTaxCnt = maxTaxonIt->second;
                }
            }

            std::set<std::pair<unsigned int, int>>::iterator minDbKeysIt = minDbKeys.begin();
            std::string taxons;
            std::string dbLength;
            unsigned int outDBKey = UINT_MAX;
            for (;minDbKeysIt != minDbKeys.end(); minDbKeysIt++){
                unsigned int dbKey = minDbKeysIt->first;
                outDBKey = std::min(dbKey, outDBKey);
                int taxId = minDbKeysIt->second;
                resultData.append(Util::parseFastaHeader(header.getDataByDBKey(dbKey, thread_idx)));
                dbLength.append(SSTR(sequences.getSeqLens(sequences.getId(dbKey))));
                const TaxonNode* node= t->taxonNode(taxId, false);
                if(node == NULL){
                    taxons.append("Undef");
                }else{
                    taxons.append(node->name);
                }
                taxons.push_back(',');
                dbLength.push_back(',');
                resultData.push_back(',');
            }
            taxons.pop_back();
            dbLength.pop_back();
            resultData.pop_back();
            resultData.push_back('\t');
            resultData.append(taxons);
            resultData.push_back('\t');
            resultData.append(dbLength);
            resultData.push_back('\t');
            resultData.append(Util::parseFastaHeader(header.getDataByDBKey(maxDbKey, thread_idx)));
            resultData.push_back('\t');
            const TaxonNode* node= t->taxonNode(maxTaxId, false);
            if(node == NULL) {
                resultData.append("Undef");
            }else{
                resultData.append(node->name);
            }
            resultData.push_back('\t');
            resultData.append(SSTR(sequences.getSeqLens(sequences.getId(maxDbKey))));
            resultData.push_back('\n');
            writer.writeData(resultData.c_str(), resultData.size(), outDBKey, thread_idx);
        }
    }
    Debug(Debug::INFO) << "\nDetected potentail conterminetaion in the following Taxons: \n" ;
    Debug(Debug::INFO)  << "Term\tCount\n";
    for(size_t i = 0; i < taxTermCount; i++){
        Debug(Debug::INFO)  << SSTR(i) << "\t" << totalContermCounter[i] << "\n";
    }
    delete t;
    writer.close();
    reader.close();
    header.close();
    sequences.close();
    return EXIT_SUCCESS;
}

