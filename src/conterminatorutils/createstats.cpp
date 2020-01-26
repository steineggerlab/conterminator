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
#include "LocalParameters.h"

#ifdef OPENMP
#include <omp.h>
#endif


int createstats(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    // bacteria, archaea, eukaryotic, virus
    par.parseParameters(argc, argv, command, true, 0, 0);

    NcbiTaxonomy * t = NcbiTaxonomy::openTaxonomy(par.db1);

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
#pragma omp parallel
    {
        KingdomExpression kingdomExpression(par.kingdoms, *t);
        size_t taxTermCount = kingdomExpression.getTaxTerms().size();
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
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getEntryLen(i);

            if (length == 1) {
                continue;
            }
            // find taxonomical information
            memset(taxaCounter, 0, taxTermCount * sizeof(size_t));
            TaxonUtils::assignTaxonomy(elements, data, mapping, *t, kingdomExpression, blackList, taxaCounter, true);
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
            int minTaxTerm = -1;
            int maxTaxId = -1;
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
            if(minTaxTerm == -1 || maxTaxId == -1){
                continue;
            }
            struct Contamination{
                unsigned int key;
                int taxId;
                unsigned int start;
                unsigned int end;
                bool operator()(const Contamination& first, const Contamination& second) const
                {
                    if(first.key < second.key )
                        return true;
                    if(second.key < first.key )
                        return false;
                    if(first.taxId < second.taxId )
                        return true;
                    if(second.taxId < first.taxId )
                        return false;
                    return false;
                }
            };
            std::set<Contamination, Contamination> minDbKeys;

            struct TaxonCnt{
                unsigned int key;
                unsigned int cnt;
            };
            std::map<int, TaxonCnt> maxTaxon;
            for(size_t i = 0; i < elements.size(); i++){
                if(elements[i].termId == minTaxTerm){
                    Contamination conterm;
                    conterm.key = elements[i].dbKey;
                    conterm.start = elements[i].start;
                    conterm.end = elements[i].end;
                    conterm.taxId = elements[i].currTaxa;
                    minDbKeys.insert(conterm);
                }
                if(elements[i].termId == maxTaxId) {
                    maxTaxon[elements[i].currTaxa].cnt++;
                    maxTaxon[elements[i].currTaxa].key = elements[i].dbKey;
                }
            }
            std::map<int, TaxonCnt>::iterator maxTaxonIt;
            unsigned int maxDbKey=UINT_MAX;
            maxTaxCnt = 0;
            for ( maxTaxonIt = maxTaxon.begin(); maxTaxonIt != maxTaxon.end(); maxTaxonIt++ )
            {
                if(maxTaxonIt->second.cnt >= maxTaxCnt){
                    maxTaxId = maxTaxonIt->first;
                    maxDbKey =  maxTaxonIt->second.key;
                    maxTaxCnt = maxTaxonIt->second.cnt;
                }
            }

            std::set<Contamination>::iterator minDbKeysIt = minDbKeys.begin();
            std::string taxons;
            std::string dbLength;
            std::string contermStartPos;
            std::string contermEndPos;

            unsigned int outDBKey = UINT_MAX;
            for (;minDbKeysIt != minDbKeys.end(); minDbKeysIt++){
                unsigned int dbKey = minDbKeysIt->key;
                outDBKey = std::min(dbKey, outDBKey);
                int taxId = minDbKeysIt->taxId;
                resultData.append(Util::parseFastaHeader(header.getDataByDBKey(dbKey, thread_idx)));
                dbLength.append(SSTR(sequences.getSeqLen(sequences.getId(dbKey))));
                contermStartPos.append(SSTR(minDbKeysIt->start));
                contermEndPos.append(SSTR(minDbKeysIt->end));
                const TaxonNode* node= t->taxonNode(taxId, false);
                if(node == NULL){
                    taxons.append("Undef");
                }else{
                    taxons.append(node->name);
                }
                taxons.push_back(',');
                dbLength.push_back(',');
                contermStartPos.push_back(',');
                contermEndPos.push_back(',');
                resultData.push_back(',');
            }
            taxons.pop_back();
            dbLength.pop_back();
            contermStartPos.pop_back();
            contermEndPos.pop_back();
            resultData.pop_back();
            resultData.push_back('\t');
            resultData.append(taxons);
            resultData.push_back('\t');
            resultData.append(contermStartPos);
            resultData.push_back('\t');
            resultData.append(contermEndPos);
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
            resultData.append(SSTR(sequences.getSeqLen(sequences.getId(maxDbKey))));
            resultData.push_back('\n');
            writer.writeData(resultData.c_str(), resultData.size(), queryKey, thread_idx);
        }
        delete [] taxaCounter;
    }
    Debug(Debug::INFO) << "\nDetected potentail conterminetaion in the following Taxons: \n" ;
    Debug(Debug::INFO)  << "Term\tCount\n";

    delete t;
    writer.close();
    reader.close();
    header.close();
    sequences.close();
    return EXIT_SUCCESS;
}

