#include "multipletaxas.h"
#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include <omptl/omptl_algorithm>
#include <set>


#ifdef OPENMP
#include <omp.h>
#endif


int createstats(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    // bacteria, archaea, eukaryotic, virus
    par.taxonList = "2,2157,2759,10239";
    // unclassified sequences , other sequences,  artifical sequences, retro virus, environmental samples
    par.blacklist = "12908,28384,81077,11632,340016,61964,48479,48510";
    par.parseParameters(argc, argv, command, true, 0, 0);

    std::string nodesFile = par.db1 + "_nodes.dmp";
    std::string namesFile = par.db1 + "_names.dmp";
    std::string mergedFile = par.db1 + "_merged.dmp";
    std::string delnodesFile = par.db1 + "_delnodes.dmp";
    if (FileUtil::fileExists(nodesFile.c_str())
        && FileUtil::fileExists(namesFile.c_str())
        && FileUtil::fileExists(mergedFile.c_str())
        && FileUtil::fileExists(delnodesFile.c_str())) {
    } else if (FileUtil::fileExists("nodes.dmp")
               && FileUtil::fileExists("names.dmp")
               && FileUtil::fileExists("merged.dmp")
               && FileUtil::fileExists("delnodes.dmp")) {
        nodesFile = "nodes.dmp";
        namesFile = "names.dmp";
        mergedFile = "merged.dmp";
        delnodesFile = "delnodes.dmp";
    } else {
        Debug(Debug::ERROR)
                << "names.dmp, nodes.dmp, merged.dmp or delnodes.dmp from NCBI taxdump could not be found!\n";
        EXIT(EXIT_FAILURE);
    }
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

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    Debug(Debug::INFO) << "Loading NCBI taxonomy...\n";
    NcbiTaxonomy t(namesFile, nodesFile, mergedFile);

    // a few NCBI taxa are blacklisted by default, they contain unclassified sequences (e.g. metagenomes) or other sequences (e.g. plasmids)
    // if we do not remove those, a lot of sequences would be classified as Root, even though they have a sensible LCA
    std::vector<std::string> taxlistStr = Util::split(par.taxonList, ",");
    std::vector<int> taxonList;
    int maxTaxa = 0;
    for (size_t i = 0; i < taxlistStr.size(); ++i) {
        taxonList.push_back(Util::fast_atoi<int>(taxlistStr[i].c_str()));
        maxTaxa = std::max(taxonList[i], maxTaxa);
    }
    int *ancestorTax2int = new int[maxTaxa + 1];
    for (size_t i = 0; i < taxonList.size(); ++i) {
        ancestorTax2int[taxonList[i]] = i;
    }

    std::vector<std::string> blacklistStr = Util::split(par.blacklist, ",");
    std::vector<int> blackList;
    for (size_t i = 0; i < blacklistStr.size(); ++i) {
        int currTaxa = Util::fast_atoi<int>(blacklistStr[i].c_str());
        blackList.push_back(currTaxa);
    }


    Debug::Progress progress(reader.getSize());
    size_t *totalContermCounter = new size_t[taxonList.size()];
    memset(totalContermCounter, 0, taxonList.size() * sizeof(size_t));
#pragma omp parallel
    {
        std::string resultData;
        resultData.reserve(4096);
        size_t *taxaCounter = new size_t[taxonList.size()];
        unsigned int thread_idx = 0;
        std::vector<Multipletaxas::TaxonInformation> elements;
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

            unsigned int queryTaxon = Multipletaxas::getTaxon(queryKey, mapping);
            if (queryTaxon == 0 || queryTaxon == UINT_MAX) {
                continue;
            }
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getSeqLens(i);

            if (length == 1) {
                continue;
            }
            // find taxonomical information
            memset(taxaCounter, 0, taxonList.size() * sizeof(size_t));
            Multipletaxas::assignTaxonomy(elements, data, mapping, t, taxonList, blackList, taxaCounter);
            // recount
            memset(taxaCounter, 0, taxonList.size() * sizeof(size_t));
            for(size_t i = 0; i < elements.size(); i++){
                char dbKeyBuffer[255 + 1];
                Util::parseKey(elements[i].data, dbKeyBuffer);
                const unsigned int key = (unsigned int) strtoul(dbKeyBuffer, NULL, 10);
                const bool existsAlready = idDetected.find(key) != idDetected.end();
                if(existsAlready == false) {
                    taxaCounter[ancestorTax2int[elements[i].ancestorTax]]++;
                    idDetected.insert(key);
                }
            }
            std::sort(elements.begin(), elements.end(), Multipletaxas::TaxonInformation::compareByTaxAndStart);
            // find max. taxa
            size_t minTaxCnt = SIZE_T_MAX;
            int minTaxId = 0;
            int maxTaxId = 0;
            size_t maxTaxCnt = 0;
            for (size_t taxId = 0; taxId < taxonList.size(); taxId++) {
                resultData.append(SSTR(taxaCounter[taxId]));
                if(taxaCounter[taxId] < minTaxCnt && taxaCounter[taxId] != 0){
                    minTaxCnt = taxaCounter[taxId];
                    minTaxId = taxonList[taxId];
                }

                if(taxaCounter[taxId] > maxTaxCnt){
                    maxTaxCnt = taxaCounter[taxId];
                    maxTaxId = taxonList[taxId];
                }
                resultData.push_back('\t');
            }
            __sync_fetch_and_add(&(totalContermCounter[ancestorTax2int[minTaxId]]), 1);

            std::set<std::pair<unsigned int, int> > minDbKeys;
            std::map<int, size_t> maxTaxon;

            for(size_t i = 0; i < elements.size(); i++){
                if(elements[i].ancestorTax == minTaxId){
                    std::pair<unsigned int, int> pair =  std::make_pair(elements[i].dbKey, elements[i].currTaxa);
                    minDbKeys.insert(pair);
                }
                if(elements[i].ancestorTax == maxTaxId) {
                    maxTaxon[elements[i].currTaxa]++;
                }
            }
            std::map<int, size_t>::iterator maxTaxonIt;
            maxTaxCnt = 0;
            for ( maxTaxonIt = maxTaxon.begin(); maxTaxonIt != maxTaxon.end(); maxTaxonIt++ )
            {
                if(maxTaxonIt->second > maxTaxCnt){
                    maxTaxId = maxTaxonIt->first;
                    maxTaxCnt = maxTaxonIt->second;
                }
            }

            std::set<std::pair<unsigned int, int>>::iterator minDbKeysIt = minDbKeys.begin();
            std::string taxons;
            for (;minDbKeysIt != minDbKeys.end(); minDbKeysIt++){
                unsigned int dbKey = minDbKeysIt->first;
                int taxId = minDbKeysIt->second;
                resultData.append(Util::parseFastaHeader(header.getDataByDBKey(dbKey, thread_idx)));
                const TaxonNode* node= t.taxonNode(taxId, false);
                if(node == NULL){
                    taxons.append("Undef");
                }else{
                    taxons.append(node->name);
                }
                taxons.push_back(',');
                resultData.push_back(',');
            }
            taxons.pop_back();
            resultData.pop_back();
            resultData.push_back('\t');
            resultData.append(taxons);
            resultData.push_back('\t');
            const TaxonNode* node= t.taxonNode(maxTaxId, false);
            if(node == NULL) {
                resultData.append("Undef");
            }else{
                resultData.append(node->name);
            }
            resultData.push_back('\n');
            writer.writeData(resultData.c_str(), resultData.size(), queryKey, thread_idx);
        }
    }
    Debug(Debug::INFO) << "\nDetected potentail conterminetaion in the following Taxons: \n" ;
    for(size_t i = 0; i < taxonList.size(); i++){
        const TaxonNode* node= t.taxonNode(taxonList[i], false);
        if(node == NULL) {
            Debug(Debug::INFO) << "Undef";
        }else{
            Debug(Debug::INFO) << node->name;
        }
        Debug(Debug::INFO) << "\t" << totalContermCounter[i] << "\n";
    }
    delete[] ancestorTax2int;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}

