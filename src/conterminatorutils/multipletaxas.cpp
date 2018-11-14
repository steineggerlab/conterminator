#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "filterdb.h"
#include <algorithm>

#ifdef OPENMP
#include <omp.h>
#endif


int multipletaxas(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.taxonList = "2,2157,2759,10239";
    par.parseParameters(argc, argv, command, 3);

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
        Debug(Debug::ERROR) << "names.dmp, nodes.dmp, merged.dmp or delnodes.dmp from NCBI taxdump could not be found!\n";
        EXIT(EXIT_FAILURE);
    }
    std::vector<std::pair<unsigned int, unsigned int>> mapping;
    if(FileUtil::fileExists(std::string(par.db1 + "_mapping").c_str()) == false){
        Debug(Debug::ERROR) << par.db1 + "_mapping" << " does not exist. Please create the taxonomy mapping!\n";
        EXIT(EXIT_FAILURE);
    }
    bool isSorted = Util::readMapping( par.db1 + "_mapping", mapping);
    if(isSorted == false){
        std::stable_sort(mapping.begin(), mapping.end(), ffindexFilter::compareFirstInt());
    }
    std::vector<std::string> ranks = Util::split(par.lcaRanks, ":");

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str());
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads);
    writer.open();

    Debug(Debug::INFO) << "Loading NCBI taxonomy...\n";
    NcbiTaxonomy t(namesFile, nodesFile, mergedFile, delnodesFile);

    Debug(Debug::INFO) << "Add taxonomy information ...\n";

    // a few NCBI taxa are blacklisted by default, they contain unclassified sequences (e.g. metagenomes) or other sequences (e.g. plasmids)
    // if we do not remove those, a lot of sequences would be classified as Root, even though they have a sensible LCA
    std::vector<std::string> list = Util::split(par.taxonList, ",");
    const size_t taxListSize = list.size();
    int* taxalist = new int[taxListSize];
    for (size_t i = 0; i < taxListSize; ++i) {
        taxalist[i] = Util::fast_atoi<int>(list[i].c_str());
    }
    size_t * taxaCounter = new size_t[taxListSize];

#pragma omp parallel
    {
        char *entry[255];
        char buffer[10000];
        unsigned int thread_idx = 0;
        std::string resultData;
        resultData.reserve(4096);
        std::vector< std::pair<unsigned int, char* > > elements;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            Debug::printProgress(i);
            memset(taxaCounter, 0, taxListSize * sizeof(size_t));
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i);
            size_t length = reader.getSeqLens(i);

            if (length == 1) {
                continue;
            }

            std::vector<int> taxa;
            while (*data != '\0') {
                bool isAncestor = false;
                const size_t columns = Util::getWordsOfLine(data, entry, 255);
                if (columns == 0) {
                    Debug(Debug::WARNING) << "Empty entry: " << i << "!";
                    data = Util::skipLine(data);
                    continue;
                }
                unsigned int id = Util::fast_atoi<unsigned int>(entry[0]);
                std::pair<unsigned int, unsigned int> val;
                val.first = id;
                std::vector<std::pair<unsigned int, unsigned int>>::iterator mappingIt = std::upper_bound(
                        mapping.begin(), mapping.end(), val, ffindexFilter::compareToFirstInt);

                if (mappingIt == mapping.end()) {
                    Debug(Debug::WARNING) << "No taxon mapping provided for id " << id << "\n";
                    data = Util::skipLine(data);
                    continue;
                }
                unsigned int taxon = mappingIt->second;
                elements.push_back(std::make_pair(taxon, data));
                // remove blacklisted taxa
                // remove blacklisted taxa
                for (int j = 0; j < taxListSize && isAncestor == false; ++j) {
                    if (taxalist[j] == 0) {
                        isAncestor = std::max(isAncestor, (taxon == 0) ? true : false);
                        continue;
                    }
                    if (taxon == 0) {
                        continue;
                    }
                    std::vector<int> tmpList;

                    bool isTaxaAncestor = t.IsAncestor(taxalist[j], taxon);
                    taxaCounter[j] += isTaxaAncestor;
                    isAncestor = std::max(isAncestor, isTaxaAncestor);
                }
                data = Util::skipLine(data);
            }
            int distinctTaxaCnt = 0;
            for(size_t i = 0; i < taxListSize; i++){
                distinctTaxaCnt += (taxaCounter[i] > 0);
            }

            if(distinctTaxaCnt > 1){
                for(int i = 0 ; i < elements.size(); i++){
                    size_t taxon = elements[i].first;
                    char * data = elements[i].second;
                    char * nextData = Util::skipLine(data);
                    size_t dataSize = nextData - data;
                    resultData.append(data, dataSize-1);
                    resultData.push_back('\t');
                    int len;

                    if(taxon != 0){
                        TaxonNode* node = t.findNode(taxon);
                        std::string lcaRanks = Util::implode(t.AtRanks(node, ranks), ':');
                        if (ranks.empty() == false) {
                            len = snprintf(buffer, 10000, "%d\t%s\t%s\n",
                                           node->taxon, node->rank.c_str(), node->name.c_str());
                        } else {
                            len = snprintf(buffer, 10000, "%d\t%s\t%s\t%s\n",
                                           node->taxon, node->rank.c_str(), node->name.c_str(), lcaRanks.c_str());
                        }
                    } else {
                        len = snprintf(buffer, 1024, "0\tno rank\tunclassified\n");
                    }
                    if(len < 0){
                        Debug(Debug::WARNING) << "Taxon record could not be written. Entry: " << i << "!\n";
                        continue;
                    }
                    resultData.append(buffer, len);
                }
            }
            writer.writeData(resultData.c_str(), resultData.size(), key, thread_idx);
            resultData.clear();
            elements.clear();
        }
    }

    Debug(Debug::INFO) << "\n";
    delete[] taxalist;

    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}
