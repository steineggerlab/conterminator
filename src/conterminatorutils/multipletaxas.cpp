#include "NcbiTaxonomy.h"
#include "Parameters.h"
#include "DBWriter.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include "filterdb.h"
#include <algorithm>
#include "Matcher.h"
#include "IntervalTree.h"
#include "IntervalArray.h"

#ifdef OPENMP
#include <omp.h>
#endif



int multipletaxas(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    // bacteria, archaea, eukaryotic, virus
    par.taxonList = "2,2157,2759,10239";
    // unclassified sequences , other sequences,  artifical sequences, retro virus
    par.blacklist = "12908,28384,81077,35268";
    par.blacklist = "12908,28384,81077,35268, 340016, 61964, 48479, 48510";

    par.parseParameters(argc, argv, command, 3);

//    int ints[][2] = {{4,6},
//                     {15, 20},
//                     {10, 25},
//                     {17, 19},
//                     {5, 6},
//                     {12, 15},
//                     {12, 16},
//                     {12, 17},
//                     {12, 18},
//                     {12, 19},
//                     {30, 39}};
//
//    IntervalArray tree;
//    tree.reset();
//    for (int i = 0; i < 2; i++)
//        tree.insert(ints[i][0], ints[i][1]);
//    tree.print();
//    std::cout << tree.doesOverlap(10,20) << std::endl;
//    std::cout << tree.doesOverlap(34,28) << std::endl;
//
//
//
//    //int ints[][2] =  {{15,20} , {4,25}, {3,30}};
//    tree.reset();
//    int n = sizeof(ints)/sizeof(ints[0]);
//    for (int i = 0; i < n; i++)
//        tree.insert(ints[i][0], ints[i][1]);
//    tree.print();
//    std::cout << tree.doesOverlap(26,27) << std::endl;
//    std::cout << tree.doesOverlap(1,2) << std::endl;
//    std::cout << tree.doesOverlap(15,20) << std::endl;
//    std::cout << tree.doesOverlap(14,20) << std::endl;
//    std::cout << tree.doesOverlap(14,19) << std::endl;
//    std::cout << tree.doesOverlap(34,28) << std::endl;
//    std::cout << "New tree" << std::endl;

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
    std::vector<std::string> taxlist = Util::split(par.taxonList, ",");
    const size_t taxListSize = taxlist.size();
    int* taxalist = new int[taxListSize];
    for (size_t i = 0; i < taxListSize; ++i) {
        taxalist[i] = Util::fast_atoi<int>(taxlist[i].c_str());
    }

    std::vector<std::string> blacklist = Util::split(par.blacklist, ",");
    const size_t blackListSize = blacklist.size();
    int* blackList = new int[blackListSize];
    for (size_t i = 0; i < blackListSize; ++i) {
        int currTaxa = Util::fast_atoi<int>(blacklist[i].c_str());
        blackList[i] = currTaxa;
    }

    struct TaxonInformation{
        TaxonInformation(int currTaxa, int ancestorTax, char * data) :
            currTaxa(currTaxa), ancestorTax(ancestorTax), data(data){}
        int currTaxa;
        int ancestorTax;
        char * data;
    };
#pragma omp parallel
    {
        size_t * taxaCounter = new size_t[taxListSize];
        char *entry[255];
        char buffer[10000];
        unsigned int thread_idx = 0;
        std::string resultData;
        resultData.reserve(4096);
        std::vector<TaxonInformation> elements;
        IntervalArray speciesRange;

#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            Debug::printProgress(i);
            resultData.clear();
            elements.clear();
            speciesRange.reset();
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

                if (mappingIt->first != val.first) {
//                    Debug(Debug::WARNING) << "No taxon mapping provided for id " << id << "\n";
                    data = Util::skipLine(data);
                    continue;
                }
                unsigned int taxon = mappingIt->second;
                if (taxon == 0) {
                    goto next;
                }
                // remove blacklisted taxa
                for (size_t j = 0; j < blackListSize; ++j) {
                    if (t.IsAncestor(blackList[j], taxon)) {
                        goto next;
                    }
                }

                size_t j;
                for ( j = 0; j < taxListSize ; ++j) {
                    bool isTaxaAncestor = t.IsAncestor(taxalist[j], taxon);
                    taxaCounter[j] += isTaxaAncestor;
                    isAncestor = std::max(isAncestor, isTaxaAncestor);
                    if(isAncestor==true){
                        break;
                    }
                }
                if(isAncestor){
                    elements.push_back(TaxonInformation(taxon, taxalist[j], data));
                }else{
                    elements.push_back(TaxonInformation(taxon, 0, data));
                }
                next:
                data = Util::skipLine(data);
            }
            int distinctTaxaCnt = 0;
            size_t maxTaxCnt = 0;
            int maxTaxId = -1;

            for(size_t i = 0; i < taxListSize; i++){
                bool hasTaxa = (taxaCounter[i] > 0);
                distinctTaxaCnt += hasTaxa;
                if(hasTaxa && taxaCounter[i] > maxTaxCnt){
                    maxTaxCnt = taxaCounter[i];
                    maxTaxId =  taxalist[i];
                }
            }

            if(distinctTaxaCnt > 1){
                if(maxTaxId == -1){
                    Debug(Debug::WARNING) << "Max Tax Id: " << maxTaxId << "!\n";
                    for(size_t i = 0; i < taxListSize; i++) {
                        Debug(Debug::WARNING) << taxaCounter[i] << "\n";
                    }
                }

                // fill up interval tree with elements
                for(size_t i = 0 ; i < elements.size(); i++) {
                    if(elements[i].ancestorTax != maxTaxId) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[i].data, true);
                        speciesRange.insert(res.qStartPos, res.qEndPos);
                    }
                }
                for(size_t i = 0 ; i < elements.size(); i++){
                    size_t taxon = elements[i].currTaxa;
                    Matcher::result_t res = Matcher::parseAlignmentRecord(elements[i].data, true);
                    char * data = elements[i].data;
                    bool overlapsWithOtherSpecie = speciesRange.doesOverlap(res.qStartPos, res.qEndPos);
                    if(taxon != 0 && overlapsWithOtherSpecie == true){
                        char * nextData = Util::skipLine(data);
                        size_t dataSize = nextData - data;
                        resultData.append(data, dataSize-1);
                        resultData.push_back('\t');
                        int len;

                        TaxonNode* node = t.findNode(taxon);
                        if(node == NULL){
                            len = snprintf(buffer, 1024, "0\tno rank\tunclassified\n");
                        }else{
                            std::string lcaRanks = Util::implode(t.AtRanks(node, ranks), ':');
                            if (ranks.empty() == false) {
                                len = snprintf(buffer, 10000, "%d\t%d\t%s\t%s\n",
                                               elements[i].ancestorTax, node->taxon, node->rank.c_str(), node->name.c_str());
                            } else {
                                len = snprintf(buffer, 10000, "%d\t%d\t%s\t%s\t%s\n",
                                               elements[i].ancestorTax, node->taxon, node->rank.c_str(), node->name.c_str(), lcaRanks.c_str());
                            }
                        }
                        if(len < 0){
                            Debug(Debug::WARNING) << "Taxon record could not be written. Entry: " << i << "!\n";
                            continue;
                        }
                        resultData.append(buffer, len);
                    } else {
                        continue;
                    }

                }
            }
            writer.writeData(resultData.c_str(), resultData.size(), key, thread_idx);
        }
        delete [] taxaCounter;
    }

    Debug(Debug::INFO) << "\n";
    delete[] taxalist;
    delete[] blackList;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}
