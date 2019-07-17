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
#include <set>

#ifdef OPENMP
#include <omp.h>
#endif

struct RangEntry {
    unsigned int range;
    unsigned int kingdom;
    unsigned int species;
    unsigned int id;

    RangEntry() {};
    RangEntry(unsigned int range, unsigned int kingdom, unsigned int species, unsigned int id) :
            range(range), kingdom(kingdom), species(species), id(id) {}

    bool operator()(const RangEntry &lhs, const RangEntry &rhs) {
        if (lhs.range < rhs.range) return true;
        if (rhs.range < lhs.range) return false;
        if (lhs.kingdom < rhs.kingdom) return true;
        if (rhs.kingdom < lhs.kingdom) return false;
        if (lhs.species < rhs.species) return true;
        if (rhs.species < lhs.species) return false;
        if (lhs.id < rhs.id) return true;
        return false;
    }
};


int multipletaxas(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    // bacteria, archaea, eukaryotic, virus
    par.taxonList = "2,2157,2759,10239";
    // unclassified sequences , other sequences,  artifical sequences, retro virus, environmental samples
    par.blacklist = "12908,28384,81077,11632,340016,61964,48479,48510";

    par.parseParameters(argc, argv, command, true, 0, 0);


//    int ints[][2] = {{1,2},
//                     {4,6},
//                     {15, 20},
//                     {22, 25},
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
//    for (int i = 0; i < 4; i++)
//        tree.insert(ints[i][0], ints[i][1]);
//    tree.print();
//    std::cout << tree.doesOverlap(10,20) << std::endl;
//    std::cout << tree.doesOverlap(34,28) << std::endl;
//    tree.buildRanges();
//    std::cout << tree.findIndex(3,5) << std::endl;
//    std::cout << tree.findIndex(5,5) << std::endl;
//    std::cout << tree.findIndex(23,24) << std::endl;
//
//    std::cout << tree.findIndex(10,20) << std::endl;
//    std::cout << tree.findIndex(16,19) << std::endl;
//    std::cout << tree.findIndex(17,21) << std::endl;
//    std::cout << tree.findIndex(7,8) << std::endl;

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

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    Debug(Debug::INFO) << "Loading NCBI taxonomy...\n";
    NcbiTaxonomy t(namesFile, nodesFile, mergedFile);

    Debug(Debug::INFO) << "Add taxonomy information ...\n";

    // a few NCBI taxa are blacklisted by default, they contain unclassified sequences (e.g. metagenomes) or other sequences (e.g. plasmids)
    // if we do not remove those, a lot of sequences would be classified as Root, even though they have a sensible LCA
    std::vector<std::string> taxlist = Util::split(par.taxonList, ",");
    const size_t taxListSize = taxlist.size();
    int* taxalist = new int[taxListSize];
    int maxTaxa = 0;
    for (size_t i = 0; i < taxListSize; ++i) {
        taxalist[i] = Util::fast_atoi<int>(taxlist[i].c_str());
        maxTaxa = std::max(taxalist[i], maxTaxa);
    }
    int* ancestorTax2int = new int[maxTaxa+1];
    for (size_t i = 0; i < taxListSize; ++i) {
        ancestorTax2int[taxalist[i]] = i;
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
                currTaxa(currTaxa), ancestorTax(ancestorTax), data(data), overlaps(false), range(-1){}
        int currTaxa;
        int ancestorTax;
        char * data;
        bool overlaps;
        int range;
        static bool compareByRange(const TaxonInformation &first, const TaxonInformation &second) {
            if (first.range < second.range)
                return true;
            if (second.range < first.range)
                return false;
            if(first.ancestorTax < second.ancestorTax )
                return true;
            if(second.ancestorTax < first.ancestorTax )
                return false;
            if(first.currTaxa < second.currTaxa )
                return true;
            if(second.currTaxa < first.currTaxa )
                return false;
            return false;
        }
    };

    size_t rangeWithSingleHit=0;
    size_t totalRanges=0;
    Debug::Progress progress(reader.getSize());

#pragma omp parallel reduction (+: rangeWithSingleHit, totalRanges)
    {
        size_t * taxaCounter = new size_t[taxListSize];
        const char *entry[255];
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
            progress.updateProgress();
            resultData.clear();
            elements.clear();
            speciesRange.reset();
            memset(taxaCounter, 0, taxListSize * sizeof(size_t));
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getSeqLens(i);

            if (length == 1) {
                continue;
            }
            // find taxonomical information
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
            // find max. taxa
            for(size_t taxId = 0; taxId < taxListSize; taxId++){
                bool hasTaxa = (taxaCounter[taxId] > 0);
                distinctTaxaCnt += hasTaxa;
                if(hasTaxa && taxaCounter[taxId] > maxTaxCnt){
                    maxTaxCnt = taxaCounter[taxId];
                    maxTaxId =  taxalist[taxId];
                }
            }
            // find max. taxa
            if(distinctTaxaCnt > 1) {
                if (maxTaxId == -1) {
                    Debug(Debug::WARNING) << "Max Tax Id: " << maxTaxId << "!\n";
                    for (size_t i = 0; i < taxListSize; i++) {
                        Debug(Debug::WARNING) << taxaCounter[i] << "\n";
                    }
                }

                // fill up interval tree with elements
                for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                    if (elements[elementIdx].ancestorTax != maxTaxId) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                        speciesRange.insert(res.qStartPos, res.qEndPos);
                    }
                }
                speciesRange.buildRanges();
                int * rangSizes = new int[speciesRange.getRangesSize()*taxListSize];
                memset(rangSizes, 0, speciesRange.getRangesSize() * taxListSize * sizeof(int));
                int * simpleRanges = new int[speciesRange.getRangesSize()];
                memset(simpleRanges, 0, speciesRange.getRangesSize() * sizeof(int));

                {
                    std::set<RangEntry, RangEntry> rangeEntry;
    //                speciesRange.print();
                    // fill range taxon array
                    for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                        size_t taxon = elements[elementIdx].currTaxa;
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                        bool overlapsWithOtherSpecie = speciesRange.doesOverlap(res.qStartPos, res.qEndPos);
                        if (taxon != 0 && overlapsWithOtherSpecie == true) {
                            elements[elementIdx].overlaps = true;
                            int rangeIndex = speciesRange.findIndex(res.qStartPos, res.qEndPos);
                            elements[elementIdx].range = rangeIndex;
                            RangEntry rangeQuery(rangeIndex, elements[elementIdx].ancestorTax, taxon, res.dbKey);
                            const bool existsAlready = rangeEntry.find(rangeQuery) != rangeEntry.end();
                            if(existsAlready == false) {
                                rangeEntry.insert(rangeQuery);
                                rangSizes[rangeIndex * taxListSize + ancestorTax2int[elements[elementIdx].ancestorTax]] += 1;
                            }
                        } else {
                            continue;
                        }
                    }
                }
                std::sort(elements.begin(), elements.end(), TaxonInformation::compareByRange);
                {
                    bool hasSizeOne = false;
                    bool firstRange = true;
                    int taxaCount = 0;
                    int totalCount = 0;
                    int range = -1;
                    for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                        if (elements[elementIdx].overlaps) {
                            if (elements[elementIdx].range != range) {
                                if (taxaCount == 2 && totalCount >= 3 && hasSizeOne == true) {
                                    rangeWithSingleHit++;
                                    simpleRanges[range] = 1;
                                }
                                taxaCount = 0;
                                range = elements[elementIdx].range;
                                hasSizeOne = false;
                                firstRange = true;
                                totalCount = 0;
                                totalRanges++;
                            }
                            for(size_t taxIdx = 0; taxIdx < taxListSize; taxIdx++){
                                int size = rangSizes[elements[elementIdx].range*taxListSize + taxIdx];
                                hasSizeOne = std::max(hasSizeOne, (firstRange && size == 1));
                                taxaCount += (firstRange && size > 0);
                                totalCount += (firstRange) ? size : 0;

                            }
                            firstRange = false;
                        }
                    }
                    if(taxaCount == 2 && totalCount >= 3 && hasSizeOne == true){
                        simpleRanges[range] = 1;
                        rangeWithSingleHit++;
                        totalRanges++;
                    }
                }

                // print stuff
                for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                    Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                    size_t taxon = elements[elementIdx].currTaxa;
                    if(elements[elementIdx].overlaps){
                        char *data = elements[elementIdx].data;
                        char *nextData = Util::skipLine(data);
                        size_t dataSize = nextData - data;
                        resultData.append(data, dataSize - 1);
                        resultData.push_back('\t');
                        resultData.append(SSTR(elements[elementIdx].range));
                        resultData.push_back('\t');
                        IntervalArray::Range range = speciesRange.getRange(elements[elementIdx].range);
                        resultData.append(SSTR(range.start));
                        resultData.push_back('\t');
                        resultData.append(SSTR(range.end));
                        resultData.push_back('\t');
                        resultData.append(SSTR(simpleRanges[elements[elementIdx].range]));
                        resultData.push_back('\t');

                        // print ranges
                        for(size_t taxIdx = 0; taxIdx < taxListSize; taxIdx++){
                            int size = rangSizes[elements[elementIdx].range*taxListSize + taxIdx];
                            resultData.append(SSTR(size));
                            resultData.push_back('\t');
                        }
                        int len;
                        const TaxonNode *node = t.taxonNode(taxon);
                        if (node == NULL) {
                            len = snprintf(buffer, 1024, "0\tno rank\tunclassified\n");
                        } else {
                            std::string lcaRanks = Util::implode(t.AtRanks(node, ranks), ':');
                            len = snprintf(buffer, 10000, "%d\t%d\t%s\t%s\n",
                                           elements[elementIdx].ancestorTax, node->taxId, node->rank.c_str(),
                                           node->name.c_str());

                        }
                        if (len < 0) {
                            Debug(Debug::WARNING) << "Taxon record could not be written. Entry: " << elementIdx << "!\n";
                            continue;
                        }
                        resultData.append(buffer, len);
                    }
                }
                delete [] rangSizes;
                delete [] simpleRanges;
            }
            writer.writeData(resultData.c_str(), resultData.size(), key, thread_idx);
        }
        delete [] taxaCounter;
    }

    Debug(Debug::INFO) << "\n";
    Debug(Debug::INFO) << "Potential conterminated ranges: " << totalRanges << "\n";
    Debug(Debug::INFO) << "Highly likely conterminated ranges: " << rangeWithSingleHit << "\n";
    delete[] taxalist;
    delete[] blackList;
    delete[] ancestorTax2int;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}
