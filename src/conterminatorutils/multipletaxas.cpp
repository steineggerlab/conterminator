#include "multipletaxas.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "Debug.h"
#include "Util.h"
#include <algorithm>
#include "Matcher.h"
#include "IntervalTree.h"
#include "IntervalArray.h"
#include <set>

#ifdef OPENMP
#include <omp.h>
#endif

int multipletaxas(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
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
    std::vector<std::string> taxlistStr = Util::split(par.taxonList, ",");
    std::vector<int> taxonList;
    int maxTaxa = 0;
    for (size_t i = 0; i < taxlistStr.size(); ++i) {
        taxonList.push_back(Util::fast_atoi<int>(taxlistStr[i].c_str()));
        maxTaxa = std::max(taxonList[i], maxTaxa);
    }
    int* ancestorTax2int = new int[maxTaxa+1];
    for (size_t i = 0; i < taxonList.size(); ++i) {
        ancestorTax2int[taxonList[i]] = i;
    }

    std::vector<std::string> blacklistStr = Util::split(par.blacklist, ",");
    std::vector<int> blackList;
    for (size_t i = 0; i <  blacklistStr.size(); ++i) {
        int currTaxa = Util::fast_atoi<int>(blacklistStr[i].c_str());
        blackList.push_back(currTaxa);
    }



    size_t rangeWithSingleHit=0;
    size_t totalRanges=0;
    Debug::Progress progress(reader.getSize());

#pragma omp parallel reduction (+: rangeWithSingleHit, totalRanges)
    {

        size_t * taxaCounter = new size_t[taxonList.size()];
        char buffer[10000];
        unsigned int thread_idx = 0;
        std::string resultData;
        resultData.reserve(4096);
        IntervalArray speciesRange;
        std::vector<Multipletaxas::TaxonInformation> elements;

#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            resultData.clear();
            elements.clear();
            speciesRange.reset();
            memset(taxaCounter, 0, taxonList.size() * sizeof(size_t));
            unsigned int key = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            size_t length = reader.getSeqLens(i);

            if (length == 1) {
                continue;
            }
            // find taxonomical information
            std::vector<int> taxa;
            Multipletaxas::assignTaxonomy(elements, data, mapping, t, taxonList, blackList, taxaCounter);

            int distinctTaxaCnt = 0;
            size_t maxTaxCnt = 0;
            int maxTaxId = -1;
            // find max. taxa
            for(size_t taxId = 0; taxId < taxonList.size(); taxId++){
                bool hasTaxa = (taxaCounter[taxId] > 0);
                distinctTaxaCnt += hasTaxa;
                if(hasTaxa && taxaCounter[taxId] > maxTaxCnt){
                    maxTaxCnt = taxaCounter[taxId];
                    maxTaxId =  taxonList[taxId];
                }
            }
            if(distinctTaxaCnt > 1) {
                if (maxTaxId == -1) {
                    Debug(Debug::WARNING) << "Max Tax Id: " << maxTaxId << "!\n";
                    for (size_t i = 0; i < taxonList.size(); i++) {
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
                int * rangSizes = new int[speciesRange.getRangesSize()*taxonList.size()];
                memset(rangSizes, 0, speciesRange.getRangesSize() * taxonList.size() * sizeof(int));
                int * simpleRanges = new int[speciesRange.getRangesSize()];
                memset(simpleRanges, 0, speciesRange.getRangesSize() * sizeof(int));

                {
                    std::set<Multipletaxas::RangEntry, Multipletaxas::RangEntry> rangeEntry;
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
                            Multipletaxas::RangEntry rangeQuery(rangeIndex, elements[elementIdx].ancestorTax, taxon, res.dbKey);
                            const bool existsAlready = rangeEntry.find(rangeQuery) != rangeEntry.end();
                            if(existsAlready == false) {
                                rangeEntry.insert(rangeQuery);
                                rangSizes[rangeIndex * taxonList.size() + ancestorTax2int[elements[elementIdx].ancestorTax]] += 1;
                            }
                        } else {
                            continue;
                        }
                    }
                }
                std::sort(elements.begin(), elements.end(), Multipletaxas::TaxonInformation::compareByRange);
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
                            for(size_t taxIdx = 0; taxIdx < taxonList.size(); taxIdx++){
                                int size = rangSizes[elements[elementIdx].range*taxonList.size() + taxIdx];
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
                        for(size_t taxIdx = 0; taxIdx < taxonList.size(); taxIdx++){
                            int size = rangSizes[elements[elementIdx].range*taxonList.size() + taxIdx];
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
    delete[] ancestorTax2int;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}