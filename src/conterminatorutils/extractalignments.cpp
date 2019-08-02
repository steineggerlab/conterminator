#include "multipletaxas.h"
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
#include <omptl/omptl_algorithm>


#ifdef OPENMP
#include <omp.h>
#endif

size_t makeIntervalPos(unsigned int key, int pos) {
    struct IntervalPos{
        unsigned int key;
        int pos;
        size_t to_num() const
        {
            size_t num = static_cast<size_t>(key) << 32;
            return num + pos;
        }
    };
    IntervalPos intervalPos;
    intervalPos.key = key;
    intervalPos.pos = pos;
    size_t retVal = intervalPos.to_num();
    return retVal;
}


int extractalignments(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    // bacteria, archaea, eukaryotic, virus
    par.taxonList = "2,2157,2759,10239";
    // unclassified sequences , other sequences,  artifical sequences, retro virus, environmental samples
    par.blacklist = "12908,28384,81077,11632,340016,61964,48479,48510";
    int reportDistinctTaxa = 2;
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

    DBReader<unsigned int> reader(par.db2.c_str(), par.db2Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db3.c_str(), par.db3Index.c_str(), 1, par.compressed, reader.getDbtype());
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

    struct Contamination{
        Contamination(unsigned int key, int start,
                      int end, float seqId,
                      unsigned int len, int range) : key(key), start(start),
                                                     end(end), seqId(seqId),
                                                     len(len), range(range) {}
        Contamination(){};
        unsigned int key;
        int start;
        int end;
        float seqId;
        unsigned int len;
        int range;

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
        std::vector<Contamination> privateContaminations;
        size_t *taxaCounter = new size_t[taxonList.size()];
        IntervalArray speciesRange;
        std::vector<Multipletaxas::TaxonInformation> elements;

        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif

#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            elements.clear();
            memset(taxaCounter, 0, taxonList.size() * sizeof(size_t));
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
            Multipletaxas::assignTaxonomy(elements, data, mapping, t, taxonList, blackList, taxaCounter);
            std::sort(elements.begin(), elements.end(), Multipletaxas::TaxonInformation::compareByTaxAndStart);
            int distinctTaxaCnt = 0;
            //size_t maxTaxCnt = 0;
            //unsigned int maxTaxAncestor = UINT_MAX;
            for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
            }
            // find max. taxa
            for (size_t taxId = 0; taxId < taxonList.size(); taxId++) {
                bool hasTaxa = (taxaCounter[taxId] > 0);
                distinctTaxaCnt += hasTaxa;
            }
            // compute how much the distinct taxon coverage. The maximal taxa is contaminated by the minimal taxa

            unsigned int maxTaxAncestor = UINT_MAX;
            {
                int prevTaxId = -1;
                size_t currTaxCnt = 0;
                size_t maxTaxCnt = 0;
                int prevEnd = 0;

                for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                    if (elements[elementIdx].ancestorTax == 0) {
                        continue;
                    }
                    if (prevTaxId != elements[elementIdx].ancestorTax) {
                        if (currTaxCnt > maxTaxCnt) {
                            maxTaxAncestor = prevTaxId;
                            maxTaxCnt = currTaxCnt;
                        }
                        prevEnd = -1;
                        currTaxCnt = 0;
                    }
                    if (elements[elementIdx].start > prevEnd) {
                        currTaxCnt += (elements[elementIdx].end - elements[elementIdx].start);
                    }
                    prevTaxId = elements[elementIdx].ancestorTax;
                    prevEnd = elements[elementIdx].end;
                }
                if (currTaxCnt > maxTaxCnt) {
                    maxTaxAncestor = prevTaxId;
                }
            }
            // find max. taxa
            // reportDistinctTaxa == 2
            if (distinctTaxaCnt == reportDistinctTaxa) {
                speciesRange.reset();
                if (maxTaxAncestor == UINT_MAX) {
                    Debug(Debug::WARNING) << "Max Tax Id: " << maxTaxAncestor << "!\n";
                    for (size_t i = 0; i < taxonList.size(); i++) {
                        Debug(Debug::WARNING) << taxaCounter[i] << "\n";
                    }
                }

                // fill up interval tree with elements
                for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                    if (static_cast<unsigned int>(elements[elementIdx].ancestorTax) != maxTaxAncestor) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                        speciesRange.insert(res.qStartPos, res.qEndPos);
                    }
                }
                speciesRange.buildRanges();

                int *maxTaxaIdForRange = new int[speciesRange.getRangesSize()];
                memset(maxTaxaIdForRange, 0, speciesRange.getRangesSize() * sizeof(int));
                // count taxons per range
                {
                    int *rangeSizes = new int[speciesRange.getRangesSize() * taxonList.size()];
                    memset(rangeSizes, 0, speciesRange.getRangesSize() * taxonList.size() * sizeof(int));
                    int *maxCntForRange = new int[speciesRange.getRangesSize()];
                    memset(maxCntForRange, 0, speciesRange.getRangesSize() * sizeof(int));

                    std::set <std::pair<int, unsigned int>> rangeEntry;
                    for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                        if (elements[elementIdx].ancestorTax == 0 || elements[elementIdx].start == -1 || elements[elementIdx].end == -1) {
                            continue;
                        }
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                        bool overlapsWithOtherSpecie = speciesRange.doesOverlap(res.qStartPos, res.qEndPos);
                        if (overlapsWithOtherSpecie == true) {
                            int rangeIndex = speciesRange.findIndex(res.qStartPos, res.qEndPos);
                            elements[elementIdx].range = rangeIndex;
                            std::pair<int, unsigned int> rangeQuery = std::make_pair(rangeIndex, res.dbKey);
                            const bool existsAlready = rangeEntry.find(rangeQuery) != rangeEntry.end();
                            if (existsAlready == false) {
                                rangeEntry.insert(rangeQuery);
                                int ancestorIdx = ancestorTax2int[elements[elementIdx].ancestorTax];
                                rangeSizes[rangeIndex * taxonList.size() + ancestorIdx] += 1;
                                if (rangeSizes[rangeIndex * taxonList.size() + ancestorIdx] > maxCntForRange[rangeIndex]) {
                                    maxCntForRange[rangeIndex] = rangeSizes[rangeIndex * taxonList.size() + ancestorIdx];
                                    maxTaxaIdForRange[rangeIndex] = elements[elementIdx].ancestorTax;
                                }
                            }
                        }
                    }
                    delete[] rangeSizes;
                    delete[] maxCntForRange;
                }

                //  pick conterminated
                for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                    if (elements[elementIdx].range != -1) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                        int rangeIndex = speciesRange.findIndex(res.qStartPos, res.qEndPos);
                        // here we write only the contermination in the private ranges
                        if(elements[elementIdx].ancestorTax != maxTaxaIdForRange[rangeIndex]){
                            unsigned int contermKey;
                            int qStartPos = std::min(res.qStartPos, res.qEndPos);
                            int qEndPos = std::max(res.qStartPos, res.qEndPos);
                            int dbStartPos = std::min(res.dbStartPos, res.dbEndPos);
                            int dbEndPos = std::max(res.dbStartPos, res.dbEndPos);
                            int contermStart, contermEnd;
                            unsigned int contermLen;
                            if (t.IsAncestor(maxTaxaIdForRange[rangeIndex], queryTaxon) == false) {
                                contermKey = queryKey;
                                contermStart = qStartPos;
                                contermEnd = qEndPos;
                                contermLen = res.qLen;
                            }else {
                                contermKey = res.dbKey;
                                contermStart = dbStartPos;
                                contermEnd = dbEndPos;
                                contermLen = res.dbLen;
                            }

                            privateContaminations.push_back(Contamination(contermKey, contermStart, contermEnd, res.seqId, contermLen, elements[elementIdx].range ));

                        }
                    }
                }
                delete [] maxTaxaIdForRange;
            }
        }


#pragma omp critical
        {
            allContaminations.insert(allContaminations.end(), privateContaminations.begin(), privateContaminations.end());
        };
        delete[] taxaCounter;
    }

    omptl::sort(allContaminations.begin(), allContaminations.end(), Contamination::compareContaminationByKeyStartEnd);
    size_t writePos = -1;
    unsigned int prevKey = UINT_MAX;

    for(size_t i = 0; i < allContaminations.size(); i++){
        if(prevKey != allContaminations[i].key){
            writePos++;
            allContaminations[writePos]=allContaminations[i];
        }else{
            if(allContaminations[i].start < allContaminations[writePos].end){
                allContaminations[writePos].end = std::max(allContaminations[i].end, allContaminations[writePos].end);
                allContaminations[writePos].seqId = std::max(allContaminations[i].seqId, allContaminations[writePos].seqId);
            }else{
                writePos++;
                allContaminations[writePos]=allContaminations[i];
            }
        }
        prevKey = allContaminations[i].key;
    }
    writePos++;
    allContaminations.resize (writePos);
//    for(size_t i = 0; i < allContaminations.size(); i++){
//        std::cout << allContaminations[i].key << "\t" << allContaminations[i].start << "\t" << allContaminations[i].end << std::endl;
//    }
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
            writer.writeData(buffer, len, 0);
        }
    }

    delete[] ancestorTax2int;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}

