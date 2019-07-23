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
    

    IntervalTree<size_t, unsigned int>::interval_vector allRanges;
    Debug::Progress progress(reader.getSize());

#pragma omp parallel reduction (+: rangeWithSingleHit, totalRanges)
    {
        IntervalTree<size_t, unsigned int>::interval_vector privateRanges;

        size_t *taxaCounter = new size_t[taxonList.size()];
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
            std::vector<int> taxa;
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
            int prevTaxId = -1;
            size_t currTaxCnt = 0;
            size_t maxTaxCnt = 0;
            int prevEnd = 0;
            unsigned int maxTaxAncestor = UINT_MAX;

            for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                if(elements[elementIdx].ancestorTax == 0){
                    continue;
                }
                if(prevTaxId != elements[elementIdx].ancestorTax){
                    if(currTaxCnt > maxTaxCnt){
                        maxTaxAncestor = prevTaxId;
                        maxTaxCnt = currTaxCnt;
                    }
                    prevEnd = -1;
                    currTaxCnt = 0;
                }
                if(elements[elementIdx].start > prevEnd) {
                    currTaxCnt += (elements[elementIdx].end - elements[elementIdx].start);
                }
                prevTaxId = elements[elementIdx].ancestorTax;
                prevEnd = elements[elementIdx].end;
            }
            if(currTaxCnt > maxTaxCnt){
                maxTaxAncestor = prevTaxId;
            }

            // find max. taxa
            // reportDistinctTaxa == 2
            if (distinctTaxaCnt == reportDistinctTaxa) {
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

                // pick conterminated ranges
                for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                    if (static_cast<unsigned int>(elements[elementIdx].ancestorTax) != maxTaxAncestor) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                        bool overlapsWithOtherSpecie = speciesRange.doesOverlap(res.qStartPos, res.qEndPos);
                        if(overlapsWithOtherSpecie == true){
                            if (t.IsAncestor(maxTaxAncestor, queryTaxon)) {
                                int qStartPos = std::min(res.qStartPos, res.qEndPos);
                                int qEndPos = std::max(res.qStartPos, res.qEndPos);
                                privateRanges.push_back(Interval<size_t, unsigned int>(makeIntervalPos(queryKey, qStartPos),
                                                                                       makeIntervalPos(queryKey, qEndPos),
                                                                                       res.dbKey));
                            } else {
                                int dbStartPos = std::min(res.dbStartPos, res.dbEndPos);
                                int dbEndPos = std::max(res.dbStartPos, res.dbEndPos);
                                privateRanges.push_back(
                                        Interval<size_t, unsigned int>(makeIntervalPos(res.dbKey, dbStartPos),
                                                                       makeIntervalPos(res.dbKey, dbEndPos),
                                                                       queryKey));
                            }
                        }
                    }
                }


            }
        }


#pragma omp critical
        {
            allRanges.insert(allRanges.end(), privateRanges.begin(), privateRanges.end());
        }
        delete[] taxaCounter;
    }

    Debug::Progress progress2(reader.getSize());

    IntervalTree<size_t, unsigned int> tree(std::move(allRanges));

    struct SmallAlignment{
        unsigned int conterminatedKey;
        unsigned int dbKey;
        int score;
        float seqId;
        float eval;
        int qStartPos;
        int qEndPos;
        unsigned int qLen;
        int dbStartPos;
        int dbEndPos;
        unsigned int dbLen;
        char strand;
        SmallAlignment(unsigned int conterminatedKey, unsigned int dbKey, int score, float seqId, float eval,
                       int qStartPos, int qEndPos, unsigned int qLen, int dbStartPos, int dbEndPos,
                       unsigned int dbLen, char strand) : conterminatedKey(conterminatedKey), dbKey(dbKey),
                                                          score(score), seqId(seqId), eval(eval), qStartPos(qStartPos), qEndPos(qEndPos),
                                                          qLen(qLen), dbStartPos(dbStartPos), dbEndPos(dbEndPos), dbLen(dbLen), strand(strand) {}

        // need for sorting the results
        static bool compareAlignment (const SmallAlignment &first, const SmallAlignment &second) {
            //return (first.eval < second.eval);
            if (first.conterminatedKey < second.conterminatedKey)
                return true;
            if (second.conterminatedKey < first.conterminatedKey)
                return false;
            if (first.strand < second.strand)
                return true;
            if (second.strand < first.strand)
                return false;
            if (first.qStartPos < second.qStartPos)
                return true;
            if (second.qStartPos < first.qStartPos)
                return false;
            if (first.dbStartPos < second.dbStartPos)
                return true;
            if (second.dbStartPos < first.dbStartPos)
                return false;
            if (first.dbKey < second.dbKey)
                return true;
            if (second.dbKey < first.dbKey)
                return false;
            return false;
        }

    };

    std::vector<SmallAlignment> overlappingAlnRes;

#pragma omp parallel
    {
        unsigned int thread_idx = 0;
        std::vector<Matcher::result_t> alnRes;
        std::vector<SmallAlignment> tmpOverlappingAlnRes;

#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress2.updateProgress();
            unsigned int queryKey = reader.getDbKey(i);

            char *data = reader.getData(i, thread_idx);
            alnRes.clear();
            Matcher::readAlignmentResults(alnRes, data, true);

            for (size_t elementIdx = 0; elementIdx < alnRes.size(); elementIdx++) {

                size_t dbStartPos = std::min(alnRes[elementIdx].dbStartPos, alnRes[elementIdx].dbEndPos);
                size_t dbEndPos   = std::max(alnRes[elementIdx].dbStartPos, alnRes[elementIdx].dbEndPos);
                size_t qStartPos = std::min(alnRes[elementIdx].qStartPos, alnRes[elementIdx].qEndPos);
                size_t qEndPos = std::max(alnRes[elementIdx].qStartPos, alnRes[elementIdx].qEndPos);

                char strand  = (alnRes[elementIdx].qStartPos > alnRes[elementIdx].qEndPos ||
                                alnRes[elementIdx].dbStartPos > alnRes[elementIdx].dbEndPos  ) ? '-' : '+';
                if(tree.findOverlapping(makeIntervalPos(alnRes[elementIdx].dbKey, dbStartPos),
                                        makeIntervalPos(alnRes[elementIdx].dbKey, dbEndPos)).size() > 0) {
                    // simple swap, evalue should not be a problem here
                    unsigned int qLen = alnRes[elementIdx].qLen;
                    alnRes[elementIdx].qStartPos = dbStartPos;
                    alnRes[elementIdx].qEndPos = dbEndPos;
                    alnRes[elementIdx].qLen = alnRes[elementIdx].dbLen;
                    alnRes[elementIdx].dbStartPos = qStartPos;
                    alnRes[elementIdx].dbEndPos = qEndPos;
                    alnRes[elementIdx].dbLen = qLen;
                    // 00 <- +
                    // 01 <- -
                    // 10 <- -
                    tmpOverlappingAlnRes.emplace_back(alnRes[elementIdx].dbKey, queryKey, alnRes[elementIdx].score,
                                                      alnRes[elementIdx].seqId, alnRes[elementIdx].eval,
                                                      alnRes[elementIdx].qStartPos, alnRes[elementIdx].qEndPos, alnRes[elementIdx].qLen,
                                                      alnRes[elementIdx].dbStartPos, alnRes[elementIdx].dbEndPos, alnRes[elementIdx].dbLen,
                                                      strand);
                };


                if(tree.findOverlapping(makeIntervalPos(queryKey, qStartPos),
                                        makeIntervalPos(queryKey, qEndPos)).size() > 0) {
                    tmpOverlappingAlnRes.emplace_back(queryKey, alnRes[elementIdx].dbKey, alnRes[elementIdx].score,
                                                      alnRes[elementIdx].seqId, alnRes[elementIdx].eval,
                                                      qStartPos, qEndPos, alnRes[elementIdx].qLen,
                                                      dbStartPos, dbEndPos, alnRes[elementIdx].dbLen,
                                                      strand);
                };
            }
        }

#pragma omp critical
        {
            overlappingAlnRes.insert(overlappingAlnRes.end(), tmpOverlappingAlnRes.begin(), tmpOverlappingAlnRes.end());
        }
    }
    omptl::sort(overlappingAlnRes.begin(), overlappingAlnRes.end(), SmallAlignment::compareAlignment);
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        SmallAlignment prevAln(UINT_MAX, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '+');
        char buffer[4096];
        for (size_t i = 0; i < overlappingAlnRes.size(); i++) {
            if(prevAln.conterminatedKey != overlappingAlnRes[i].conterminatedKey){
                if(i > 1){
                    writer.writeEnd(prevAln.conterminatedKey, thread_idx);
                }
                writer.writeStart(thread_idx);
            }

            if (prevAln.conterminatedKey == overlappingAlnRes[i].conterminatedKey &&
                prevAln.dbKey == overlappingAlnRes[i].dbKey &&
                prevAln.qStartPos == overlappingAlnRes[i].qStartPos &&
                prevAln.qEndPos == overlappingAlnRes[i].qEndPos &&
                prevAln.dbStartPos == overlappingAlnRes[i].dbStartPos &&
                prevAln.dbEndPos == overlappingAlnRes[i].dbEndPos) {
                prevAln = overlappingAlnRes[i];
                continue;
            }
            int alnLen = overlappingAlnRes[i].qEndPos - overlappingAlnRes[i].qStartPos;
            std::string backtrace = SSTR(alnLen + 1) + "M";
            Matcher::result_t res(overlappingAlnRes[i].dbKey, overlappingAlnRes[i].score,
                                  0.0, 0.0,
                                  overlappingAlnRes[i].seqId, overlappingAlnRes[i].eval,
                                  0,
                                  overlappingAlnRes[i].qStartPos,
                                  overlappingAlnRes[i].qEndPos,
                                  overlappingAlnRes[i].qLen,
                                  overlappingAlnRes[i].dbStartPos,
                                  overlappingAlnRes[i].dbEndPos,
                                  overlappingAlnRes[i].dbLen, backtrace);
            size_t len = Matcher::resultToBuffer(buffer, res, true, false);
            //res.dbKey=overlappingAlnRes.
            writer.writeAdd(buffer, len, thread_idx);
            prevAln = overlappingAlnRes[i];
        }
        if(prevAln.conterminatedKey!=UINT_MAX){
            writer.writeEnd(prevAln.conterminatedKey, thread_idx);
        }
    }


    delete[] ancestorTax2int;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}

