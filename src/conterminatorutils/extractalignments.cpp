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


int extractalignments(int argc, const char **argv, const Command& command) {
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
        Contamination(unsigned int key, int start, int end, unsigned int len)
        : key(key), start(start), end(end), len(len) {}
        Contamination(){};
        unsigned int key;
        int start;
        int end;
        unsigned int len;

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
        IntervalArray ** speciesRanges = new IntervalArray*[taxonList.size()];
        for(size_t i = 0; i < taxonList.size(); i++){
            speciesRanges[i] = new IntervalArray();
        }
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
            unsigned int queryLen = reader.getSeqLens(i);

            unsigned int queryTaxon = Multipletaxas::getTaxon(queryKey, mapping);
            if(queryTaxon == 0 || queryTaxon == UINT_MAX ){
                continue;
            }
            unsigned int queryAncestorTaxon = UINT_MAX;

            for (size_t j = 0; j < taxonList.size(); ++j) {
                if (t.IsAncestor(taxonList[j], queryTaxon)) {
                    queryAncestorTaxon = taxonList[j];
                    break;
                }
            }

            if (queryAncestorTaxon == UINT_MAX) {
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

            // find max. taxa
            for (size_t taxId = 0; taxId < taxonList.size(); taxId++) {
                bool hasTaxa = (taxaCounter[taxId] > 0);
                distinctTaxaCnt += hasTaxa;
            }
            if (distinctTaxaCnt > 1) {
                for (size_t i = 0; i < taxonList.size(); i++) {
                    speciesRanges[i]->reset();
                }

                // fill up interval tree with elements
                for (size_t elementIdx = 0; elementIdx < elements.size(); elementIdx++) {
                    if (static_cast<unsigned int>(elements[elementIdx].ancestorTax) != queryAncestorTaxon) {
                        Matcher::result_t res = Matcher::parseAlignmentRecord(elements[elementIdx].data, true);
                        speciesRanges[ancestorTax2int[elements[elementIdx].ancestorTax]]->insert(res.qStartPos,
                                                                                                 res.qEndPos);
                    }
                }
                for (size_t i = 0; i < taxonList.size(); i++) {
                    speciesRanges[i]->buildRanges();
                }

                for (size_t i = 0; i < taxonList.size(); i++) {
                    for (size_t j = 0; j < speciesRanges[i]->getRangesSize(); j++) {
                        IntervalArray::Range range = speciesRanges[i]->getRange(j);
                        privateContaminations.push_back(Contamination(queryKey, range.start, range.end, queryLen));
                    }
                }
            }
        }

#pragma omp critical
        {
            allContaminations.insert(allContaminations.end(), privateContaminations.begin(), privateContaminations.end());
        };
        delete[] taxaCounter;

        for(size_t i = 0; i < taxonList.size(); i++){
            delete speciesRanges[i];
        }
        delete [] speciesRanges;
    }

    omptl::sort(allContaminations.begin(), allContaminations.end(), Contamination::compareContaminationByKeyStartEnd);
    size_t writePos = -1;
    unsigned int prevKey = UINT_MAX;

    for(size_t i = 0; i < allContaminations.size(); i++){
        if(prevKey != allContaminations[i].key){
            writePos++;
            allContaminations[writePos]=allContaminations[i];
        }else{
            if(allContaminations[i].start <= (allContaminations[writePos].end+1)){
                allContaminations[writePos].end = std::max(allContaminations[i].end, allContaminations[writePos].end);
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
            writer.writeData(buffer, len, allContaminations[i].key, 0);
        }
    }

    delete[] ancestorTax2int;
    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}

