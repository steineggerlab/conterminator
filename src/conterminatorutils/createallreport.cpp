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


void findLeftAndRightPos(int startPos, int endPos, char *seq, int &leftNPos, int &rightNPos) {
    for(int pos = startPos; pos > -1; pos--){
        leftNPos = pos;
        if(seq[pos] == 'N'){
            break;
        }
    }
    for(int pos = endPos; seq[pos] != '\n'; pos++) {
        rightNPos = pos;
        if(seq[pos] == 'N'){
            break;
        }
    }
}


int createallreport(int argc, const char **argv, const Command& command) {
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

    DBReader<unsigned int> header(par.hdr1.c_str(), par.hdr1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    header.open(DBReader<unsigned int>::NOSORT);
    header.readMmapedDataInMemory();

    DBReader<unsigned int> sequences(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_INDEX|DBReader<unsigned int>::USE_DATA);
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
            TaxonUtils::assignTaxonomy(elements, data, mapping, *t, taxonomyExpression, blackList, taxaCounter, true);
            std::sort(elements.begin(), elements.end(), TaxonUtils::TaxonInformation::compareByDbKeyAndStart);

            size_t writePos = -1;
            unsigned int prevKey = UINT_MAX;

            for (size_t i = 0; i < elements.size(); i++) {
                if (prevKey != elements[i].dbKey) {
                    writePos++;
                    elements[writePos] = elements[i];
                } else {
                    if (elements[i].start <= (elements[writePos].end + 1)) {
                        elements[writePos].end = std::max(elements[i].end, elements[writePos].end);
                    } else {
                        writePos++;
                        elements[writePos] = elements[i];
                    }
                }
                prevKey = elements[i].dbKey;
            }

            std::sort(elements.begin(), elements.begin() + writePos,
                      TaxonUtils::TaxonInformation::compareByTaxAndStart);

            // recount
            for (size_t i = 0; i < writePos; i++) {
                const unsigned int dbkey = elements[i].dbKey;
                const bool existsAlready = idDetected.find(dbkey) != idDetected.end();
                if (existsAlready == false) {
                    idDetected.insert(dbkey);
                    resultData.append(Util::parseFastaHeader(header.getDataByDBKey(dbkey, thread_idx)));
                    resultData.push_back('\t');
                    resultData.append(SSTR(elements[i].start));
                    resultData.push_back('\t');
                    resultData.append(SSTR(elements[i].end));
                    resultData.push_back('\t');
//                    unsigned int id = sequences.getId(dbkey);
//                    char  * seq = sequences.getData(id, thread_idx);
//                    int leftNPos = -1;
//                    int rightNPos = -1;
//                    findLeftAndRightPos(std::min(elements[i].start, elements[i].end), std::max(elements[i].start, elements[i].end),
//                                        seq, leftNPos, rightNPos);
//                    int length = (rightNPos - leftNPos) + 1;
//                    resultData.append(SSTR(length));
//                    resultData.push_back('\t');
                    resultData.append(SSTR(sequences.getSeqLens(sequences.getId(dbkey))-2));
                    resultData.push_back('\t');
                    resultData.append(SSTR(elements[i].termId));
                    resultData.push_back('\t');
                    const TaxonNode *node = t->taxonNode(elements[i].currTaxa, false);
                    if (node == NULL) {
                        resultData.append("Undef");
                    } else {
                        resultData.append(node->name);
                    }
                    resultData.push_back('\n');
                }
            }
            writer.writeData(resultData.c_str(), resultData.size(), queryKey, thread_idx);

        }
    }
    Debug(Debug::INFO) << "\nDetected potentail conterminetaion in the following Taxons: \n" ;

    delete t;
    writer.close();
    reader.close();
    header.close();
    sequences.close();
    return EXIT_SUCCESS;
}

