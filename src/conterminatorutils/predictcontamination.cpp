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


int predictcontamination(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    // bacteria, archaea, eukaryotic, virus
    par.parseParameters(argc, argv, command, true, 0, 0);

    DBReader<unsigned int> reader(par.db1.c_str(), par.db1Index.c_str(), par.threads,
                                  DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX);
    reader.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    DBWriter writer(par.db2.c_str(), par.db2Index.c_str(), par.threads, par.compressed, reader.getDbtype());
    writer.open();

    int lenThreshold = 20000;

    Debug::Progress progress(reader.getSize());
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = (unsigned int) omp_get_thread_num();
#endif
        const char * entry[255];
        std::string resultData;
        int * termLen = new int[256];
        int * termCount = new int[256];
        std::map<int, std::string> longestId;
        std::map<int, std::string> longestSpeciesName;
        std::set<std::string> idDetected;
#pragma omp for schedule(dynamic, 10)
        for (size_t i = 0; i < reader.getSize(); ++i) {
            progress.updateProgress();
            memset(termLen, 0, sizeof(int) * 256);
            memset(termCount, 0, sizeof(int) * 256);
            longestId.clear();
            longestSpeciesName.clear();
            idDetected.clear();
            resultData.clear();
            unsigned int queryKey = reader.getDbKey(i);
            char *data = reader.getData(i, thread_idx);
            char *firstData = data;
            int maxTermId = -1;
            while (*firstData != '\0') {
                // NZ_BFCN01000401.1       244     392     453     453     0       Escherichia coli
                const size_t columns = Util::getWordsOfLine(firstData, entry, 255);
                if (columns == 0) {
                    firstData = Util::skipLine(firstData);
                    continue;
                }
                size_t len = entry[1] - entry[0];
                const std::string fastaId(entry[0], len-1);
                int startPos = Util::fast_atoi<int>(entry[1]);
                int endPos = Util::fast_atoi<int>(entry[2]);
                int nLen = Util::fast_atoi<int>(entry[3]);
                int entryLen = Util::fast_atoi<int>(entry[4]);
                int termId = Util::fast_atoi<int>(entry[5]);

                int adjustedLen = ((entryLen - endPos) < lenThreshold) ? (entryLen - std::min(startPos, endPos)) : nLen;
                adjustedLen = ((entryLen - endPos) < lenThreshold) ? std::max(startPos, endPos) : adjustedLen;
                if (termLen[termId] == 0 || adjustedLen > termLen[termId]) {
                    termLen[termId] = adjustedLen;
                    longestId[termId] = fastaId;
                    char * nextLine = Util::skipLine(firstData);
                    std::string speciesName(entry[6], (nextLine - entry[6] - 1) );
                    longestSpeciesName[termId] = speciesName;
                }
                const bool existsAlready = idDetected.find(fastaId) != idDetected.end();
                if (existsAlready == false) {
                    idDetected.insert(fastaId);
                    termCount[termId]++;
                }
                maxTermId = std::max(termId, maxTermId);
                firstData = Util::skipLine(firstData);
            }

            // predict direction of contamination based on entry length
            int contermCnt = 0;
            int notContermCnt = 0;
            int notContermId = -1;
            for (int term = 0; term <= maxTermId; term++) {
                if(termLen[term] > lenThreshold){
                    notContermCnt++;
                    notContermId = term;
                } else if (termLen[term] != 0 && termLen[term] <= lenThreshold){
                    contermCnt++;
                }
            }
            if (notContermCnt == 1 && contermCnt >= 1) {
                char *secondData = data;
                while (*secondData != '\0') {
                    const size_t columns = Util::getWordsOfLine(secondData, entry, 255);
                    if (columns == 0) {
                        secondData = Util::skipLine(secondData);
                        continue;
                    }

                    int termId = Util::fast_atoi<int>(entry[5]);

                    // conterm
                    if(termId != notContermId){
                        size_t len = entry[1] - entry[0];
                        std::string fastaId(entry[0], len-1);
                        int startPos = Util::fast_atoi<int>(entry[1]);
                        int endPos = Util::fast_atoi<int>(entry[2]);
                        int nLen = Util::fast_atoi<int>(entry[3]);
                        int entryLen = Util::fast_atoi<int>(entry[4]);
                        char * nextLine = Util::skipLine(secondData);
                        std::string speciesName(entry[6], (nextLine - entry[6] - 1) );

                        int adjustedLen = ((entryLen - endPos) < lenThreshold) ? (entryLen - std::min(startPos, endPos)) : nLen;
                        adjustedLen = ((entryLen - endPos) < lenThreshold) ? std::max(startPos, endPos) : adjustedLen;
                        resultData.append(fastaId);
                        resultData.push_back('\t');
                        resultData.append(SSTR(termId));
                        resultData.push_back('\t');
                        resultData.append(speciesName);
                        resultData.push_back('\t');
                        resultData.append(SSTR(startPos));
                        resultData.push_back('\t');
                        resultData.append(SSTR(endPos));
                        resultData.push_back('\t');
                        resultData.append(SSTR(adjustedLen));
                        resultData.push_back('\t');
                        resultData.append(SSTR(longestId[notContermId]));
                        resultData.push_back('\t');
                        resultData.append(SSTR(notContermId));
                        resultData.push_back('\t');
                        resultData.append(longestSpeciesName[notContermId]);
                        resultData.push_back('\t');
                        resultData.append(SSTR(termLen[notContermId]));
                        resultData.push_back('\t');
                        resultData.append(SSTR(termCount[notContermId]));
                        resultData.push_back('\n');
                    }
                    secondData = Util::skipLine(secondData);
                }
                if(resultData.size() > 0){
                    writer.writeData(resultData.c_str(), resultData.size(), queryKey, thread_idx);
                }

            }
        }
        delete [] termLen;
    }

    writer.close();
    reader.close();
    return EXIT_SUCCESS;
}


