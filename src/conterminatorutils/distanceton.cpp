//
// Created by Martin Steinegger on 28.10.18.
//

#include "Parameters.h"
#include "Util.h"
#include "Debug.h"
#include "DBReader.h"
#include "DBWriter.h"
#include "SubstitutionMatrix.h"
#include "MathUtil.h"

#include <fstream>
#include <iomanip>
#include "Matcher.h"

void findLeftAndRightPos(int startPos, int endPos, char *seq, int &leftNPos, int &rightNPos);

#ifdef OPENMP
#include <omp.h>
#endif


int measureDistanceToN(Parameters &par) {
    Debug(Debug::INFO) << "Query file: " << par.db1 << "\n";
    DBReader<unsigned int> qdbr(par.db1.c_str(), par.db1Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    qdbr.open(DBReader<unsigned int>::NOSORT);
    if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
        qdbr.readMmapedDataInMemory();
    }

    bool sameDB = false;
    Debug(Debug::INFO) << "Target file: " << par.db2 << "\n";
    DBReader<unsigned int> *tdbr = NULL;
    if (par.db1.compare(par.db2) == 0) {
        sameDB = true;
        tdbr = &qdbr;
    } else {
        tdbr = new DBReader<unsigned int>(par.db2.c_str(), par.db2Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA | DBReader<unsigned int>::USE_INDEX );
        tdbr->open(DBReader<unsigned int>::NOSORT);
        if (par.preloadMode != Parameters::PRELOAD_MODE_MMAP) {
            qdbr.readMmapedDataInMemory();
        }
    }

    Debug(Debug::INFO) << "Alignment database: " << par.db3 << "\n";
    DBReader<unsigned int> alndbr(par.db3.c_str(), par.db3Index.c_str(), par.threads, DBReader<unsigned int>::USE_DATA|DBReader<unsigned int>::USE_INDEX);
    alndbr.open(DBReader<unsigned int>::LINEAR_ACCCESS);

    Debug(Debug::INFO) << "Start writing file to " << par.db4 << "\n";
    DBWriter dbw(par.db4.c_str(), par.db4Index.c_str(), static_cast<unsigned int>(par.threads), par.compressed, alndbr.getDbtype());
    dbw.open();

    const char newline = '\n';
#pragma omp parallel
    {
        unsigned int thread_idx = 0;
#ifdef OPENMP
        thread_idx = static_cast<unsigned int>(omp_get_thread_num());
#endif
        char buffer[1024+32768];
        std::vector<Matcher::result_t> results;
        results.reserve(300);
#pragma omp for schedule(dynamic, 1000)
        for (size_t i = 0; i < alndbr.getSize(); i++) {
            Debug::printProgress(i);

            unsigned int queryKey = alndbr.getDbKey(i);
            char *qSeq = NULL;
            size_t qSeqLen;
            if (par.extractMode == Parameters::EXTRACT_QUERY) {
                size_t id = qdbr.getId(queryKey);
                qSeqLen = qdbr.getSeqLens(id);
                qSeq = qdbr.getData(id, thread_idx);
            }
            dbw.writeStart(thread_idx);

            char *data = alndbr.getData(i, thread_idx);
            Matcher::readAlignmentResults(results, data);
            for (size_t j = 0; j < results.size(); j++) {
                Matcher::result_t &res = results[j];
                size_t length = 0;
                char *seq = NULL;
                int leftNPos = -1;
                int rightNPos = -1;
                if (qSeq) {
                    findLeftAndRightPos(std::min(res.qStartPos, res.qEndPos), std::max(res.qEndPos, res.qStartPos),
                            qSeq, leftNPos, rightNPos);
                    res.dbStartPos = leftNPos;
                    res.dbEndPos = rightNPos;
                    res.dbLen = (rightNPos == -1 ) ?  qSeqLen - std::max(0, leftNPos) : rightNPos - std::max(0, leftNPos);
                    res.qLen = qSeqLen;
                } else if (par.extractMode == Parameters::EXTRACT_TARGET) {
                    size_t id = tdbr->getId(res.dbKey);
                    size_t seqLen = tdbr->getSeqLens(id);
                    seq = tdbr->getData(id, thread_idx);
                    findLeftAndRightPos(std::min(res.dbStartPos, res.dbEndPos), std::max(res.dbStartPos, res.dbEndPos),
                            seq, leftNPos, rightNPos);
                    res.qStartPos = leftNPos;
                    res.qEndPos = rightNPos;
                    res.qLen = (rightNPos == -1 ) ?  seqLen - std::max(0, leftNPos) : rightNPos - std::max(0, leftNPos);
                    res.dbLen = seqLen;
                } else {
                    Debug(Debug::ERROR) << "Missing extraction type!\n";
                    EXIT(EXIT_FAILURE);
                }
                size_t len = Matcher::resultToBuffer(buffer, res, false);
                dbw.writeAdd(buffer, len, thread_idx);
            }
            dbw.writeEnd(queryKey, thread_idx);
            results.clear();
        }
    }
    dbw.close();
    alndbr.close();
    qdbr.close();
    if (sameDB == false) {
        tdbr->close();
        delete tdbr;
    }

    return EXIT_SUCCESS;
}

void findLeftAndRightPos(int startPos, int endPos, char *seq, int &leftNPos, int &rightNPos) {
    for(int pos = startPos; pos > 0; pos--){
        if(seq[pos] == 'N'){
            leftNPos = pos;
            break;
        }
    }
    for(int pos = endPos; seq[pos] != '\n'; pos++) {
        if(seq[pos] == 'N'){
            rightNPos = pos;
            break;
        }
    }
}


int distanceton(int argc, const char **argv, const Command& command) {
    Parameters& par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, 4);

    // never allow deletions
    par.allowDeletion = false;

    int retCode = measureDistanceToN(par);

    return retCode;
}
