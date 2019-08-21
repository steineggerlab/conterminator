#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

#include "conterminator.sh.h"

void setConterminatorWorkflowDefaults(LocalParameters *p) {
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->addBacktrace = true;
    p->sensitivity = 5.7;
    p->evalThr = 0.001;
    //p->orfLongest = true;
    p->orfStartMode = 1;
    p->orfMinLength = 30;
    p->seqIdThr = 0.9;
    p->maskMode = 0;
    p->kmerSize = 15;
    p->alnLenThr = 100;
    p->orfMaxLength = 32734;
    p->forwardFrames= "1";
    p->reverseFrames= "1";
    p->maxSeqLen = 1000;
    p->covThr = 0.0;
    p->sequenceOverlap = 0;
    p->kmersPerSequence = 100;
    p->rescoreMode = 2;
    p->minDiagScoreThr = p->alnLenThr;
    // leave ungapped alignment untouched
    if(p->alignmentMode != Parameters::ALIGNMENT_MODE_UNGAPPED){
        p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    }
    //p->orfLongest = true;
    p->exactKmerMatching = true;
//    if ( p->PARAM_DIAGONAL_SCORING.wasSet == false) {
//        p->diagonalScoring = 0;
//    }
    if ( p->PARAM_STRAND.wasSet == false) {
        p->strand = 2;
    }
    if ( p->PARAM_K.wasSet == false) {
        p->kmerSize = 15;
    }
    if (  p->PARAM_MAX_SEQ_LEN.wasSet == false) {
        p->maxSeqLen = 10000;
    }
    if( p->PARAM_GAP_OPEN.wasSet == false){
        p->gapOpen = 5;
    }
    if( p->PARAM_GAP_EXTEND.wasSet  == false){
        p->gapExtend = 2;
    }
    // (Bacteria, Archaea), Fungi, Animalia, Plantae, Rest of Eukaryota
    p->taxonList = "(2|2157),4751,33208,33090,(2759&!4751&!33208&!33090)";
    // virus, unclassified sequences, other sequences,  artifical sequences, retro virus, environmental samples
    p->blacklist = "10239,12908,28384,81077,11632,340016,61964,48479,48510";
}

int conterminatorworkflow(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setConterminatorWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_COMMON);

    CommandCaller cmd;
    std::string tmpPath = par.filenames.back();
    if (FileUtil::directoryExists(tmpPath.c_str()) == false) {
        Debug(Debug::INFO) << "Temporary folder " << tmpPath << " does not exist or is not a directory.\n";
        if (FileUtil::makeDir(tmpPath.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not crate tmp folder " << tmpPath << ".\n";
            return EXIT_FAILURE;
        } else {
            Debug(Debug::INFO) << "Created directory " << tmpPath << "\n";
        }
    }
    size_t hash = par.hashParameter(par.filenames, par.conterminatorworkflow);
    std::string tmpDir = tmpPath + "/" + SSTR(hash);
    if (FileUtil::directoryExists(tmpDir.c_str()) == false) {
        if (FileUtil::makeDir(tmpDir.c_str()) == false) {
            Debug(Debug::ERROR) << "Could not create sub folder in temporary directory " << tmpDir << ".\n";
            return EXIT_FAILURE;
        }
    }
    par.filenames.pop_back();
    FileUtil::symlinkAlias(tmpDir, "latest");
    char *p = realpath(tmpDir.c_str(), NULL);
    if (p == NULL) {
        Debug(Debug::ERROR) << "Could not get real path of " << tmpDir << "!\n";
        EXIT(EXIT_FAILURE);
    }
    par.filenames.push_back(p);
    free(p);


    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("EXTRACTALIGNMENTS_PAR", par.createParameterString(par.extractalignments).c_str());
    cmd.addVariable("CREATESTATS_PAR", par.createParameterString(par.createstats).c_str());

    cmd.addVariable("SPLITSEQ_PAR", par.createParameterString(par.splitsequence).c_str());
    cmd.addVariable("RESCORE_DIAGONAL_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("EXTRACT_FRAMES_PAR", par.createParameterString(par.extractframes).c_str());


    par.kmerSize = 24;
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());
    par.kmerSize = 15;
    par.maxSeqLen = 1000000;
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());

//    par.alphabetSize = alphabetSize;
//    par.kmerSize = kmerSize;
//
//    cmd.addVariable("ALIGNMENT_PAR", par.createParameterString(par.align).c_str());
//    cmd.addVariable("CLUSTER_PAR", par.createParameterString(par.clust).c_str());


    FileUtil::writeFile(tmpDir + "/conterminator.sh", conterminator_sh, conterminator_sh_len);
    std::string program(tmpDir + "/conterminator.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
