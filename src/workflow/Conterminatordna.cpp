#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "conterminatordna.sh.h"

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

int conterminatordna(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setConterminatorWorkflowDefaults(&par);
    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_COMMON);

    CommandCaller cmd;
    std::string tmpDir = par.db3;
    std::string hash = SSTR(par.hashParameter(par.filenames, par.linclustworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);

    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);

    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("EXTRACTALIGNMENTS_PAR", par.createParameterString(par.extractalignments).c_str());
    cmd.addVariable("CREATESTATS_PAR", par.createParameterString(par.createstats).c_str());
    int prevCompressed = par.compressed;
    par.compressed = 1;
    cmd.addVariable("SPLITSEQ_PAR", par.createParameterString(par.splitsequence).c_str());
    par.compressed = prevCompressed;
    cmd.addVariable("RESCORE_DIAGONAL1_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("THREADS_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("EXTRACT_FRAMES_PAR", par.createParameterString(par.extractframes).c_str());

    par.kmerSize = 24;
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());
    par.kmerSize = 15;
    par.maxSeqLen = 1000000;
    par.maskMode = 1;
    cmd.addVariable("PREFILTER_PAR", par.createParameterString(par.prefilter).c_str());
    float tmpSeqIdThr = par.seqIdThr;
    par.seqIdThr = sqrt(par.seqIdThr);
    cmd.addVariable("RESCORE_DIAGONAL2_PAR", par.createParameterString(par.rescorediagonal).c_str());
    par.seqIdThr = tmpSeqIdThr;

    FileUtil::writeFile(tmpDir + "/conterminatordna.sh", conterminatordna_sh, conterminatordna_sh_len);
    std::string program(tmpDir + "/conterminatordna.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
