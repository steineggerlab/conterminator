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
    p->alnLenThr = 100;
    p->orfMaxLength = 32734;
    p->forwardFrames= "1";
    p->reverseFrames= "1";
    p->maxSeqLen = 1000;
    p->covThr = 0.0;
    p->sequenceOverlap = 0;
    p->kmersPerSequence = 30;
    p->rescoreMode = 2;
}

int conterminator(int argc, const char **argv, const Command &command) {
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

    cmd.addVariable("SPLITSEQ_PAR", par.createParameterString(par.splitsequence).c_str());
    cmd.addVariable("RESCORE_DIAGONAL_PAR", par.createParameterString(par.rescorediagonal).c_str());
    cmd.addVariable("OFFSETALIGNMENT_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("MULTIPLETAXA_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("DISTANCETON_PAR", par.createParameterString(par.onlythreads).c_str());
    cmd.addVariable("KMERMATCHER_PAR", par.createParameterString(par.kmermatcher).c_str());
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
