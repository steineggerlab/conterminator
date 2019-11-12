#include "DBReader.h"
#include "Util.h"
#include "CommandCaller.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "conterminatorprotein.sh.h"

void setConterminatorProteinDefaults(LocalParameters *p) {
    //p->orfLongest = true;
    p->covThr = 0.8;
    p->maskMode = 0;
    p->evalThr = 0.001;
    p->seqIdThr = 0.95;
    p->removeTmpFiles = true;
    p->covThr = 0.95;
    p->alignmentMode = Parameters::ALIGNMENT_MODE_SCORE_COV_SEQID;
    p->covMode = Parameters::COV_MODE_BIDIRECTIONAL;
    p->kmersPerSequence = 100;
    // (Bacteria, Archaea), Fungi, Animalia, Plantae, Rest of Eukaryota
    p->taxonList = "(2|2157),2759";
    // virus, unclassified sequences, other sequences,  artifical sequences, retro virus, environmental samples
    p->blacklist = "10239,12908,28384,81077,11632,340016,61964,48479,48510";
}

int conterminatorprotein(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setConterminatorProteinDefaults(&par);
    for(size_t i = 0; i < par.conterminatorprotein.size(); i++){
        par.conterminatorprotein[i]->category |= MMseqsParameter::COMMAND_EXPERT;
    }
    par.PARAM_TAXON_LIST.category = MMseqsParameter::COMMAND_MISC;
    par.PARAM_BLACKLIST.category = MMseqsParameter::COMMAND_MISC;
    par.PARAM_KMER_PER_SEQ.category = MMseqsParameter::COMMAND_PREFILTER;
    par.PARAM_MIN_SEQ_ID.category = MMseqsParameter::COMMAND_ALIGN;
    par.PARAM_C.category = MMseqsParameter::COMMAND_ALIGN;
    par.PARAM_COV_MODE.category = MMseqsParameter::COMMAND_ALIGN;
    par.PARAM_THREADS.category = MMseqsParameter::COMMAND_COMMON;
    par.PARAM_V.category = MMseqsParameter::COMMAND_COMMON;
    par.PARAM_SPLIT_MEMORY_LIMIT.category = MMseqsParameter::COMMAND_COMMON;

    par.parseParameters(argc, argv, command, true, 0, MMseqsParameter::COMMAND_COMMON);

    CommandCaller cmd;
    std::string tmpDir = par.db4;
    std::string hash = SSTR(par.hashParameter(par.filenames, par.linclustworkflow));
    if (par.reuseLatest) {
        hash = FileUtil::getHashFromSymLink(tmpDir + "/latest");
    }
    tmpDir = FileUtil::createTemporaryDirectory(tmpDir, hash);

    par.filenames.pop_back();
    par.filenames.push_back(tmpDir);
    cmd.addVariable("CREATEDB_PAR", par.createParameterString(par.createdb).c_str());
    cmd.addVariable("TAXMAPPINGFILE", par.db2.c_str());
    cmd.addVariable("REMOVE_TMP", par.removeTmpFiles ? "TRUE" : NULL);
    cmd.addVariable("RUNNER", par.runner.c_str());
    cmd.addVariable("LINCLUST_PAR", par.createParameterString(par.linclustworkflow).c_str());
    cmd.addVariable("CROSSTAXA_PAR", par.createParameterString(par.extractalignments).c_str());
    cmd.addVariable("CREATESTATS_PAR", par.createParameterString(par.createstats).c_str());
    cmd.addVariable("ALN_PAR", par.createParameterString(par.align).c_str());

    FileUtil::writeFile(tmpDir + "/conterminatorprotein.sh", conterminatorprotein_sh, conterminatorprotein_sh_len);
    std::string program(tmpDir + "/conterminatorprotein.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
