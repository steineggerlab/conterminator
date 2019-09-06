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
    cmd.addVariable("LINCLUST_PAR", par.createParameterString(par.linclustworkflow).c_str());
    cmd.addVariable("CROSSTAXA_PAR", par.createParameterString(par.extractalignments).c_str());
    cmd.addVariable("CREATESTATS_PAR", par.createParameterString(par.createstats).c_str());
    cmd.addVariable("ALN_PAR", par.createParameterString(par.align).c_str());

    FileUtil::writeFile(tmpDir + "/conterminatorprotein.sh", conterminatorprotein_sh, conterminatorprotein_sh_len);
    std::string program(tmpDir + "/conterminatorprotein.sh");
    cmd.execProgram(program.c_str(), par.filenames);

    return EXIT_SUCCESS;
}
