#include "Command.h"
#include "LocalCommandDeclarations.h"
#include "LocalParameters.h"

const char* binary_name = "conterminator";
const char* tool_name = "conterminator";
const char* tool_introduction = "Search for Contamination.";
const char* main_author = "Martin Steinegger (themartinsteinegger@gmail.com)";
const char* show_extended_help = "1";
const char* show_bash_info = "1";
bool hide_base_commands = true;

LocalParameters& localPar = LocalParameters::getLocalInstance();
std::vector<struct Command> commands = {
        // Main tools  (for non-experts)
        {"dna",             conterminatordna,             &localPar.conterminatordna,             COMMAND_MAIN,
                "Searches for cross taxon contamination in DNA sequences",
                "Searches for cross taxon contamination in DNA sequences",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:fasta/q> <i:mappingFile> <o:result> <tmpDir>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                 {"mappingFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                 {"result", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                 {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"protein",             conterminatorprotein,             &localPar.conterminatorprotein,             COMMAND_MAIN,
                "Searches for cross taxon contamination in protein sequences",
                "Searches for cross taxon contamination in protein sequences",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:fasta/q> <i:mappingFile> <o:result> <tmpDir>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                 {"mappingFile", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                 {"result", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::flatfile },
                 {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"predictcontamination",          predictcontamination,          &localPar.createstats,         COMMAND_HIDDEN,
                "Predict contaminated taxon",
                "Predict contaminated taxon",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:resultDB> <o:resultDB>", CITATION_MMSEQS2,
                {{"resultDB",   DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                 {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"createstats",          createstats,          &localPar.createstats,         COMMAND_HIDDEN,
                "Create taxon statistic",
                "Create taxon statistic",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequeceDB> <i:resultDB> <o:resultDB>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                        {"alnDB",   DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                        {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"createallreport",          createallreport,          &localPar.createstats,         COMMAND_HIDDEN,
                "Create taxon statistic",
                "Create taxon statistic",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequeceDB> <i:resultDB> <o:resultDB>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                        {"alnDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                        {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"crosstaxonfilterorf",          crosstaxonfilterorf,          &localPar.crosstaxonfilterorf,         COMMAND_HIDDEN,
                "Extract cluster with n taxas with orf indirection",
                "Extract cluster with n taxas with orf indirection",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequenceDB> <i:clusterDB> <o:clusterDB>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                 {"orfHeader", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA, &DbValidator::genericDb },
                 {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb },
                 {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::resultDb }}},
        {"crosstaxonfilter",          crosstaxonfilter,          &localPar.extractalignments,         COMMAND_HIDDEN,
                "Extract cluster with n taxas",
                "Extract cluster with n taxas",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequenceDB> <i:clusterDB> <o:clusterDB>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                  {"clusterDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb },
                  {"clusterDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::clusterDb }}},
        {"extractalignments",          extractalignments,          &localPar.extractalignments,         COMMAND_HIDDEN,
                "Extract alignments containing n taxas",
                "Extract alignments containing n taxas",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequenceDB> <i:alnDB> <o:alnDB>",CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                 {"alnDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                 {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}}

};
