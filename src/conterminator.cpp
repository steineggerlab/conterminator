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
        {"conterminator",             conterminatorworkflow,             &localPar.conterminatorworkflow,             COMMAND_MAIN,
                "Searches for contamination and creates a masking file",
                "Searches for contamination and creates a masking file.",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequenceDB> <o:alnDB> <tmpDir>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                {"alnDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                {"tmpDir", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::directory }}},
        {"createstats",          createstats,          &localPar.createstats,         COMMAND_TAXONOMY,
                "Create taxon statistic",
                "Create taxon statistic",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequeceDB> <i:resultDB> <o:resultDB>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                 {"alnDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                 {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"createallreport",          createallreport,          &localPar.createstats,         COMMAND_TAXONOMY,
                "Create taxon statistic",
                "Create taxon statistic",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:sequeceDB> <i:resultDB> <o:resultDB>", CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                        {"alnDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                        {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}},
        {"extractalignments",          extractalignments,          &localPar.extractalignments,         COMMAND_TAXONOMY,
                "Extract alignments containing n taxas",
                "Extract alignments containing n taxas",
                "Martin Steinegger <martin.steinegger@mpibpc.mpg.de>",
                "<i:resultDB> <o:resultDB>",CITATION_MMSEQS2,
                {{"sequenceDB", DbType::ACCESS_MODE_INPUT, DbType::NEED_DATA|DbType::NEED_TAXONOMY, &DbValidator::taxSequenceDb },
                {"alnDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb },
                {"resultDB",   DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::alignmentDb }}}

};
