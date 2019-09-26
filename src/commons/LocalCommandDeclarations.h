#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int conterminatordna(int argc, const char** argv, const Command &command);
extern int conterminatorprotein(int argc, const char** argv, const Command &command);
extern int predictcontamination(int argc, const char** argv, const Command &command);
extern int extractalignments(int argc, const char** argv, const Command &command);
extern int crosstaxonfilter(int argc, const char **argv, const Command &command);
extern int createstats(int argc, const char** argv, const Command &command);
extern int createallreport(int argc, const char** argv, const Command &command);
#endif
