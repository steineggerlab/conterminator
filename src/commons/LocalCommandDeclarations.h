#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int conterminatorworkflow(int argc, const char** argv, const Command &command);
extern int extractalignments(int argc, const char** argv, const Command &command);
extern int createstats(int argc, const char** argv, const Command &command);
extern int createallreport(int argc, const char** argv, const Command &command);
#endif
