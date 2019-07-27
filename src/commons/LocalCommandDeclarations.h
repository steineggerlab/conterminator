#ifndef LOCALCOMMANDDECLARATIONS_H
#define LOCALCOMMANDDECLARATIONS_H

#include "Command.h"

extern int distanceton(int argc, const char** argv, const Command &command);
extern int conterminator(int argc, const char** argv, const Command &command);
extern int multipletaxas(int argc, const char** argv, const Command &command);
extern int extractalignments(int argc, const char** argv, const Command &command);
extern int createstats(int argc, const char** argv, const Command &command);
#endif
