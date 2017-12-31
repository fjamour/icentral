#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/_bc.o \
	${OBJECTDIR}/src/_bc_mem.o \
	${OBJECTDIR}/src/bc.o \
	${OBJECTDIR}/src/bcc_delta.o \
	${OBJECTDIR}/src/bcc_scratch.o \
	${OBJECTDIR}/src/bicon.o \
	${OBJECTDIR}/src/experiments.o \
	${OBJECTDIR}/src/graph_hash_t.o \
	${OBJECTDIR}/src/graph_t.o \
	${OBJECTDIR}/src/iter_info_t.o \
	${OBJECTDIR}/src/main.o \
	${OBJECTDIR}/src/mcb_find.o \
	${OBJECTDIR}/src/paper_exp.o \
	${OBJECTDIR}/src/qube.o \
	${OBJECTDIR}/src/subgraph_t.o \
	${OBJECTDIR}/src/utility.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-w -pthread -std=c++11
CXXFLAGS=-w -pthread -std=c++11

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fjbc_mpi

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fjbc_mpi: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fjbc_mpi ${OBJECTFILES} ${LDLIBSOPTIONS} -pthread -std=c++11

${OBJECTDIR}/src/_bc.o: src/_bc.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/_bc.o src/_bc.cc

${OBJECTDIR}/src/_bc_mem.o: src/_bc_mem.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/_bc_mem.o src/_bc_mem.cc

${OBJECTDIR}/src/bc.o: src/bc.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bc.o src/bc.cc

${OBJECTDIR}/src/bcc_delta.o: src/bcc_delta.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bcc_delta.o src/bcc_delta.cc

${OBJECTDIR}/src/bcc_scratch.o: src/bcc_scratch.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bcc_scratch.o src/bcc_scratch.cc

${OBJECTDIR}/src/bicon.o: src/bicon.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bicon.o src/bicon.cc

${OBJECTDIR}/src/experiments.o: src/experiments.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/experiments.o src/experiments.cc

${OBJECTDIR}/src/graph_hash_t.o: src/graph_hash_t.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/graph_hash_t.o src/graph_hash_t.cc

${OBJECTDIR}/src/graph_t.o: src/graph_t.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/graph_t.o src/graph_t.cc

${OBJECTDIR}/src/iter_info_t.o: src/iter_info_t.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/iter_info_t.o src/iter_info_t.cc

${OBJECTDIR}/src/main.o: src/main.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/main.o src/main.cc

${OBJECTDIR}/src/mcb_find.o: src/mcb_find.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/mcb_find.o src/mcb_find.cc

${OBJECTDIR}/src/paper_exp.o: src/paper_exp.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/paper_exp.o src/paper_exp.cc

${OBJECTDIR}/src/qube.o: src/qube.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/qube.o src/qube.cc

${OBJECTDIR}/src/subgraph_t.o: src/subgraph_t.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/subgraph_t.o src/subgraph_t.cc

${OBJECTDIR}/src/utility.o: src/utility.cc 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/utility.o src/utility.cc

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fjbc_mpi

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
