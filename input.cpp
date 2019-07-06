/*
 * Copyright (c) 2017-2019, Triad National Security, LLC.
 * All rights Reserved.
 * 
 * This is the code released under LANL Copyright Disclosure C17041/LA-CC-17-041
 * Copyright 2017-2019.  Triad National Security, LLC. This material was produced
 * under U.S. Government contract 89233218CNA000001 for Los Alamos National
 * Laboratory (LANL), which is operated by Triad National Security, LLC
 * for the U.S. Department of Energy. See LICENSE file for details.
 *
 * Released under the New BSD License
 *
 * Bob Robey brobey@lanl.gov and Rao Garimella rao@lanl.gov
 */

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

//  Global variables.
const char *PACKAGE_VERSION = "1.0";
char progName[13];      //  Program name.
char progVers[8];       //  Program version.

//  External global variables.
extern bool memory_verbose,
            verbose;
extern int  ncells,
            nmats,
            itermax;

void outputHelp()
{   cout << "MultiMatTest is a testbed for examining different multimaterial data structures." << endl
         << "Version is " << PACKAGE_VERSION << endl << endl
         << "Usage:  " << progName << " [options]..." << endl
         << "  -d                detailed memory output;" << endl
         << "  -h                display this help message;" << endl
         << "  -c <N>            number of cells;" << endl
         << "  -m <N>            number of materials;" << endl
         << "  -n <N>            number of interations;" << endl
         << "  -V                use verbose output;" << endl
         << "  -v                display version information." << endl; }

void outputVersion()
{   cout << progName << " " << progVers << endl; }

/*  parseInput(const int argc, char** argv)
 *  
 *  Interpret the command line input.
 */
void parseInput(const int argc, char** argv)
{   strcpy(progName, "MultiMatTest");
    strcpy(progVers, PACKAGE_VERSION);
    
    //	Reconstruct command line argument as a string.
    char progCL[256];       //  Complete program command line.
    strcpy(progCL, argv[0]);
    for (int i = 1; i < argc; i++)
    {   strcat(progCL, " ");
        strcat(progCL, argv[i]); }
    
    //  Set variables to defaults, which may be overridden by CLI.
    //verbose            = false;
    
    char   *val;
    if (argc > 1)
    {   int i = 1;
        val = strtok(argv[i++], " ,.-");
        while (val != NULL)
        {   switch (val[0])
            {
                case 'd':   //  detailed memory output.
                    memory_verbose = true;
                    break;
                    
                case 'h':   //  Output help.
                    outputHelp();
                    cout.flush();
                    exit(EXIT_SUCCESS);
                    break;
                    
                case 'c':   //  Number of cells.
                    val = strtok(argv[i++], " ,");
                    ncells = atoi(val);
                    break;
                    
                case 'm':   //  Number of materials.
                    val = strtok(argv[i++], " ,");
                    nmats = atoi(val);
                    break;
                    
                case 'n':   //  Number of iterations.
                    val = strtok(argv[i++], " ,");
                    itermax = atoi(val);
                    break;
                    
                case 'V':   //  Verbose output desired.
                    verbose = true;
                    break;
                    
                case 'v':   //  Version.
                    outputVersion();
                    cout.flush();
                    exit(EXIT_SUCCESS);
                    break;
                    
                default:    //  Unknown parameter encountered.
                    cout << "⚠ Unknown input parameter " << val << endl;
                    outputHelp();
                    cout.flush();
                    exit(EXIT_FAILURE);
                    break; }
            
            val = strtok(argv[i++], " ,.-"); } } }
