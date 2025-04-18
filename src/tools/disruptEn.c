#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>

/*
  disruptEn -- version 1.1
 * 2010-03-17 Revision by Adam Wolenc
 * 1. Provide pardir to findfe77 to match features of sfold and enable
 * 2. data param files opened once only (Ivan Auger)
 * 3. Correct calloc of sampleFE_a and sampleFE_b which gave rise to segfault
 *when #structures exceeded #samples. Resolved segfault.
 */

/*
  disruptEn -- version 1
  01/05/05
  modified from a program of yu shao (Nov 12, 2004)
  Goal: Calculate the disruption energies of RNAi binding sites. 
        Features include:
        1. Parameters: 
           a. RNAi binding site positions (start and end)
           b. use fe.out if present
           c. number of structures (sampleSize)
           d. a switch for different types of disruption energy
              calculation.
        2. different types of disruption energy calculation
           a. whole binding site disruption energy (basic)
           b. a walk through of partial disruption energy from both ends
              of the binding sites.
           c. partial disruption energy for all the possible segments within
              each binding site.
 */

/* 
  Updated 12/29/2004:
  1. use getopt function to get the command line options
  2. converted all the file names to their absolute version
     inside the program (needed to run on cluster) within
     the program, so that relative path names can be used as
     the command line parameters now.
  3. chdir within the program to the sfold bin dir (defined
     by environment varialbe SFOLDBIN -- now the program can
     run from locations.
  4. added an option for the output file (optional), if not specified,
     the output will go to stdout.
 */

#define DEF_FORMAT          1

#ifdef _CRAY
#include <fortran.h>
#define f77findfe         F77FINDFE
#else
#if !defined(_AIX) && !defined(__hpux)
#define f77findfe       f77findfe_
#endif
#define _fcd          char *
#define _cptofcd(a,b) (a)
#define _fcdlen(a)    strlen(a)
#endif

struct ProgArgs {
    int format;
    char *strfile;
    char *seqfile;
};

// constants 
const int maxRNAiLen = 70; // max. length of binding sites on mRNAs

/* function prototype */
int check_fbp(int n, int *fbp);
//2010-02-16 (Adam): This is supposed to provide , char *pardir, int *pardirlen as well.
void f77findfe(int *n, int *numseq, int *fbp, int *h5, int *h3, int *hnum,
        char *seq, double *eall, char *pardir, int *pardirlen);
FILE * file_open(char *p, char *m);
int fillvar_bp(struct ProgArgs *args, int n, char *s, int *fbp);
int fillvar_ct(struct ProgArgs *args, int n, char *s, int *fbp);
int fillvar_hx(struct ProgArgs *args, int n, char *s, int *fbp);
int getctlen(char *fn);
int getseqlen(char *fn);
int readFASTA(char *fn, int n, char *s);
struct ProgArgs * setArgs(int argc, char *argv[]);

int readSampleBP(char * bpf, int mRNALen, int * sampleBP, int sampleSize);
int sampleBPisValid(int mRNALen, int * sampleBP, int sampleSize);
int totalStruct(char * fn);

/* Take command line arguments from C and pass them to the fortran subroutine */
int main(int argc, char *argv[]) {


    int i, j, count, coop,
            sampleSize, *siRNA, *RNAiLen, k, ii, jj, type, totalsiRNA;

    int mRNALen, *numseq, *sampleBP, *tempBP, *h5, *h3, *hnum, pardirlen;

    char mRNAf[10240], mRNAbpf[10240], siRNAf[10240],
            * mRNA, *ch, buf[10240], *tok, fef[10240], *tok1,
            * tok2, outf[10240], pardir[10240];

    FILE *siRNAfh, *fefh, *outfh;

    double eall, *sampleFE_a, *sampleFE_b,
            fe_b, fe_b_ave, fe_a, *fe_allseg;

    int sampleFE_a_ubound, sampleFE_b_ubound, fe_allseg_ubound;

    /* default settings */
    strcpy(mRNAf, "");
    strcpy(siRNAf, "");
    strcpy(mRNAbpf, "");
    coop = 0;
    sampleSize = 1000;
    type = 1;
    strcpy(fef, "");
    strcpy(outf, "");


    while ((i = getopt(argc, argv, "m:p:s:c:f:n:t:o:")) != EOF) {
        switch (i) {
            case 'm':
                strcpy(mRNAf, optarg);
                break;
            case 'p':
                strcpy(siRNAf, optarg);
                break;
            case 's':
                strcpy(mRNAbpf, optarg);
                break;
            case 'c':
                coop = atoi(optarg);
                break;
            case 'f':
                strcpy(fef, optarg);
                break;
            case 'n':
                sampleSize = atoi(optarg);
                break;
            case 'o':
                strcpy(outf, optarg);
                break;
            case 't':
                type = atoi(optarg);
                break;
            case '?':
                exit(-1);
                break;
        }
    }

    if (strcmp(mRNAf, "") == 0 ||
            strcmp(siRNAf, "") == 0 ||
            strcmp(mRNAbpf, "") == 0 ||
            ((coop == 1) && (type != 1))) {

        fprintf(stdout,
                "disruptEn: calculate disruption energy for given RNAi sites.\n");
        fprintf(stdout, "Options:\n");
        fprintf(stdout, "    -m mRNA file name\n");
        fprintf(stdout, "    -p binding site position list file name\n");
        fprintf(stdout, "    -s sample structure file name\n");
        fprintf(stdout, "    -c cooperatively openness (default: 0)\n");
        fprintf(stdout, "       0: open each individual site\n");
        fprintf(stdout, "       1: open all the sites at the same times (only for type 1)\n");
        fprintf(stdout,
                "    -f sample structure free energy file fe.out (defualt: "")\n");
        fprintf(stdout, "    -n sample size (default: 1000)\n");
        fprintf(stdout, "    -o output file (defualt: "")\n");
        fprintf(stdout,
                "    -t type of disruption energy calculation (default: 1)\n");
        fprintf(stdout, "       1: whole site disruption energy (one number)\n");
        fprintf(stdout,
                "       2: partial disruption energy for all possible sub-sites\n");
        fprintf(stdout, "Version 1.0, Dang Long and Yu Shao, JAN 28 2005\n");
        fprintf(stdout, "Version 1.1, Adam Wolenc and Ivan Auger, MAR 17 2010\n");
        return -1;

    }

    /* convert all filenames to absolute file names */
    if (strcmp(mRNAf, "")) {
        realpath(mRNAf, buf);
        strcpy(mRNAf, buf);
    }

    if (strcmp(siRNAf, "")) {
        realpath(siRNAf, buf);
        strcpy(siRNAf, buf);
    }

    if (strcmp(mRNAbpf, "")) {
        realpath(mRNAbpf, buf);
        strcpy(mRNAbpf, buf);
    }

    if (strcmp(fef, "")) {
        realpath(fef, buf);
        strcpy(fef, buf);
    }

    /* for output file, create it if it doesn't exist yet -- since
       realpath only works for existing files. */
    if (strcmp(outf, "")) {
        if (!(outfh = fopen(outf, "w"))) {
            fprintf(stderr, "can't create output file: %s.\n", outf);
            exit(-1);
        }
        fclose(outfh);

        realpath(outf, buf);
        strcpy(outf, buf);
    }


    /*
       -- Change working directory to sfold bin so the program
          can run from anywhere.
       -- Define the sfold bin directory as an environment variable
          (SFOLDBIN) to avoid the hardwired directory name.
     */
    if (!(ch = getenv("SFOLDBIN"))) {
        fprintf(stderr, "Error: environment variable SFOLDBIN not defined.\n");
        exit(-1);
    }
    if (chdir(ch) < 0) {
        fprintf(stderr, "chdir to %s ($SFOLDBIN) failed.\n", ch);
        exit(-1);
    }
    /* else {
      fprintf (stderr, "chdir to %s ($SFOLDBIN) succeeded.\n", ch);
    }
     */

    strncpy(pardir, getenv("SFOLDBIN"), 500);
    strcat(pardir, "../param");
    pardirlen = strlen(pardir);

    mRNALen = getseqlen(mRNAf);

    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    /*   the memory size should include the end of string  */
    /*   character \0 which is used by readFASTA function. */
    /*                      OR                             */
    /*   DO NOT add \0 to the end of the arry -- just like */
    /*   what the original code did.....                   */
    /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
    mRNA = (char *) calloc(((size_t) mRNALen), sizeof (char));
    if (mRNALen != readFASTA(mRNAf, mRNALen, mRNA)) {
        fprintf(stderr, "inconsistent mRNA length: %d : %d\n", mRNALen,
                readFASTA(mRNAf, mRNALen, mRNA));
        return -1;
    };

    /* samplesize is specified through the command line */
    /*
      Before allocating memory, find out # of structure available N first,
      if N < sampleSize, set sampleSize = N and print out a warning message.
     */
    count = totalStruct(mRNAbpf);
    if (count < sampleSize) {
        sampleSize = count;
        fprintf(stderr, "Warning: sampleSize has been changed to %d", count);
        fprintf(stderr, "(only %d structure available in %s.\n", count, mRNAbpf);
    }

    /* read in the positions of the RNAi's */
    /* -- scan the file twice to get the total number of RNAi's first */
    if ((siRNAfh = file_open(siRNAf, "r")) == NULL) {
        fprintf(stderr, "%s: file open failed.\n", siRNAf);
        exit(-1);
    }
    i = 0;
    while ((ch = fgets(buf, 1024, siRNAfh)) != NULL) {
        tok = strtok(ch, " ");
        if (tok == NULL) {
            fprintf(stderr, "Wrong siRNA summary file format!\n");
            exit(-1);
        }
        i++;
    }
    totalsiRNA = i;
    if ((siRNA = calloc(((size_t) totalsiRNA), sizeof (int))) == NULL) {
        fprintf(stderr, "memory allocation failed for RNAi sites.\n");
        exit(-1);
    }
    if ((RNAiLen = calloc(((size_t) totalsiRNA), sizeof (int))) == NULL) {
        fprintf(stderr, "memory allocation failed for RNAi site lengths.\n");
        exit(-1);
    }
    rewind(siRNAfh);
    i = 0;
    while ((ch = fgets(buf, 1024, siRNAfh)) != NULL) {
        tok1 = strtok(ch, "\t");
        siRNA[i] = atoi(tok1);
        tok2 = strtok(NULL, "\n");
        RNAiLen[i] = atoi(tok2) - siRNA[i] + 1;
        i++;
    }
    fclose(siRNAfh);

    /* allocate space */
    numseq = (int *) calloc(((size_t) mRNALen), sizeof (int));
    h5 = (int *) calloc(((size_t) mRNALen), sizeof (int));
    h3 = (int *) calloc(((size_t) mRNALen), sizeof (int));
    hnum = (int *) calloc(((size_t) mRNALen), sizeof (int));
    if (numseq == NULL || h5 == NULL ||
            h3 == NULL || hnum == NULL) {
        fprintf(stderr, "Failed to allocate space for 1-d integer arrays\n");
        return -2;
    }

    fe_allseg_ubound=((maxRNAiLen + 1) * maxRNAiLen / 2)-1;
    fe_allseg = (double *) calloc(((size_t) (maxRNAiLen + 1) * maxRNAiLen / 2),
            sizeof (double));

    /* n -- length of mRNA
       numseq --
       fbp -- structure in base-pair format
       h5 --
       h3 --
       hnum --
       seq -- mRNA sequence
     */

    /* allocate memory for free energy */
    sampleFE_b_ubound=sampleSize-1;
    if ((sampleFE_b = calloc(((size_t) sampleSize), sizeof (double))) == NULL) {
        fprintf(stderr, "Failed to allocate space for 1-d double arrays\n");
        exit(-1);
    }

    //2010-03-17 (Adam): sampleFE_a should be allocated depending on totalsiRNA
    // NOT sampleSize. Segfaults if totalsiRNA>sampleSize.
    sampleFE_a_ubound=totalsiRNA-1;
    if ((sampleFE_a = calloc(((size_t) totalsiRNA), sizeof (double))) == NULL) {
        fprintf(stderr, "Failed to allocate space for 1-d double arrays\n");
        exit(-1);
    }

    //2010-03-16 (Adam): change to single-stage contiguous array allocation
    if ((sampleBP = (int *) calloc((size_t) (sampleSize * mRNALen), sizeof (int)))
            == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return -1;
    }

    if ((tempBP = (int *) calloc(((size_t) mRNALen), sizeof (int))) == NULL) {
        fprintf(stderr, "Memory allocation failed.\n");
        return -1;
    }

    count = readSampleBP(mRNAbpf, mRNALen, sampleBP, sampleSize);
    if (count != (sampleSize - 1)) {
        fprintf(stderr, "incorrect number of structures (%d) read from %s!\n",
                sampleSize, mRNAbpf);
        return -1;
    }

    /* read/calcuate the free energy for the sample structures */
    fe_b_ave = 0;

    if (fefh = fopen(fef, "r")) {
        i = 0;
        while (fgets(buf, 1024, fefh)) {
            tok1 = strtok(buf, " ");
            tok2 = strtok(NULL, " ");
            if (atoi(tok1) - 1 > sampleFE_b_ubound) printf("j=%d out of bounds %s:%d\n", atoi(tok1) - 1, __FILE__, __LINE__);
            sampleFE_b[atoi(tok1) - 1] = atof(tok2);
            fe_b_ave += atof(tok2);
            i++;
            if (i >= sampleSize) break;
        }
        fclose(fefh);
    } else {
        for (i = 0; i < sampleSize; i++) {
            f77findfe(&mRNALen, numseq, &(sampleBP[i * mRNALen]), h5, h3, hnum, mRNA, &fe_b, pardir, &pardirlen);
            if (i > sampleFE_b_ubound) printf("i=%d out of bounds %s:%d\n", i, __FILE__, __LINE__);
            sampleFE_b [i] = fe_b;
            fe_b_ave += fe_b;
            pardirlen = 0;
        }
    }

    fe_b_ave /= sampleSize;

    outfh = fopen(outf, "w");
    if (!outfh) {
        outfh = stdout;
    }

    if (coop == 0) // calc. for individual sites
    {
        for (j = 0; j < totalsiRNA; j++) {
            /* calculate the average energy of the mRNA after
               opening the binding site */
            fe_a = 0;

            for (ii = 0; ii < (RNAiLen[j] + 1) * RNAiLen[j] / 2; ii++) {
                if (ii > fe_allseg_ubound) printf("ii=%d out of bounds %s:%d\n", ii, __FILE__, __LINE__);
                fe_allseg[ii] = 0;
            }

            for (i = 0; i < sampleSize; i++) {
                if (type == 2) {
                    /* disrupt all the possible segments within the binding
                       site of this siRNA
                     */
                    /*
                       siRNA[j] -- starts from 1
                       sampleBP[i*mRNALen+k] -- starts from 1
                       tempBP[k] -- starts from 1
                     */
                    count = 0;
                    for (ii = siRNA[j]; ii < siRNA[j] + RNAiLen[j]; ii++) {
                        for (jj = ii; jj < siRNA[j] + RNAiLen[j]; jj++) {
                            /* for the current segment:
                               start -- ii (starts from ii)
                               end   -- jj (starts from jj)
                             */
                            for (k = 0; k < mRNALen; k++) {
                                if ((k >= (ii - 1) && k <= (jj - 1)) ||
                                        (sampleBP[i * mRNALen + k] >= ii && sampleBP[i * mRNALen + k] <= jj)) {
                                    tempBP[k] = 0;
                                } else tempBP[k] = sampleBP[i * mRNALen + k];
                            }
                            f77findfe(&mRNALen, numseq, tempBP, h5, h3, hnum, mRNA, &eall, pardir, &pardirlen);
                            pardirlen = 0;
                            if (count > fe_allseg_ubound) printf("count=%d out of bounds %s:%d\n", count, __FILE__, __LINE__);
                            fe_allseg[count] += eall;
                            count++;
                        }
                    }
                } // end if type ==2

                /* disrupt the whole siRNA region type =1*/
                for (k = 0; k < mRNALen; k++) {
                    /* remove the base pairs from the siRNA binding site */
                    /* the new structure is saved in tempBP */
                    if ((k >= (siRNA[j] - 1) && k <= (siRNA[j] + RNAiLen[j] - 1)) ||
                            (sampleBP[i * mRNALen + k] >= siRNA[j] &&
                            (sampleBP[i * mRNALen + k] <= siRNA[j] + RNAiLen[j] - 1))) {
                        tempBP[k] = 0;
                    } else {
                        tempBP[k] = sampleBP[i * mRNALen + k];
                        //printf("%3d ", tempBP[k]);
                    }

                } // end for k

                f77findfe(&mRNALen, numseq, tempBP, h5, h3, hnum, mRNA, &eall, pardir, &pardirlen);
                pardirlen = 0;
                fe_a += eall;
            } // end for i from 0 to sampleSize

            fe_a /= sampleSize;
            //if (!sampleBPisValid(mRNALen, sampleBP, sampleSize)) printf("j=%d %s:%d\n", j, __FILE__, __LINE__);
            if (j > sampleFE_a_ubound) printf("j=%d out of bounds %s:%d\n", j, __FILE__, __LINE__);
            sampleFE_a[j] = fe_a;   //manifestation of 2010-03-17 segfault
            //if (!sampleBPisValid(mRNALen, sampleBP, sampleSize)) printf("j=%d %s:%d\n", j, __FILE__, __LINE__);
            fprintf(outfh, "%8.3lf  ", fe_b_ave); // for testing
            
            if (type == 2) {
                for (ii = 0; ii < (RNAiLen[j] + 1) * RNAiLen[j] / 2; ii++) {
                    if (ii > fe_allseg_ubound) printf("ii=%d out of bounds %s:%d\n", ii, __FILE__, __LINE__);
                    fe_allseg[ii] /= sampleSize;
                }
                for (ii = 0; ii < (RNAiLen[j] + 1) * RNAiLen[j] / 2; ii++) {
                    fprintf(outfh, "%8.3lf  ", fe_allseg[ii] - fe_b_ave);
                }
            }

            fprintf(outfh, "%8.3lf  ", fe_a - fe_b_ave);
            fprintf(outfh, "\n");
        } //end for j from 0 to totalsiRNA
    } else
        // open all the sites at one time, only for type1 (whole site disruption ener.)
    {
        fe_a = 0;
        for (i = 0; i < sampleSize; i++) {
            for (k = 0; k < mRNALen; k++) {
                tempBP[k] = sampleBP[i * mRNALen + k];
            }
            for (j = 0; j < totalsiRNA; j++) {
                /* disrupt the whole siRNA regions */
                /* the new structure is saved in tempBP */
                for (k = 0; k < mRNALen; k++) {
                    /* remove the base pairs from the siRNA binding site */
                    if ((k >= (siRNA[j] - 1) && k <= (siRNA[j] + RNAiLen[j] - 1)) || (sampleBP[i * mRNALen + k] >= siRNA[j] &&
                            (sampleBP[i * mRNALen + k] <= siRNA[j] + RNAiLen[j] - 1))) {
                        tempBP[k] = 0;
                    } else {
                        tempBP[k] = tempBP[k];
                    }
                } // end for k
            } //end for j from 0 to totalsiRNA
            f77findfe(&mRNALen, numseq, tempBP, h5, h3, hnum, mRNA, &eall, pardir, &pardirlen);
            pardirlen = 0;
            fe_a += eall;
        } // end for i from 0 to sampleSize
        fe_a /= sampleSize;
        fprintf(outfh, "%8.3lf  ", fe_b_ave); // for testingg

        fprintf(outfh, "%8.3lf  ", fe_a - fe_b_ave);
        fprintf(outfh, "\n");
    }

    fclose(outfh);



    return 0;
}

/*
 * Check the array fbp to make sure it has neither inconsistent base
 * pairs nor pseudoknots
 */
int check_fbp(int n, int *fbp) {
    int i, j;

    for (i = 0; i < n; i++) {
        /* Check for invalid fbp value */
        if (fbp[i] < 0 || fbp[i] > n) {
            fprintf(stderr, "Invalid base pair (%d, %d)\n", i + 1, fbp[i]);
            return -1;
        }

        /* check for consistent pairs */
        if (fbp[i] > 0 && fbp[i] <= n) {
            if (fbp[fbp[i] - 1] != i + 1) {
                fprintf(stderr, "Inconsistent base pair (%d, %d) and (%d, %d)\n",
                        i + 1, fbp[i], fbp[i], fbp[fbp[i] - 1]);
                return -2;
            }

            /* check for pseudoknots */
            if ((i + 1) < fbp[i])
                for (j = i + 1; j < fbp[i] - 1; j++)
                    if ((fbp[j] < i + 1 || fbp[j] > fbp[i]) && fbp[j] != 0) {
                        fprintf(stderr, "Found pseudoknot (%d, %d) x (%d, %d)\n",
                                i + 1, fbp[i], j + 1, fbp[j]);
                        return -3;
                    }
        }
    }

    return 0;
}

/*
 * A wrapper function for fopen to determine where we should read from.
 * If STDIN is the source of input, we can call fdopen instead of fopen
 */
FILE *file_open(char *p, char *m) {
    FILE *fp;

    if (strcasecmp(p, "STDIN") == 0) {
        /* we are reading from the standard input (0) */
        fp = fdopen(0, m);
    } else {
        /* do a normal file open */
        fp = fopen(p, m);
    }

    return fp;
}

/* read "sampleSize" number of sample structures from the bp.out */
int readSampleBP(char * bpf, int mRNALen, int * sampleBP, int sampleSize) {
    int count;
    FILE *fp;
    char *ch, *tok;
    char buf[256];
    int i, j;

    if ((fp = file_open(bpf, "r")) == NULL) {
        fprintf(stderr, "Unable to read from file '%s'\n", bpf);
        return -1;
    }

    /*  data start here */
    count = -1;
    while ((ch = fgets(buf, 256, fp)) != NULL) {
        i = j = 0;

        /* if (strncmp (ch, " Structure", 10) == 0) { */
        if (strstr(ch, "Structure") != NULL) {
            /* if (count >= 9) break; */
            count++;
            /* fprintf (stderr, "structure: %d\n", count); */
            if (count >= sampleSize) {
                fprintf(stderr,
                        "Warning: more sample structures than allocated (%d)!\n",
                        sampleSize);
                /* exit (-1); */
                return count - 1;
            }
            continue;
        }
        /* 1st column: position i */
        tok = strtok(ch, " ");
        if (tok != NULL)
            i = atoi(tok);
        /* 2nd column: position j */
        tok = strtok(NULL, " ");
        if (tok != NULL)
            j = atoi(tok);
        /* Check to see if it has the 3rd column; if so, possibly a different format */
        tok = strtok(NULL, " ");
        if (tok != NULL) {
            fprintf(stderr, "Structure file not in bp format\n");
            return -3;
        }

        if (i <= 0 || j <= 0 || i > mRNALen || j > mRNALen) {
            fprintf(stderr, "inappropriate base pair: %d, %d (%d),%s\n", i, j, mRNALen, ch);
            exit(-1);
        }

        if (i > 0 && i <= mRNALen && j > 0 && j <= mRNALen && i != j) {
            sampleBP[count * mRNALen + i - 1] = j;
            sampleBP[count * mRNALen + j - 1] = i;
        }

    }

    fclose(fp);

    return count;
}

int sampleBPisValid(int mRNALen, int * sampleBP, int sampleSize) {
    // 2010-03-17 (Adam): Integrity Check
    int i;
    for (i = 0; i < sampleSize * mRNALen; i++) {
        if (sampleBP[i] > mRNALen || sampleBP[i] < 0) {
            printf("Invalid value in sampleBP[%d]=%d\n", i, sampleBP[i]);
            return 0;
        }
    }
    return 1;
}

/*
 * Fill sequence data and structure from bp file into arrays
 */
int fillvar_bp(struct ProgArgs *args, int n, char *s, int *fbp) {
    int count = 0;
    FILE *fp;
    char *ch, *tok;
    char buf[256];
    int i, j;

    count = readFASTA(args->seqfile, n, s);
    if (count != n) {
        fprintf(stderr, "Inconsistent sequence length\n");
        return -2;
    }

    if ((fp = file_open(args->strfile, "r")) == NULL) {
        fprintf(stderr, "Unable to read from file '%s'\n", args->strfile);
        return -1;
    }

    /*  data start here */
    do {
        i = j = 0;

        ch = fgets(buf, 256, fp);
        /* 1st column: position i */
        tok = strtok(ch, " ");
        if (tok != NULL)
            i = atoi(tok);
        /* 2nd column: position j */
        tok = strtok(NULL, " ");
        if (tok != NULL)
            j = atoi(tok);
        /* Check to see if it has the 3rd column; if so, possibly a different format */
        tok = strtok(NULL, " ");
        if (tok != NULL) {
            fprintf(stderr, "Structure file not in bp format\n");
            return -3;
        }

        if (i > 0 && i <= n && j > 0 && j <= n && i != j) {
            fbp[i - 1] = j;
            fbp[j - 1] = i;
        }

    } while (ch != NULL);

    fclose(fp);

    return count;
}

/*
 * Fill sequence data and structure from GCG connect file into arrays
 */
int fillvar_ct(struct ProgArgs *args, int n, char *s, int *fbp) {
    int count = 0;
    FILE *fp;
    char buf[256];
    char *ch, *tok;
    int i, j;
    char nuc;

    if ((fp = file_open(args->strfile, "r")) == NULL) {
        fprintf(stderr, "Unable to read from file '%s'\n", args->strfile);
        return -1;
    }

    /* sequence data start here */
    do {
        i = j = 0;
        nuc = ' ';

        ch = fgets(buf, 256, fp);
        /* 1st column: position i */
        tok = strtok(ch, " ");
        if (tok != NULL)
            i = atoi(tok);
        /* 2nd column: Sequence alphabet */
        tok = strtok(NULL, " ");
        if (tok != NULL)
            nuc = tok[0];
        /* Skip column 3 and 4 */
        tok = strtok(NULL, " ");
        tok = strtok(NULL, " ");
        /* 5th column: position j */
        tok = strtok(NULL, " ");
        if (tok != NULL)
            j = atoi(tok);

        if (i > 0 && i <= n && j >= 0 && j <= n && i != j &&
                (nuc == 'A' || nuc == 'T' || nuc == 'C' || nuc == 'G' || nuc == 'U' ||
                nuc == 'a' || nuc == 't' || nuc == 'c' || nuc == 'g' || nuc == 'u')) {
            fbp[i - 1] = j;
            fbp[j - 1] = i;
            nuc = toupper(nuc);
            if (nuc == 'T')
                nuc = 'U';
            s[count++] = nuc;
        }
    } while (ch != NULL && count < n);

    fclose(fp);

    if (count != n) {
        fprintf(stderr, "Inconsistent sequence length\n");
        return -2;
    }

    return count;
}

/*
 * Fill sequence data and structure from helix format file into arrays
 */
int fillvar_hx(struct ProgArgs *args, int n, char *s, int *fbp) {
    int count = 0;
    FILE *fp;
    char *ch, *tok;
    char buf[256];
    int i, j, k, x;

    count = readFASTA(args->seqfile, n, s);
    if (count != n) {
        fprintf(stderr, "Inconsistent sequence length\n");
        return -2;
    }

    if ((fp = file_open(args->strfile, "r")) == NULL) {
        fprintf(stderr, "Unable to read from file '%s'\n", args->strfile);
        return -1;
    }

    /*  data start here */
    do {
        i = j = k = 0;

        ch = fgets(buf, 256, fp);
        /* 1st column: position i */
        tok = strtok(ch, " ");
        if (tok != NULL)
            i = atoi(tok);
        /* 2nd column: position j */
        tok = strtok(NULL, " ");
        if (tok != NULL)
            j = atoi(tok);
        /* 3rd column: count k */
        tok = strtok(NULL, " ");
        if (tok != NULL)
            k = atoi(tok);
        if (i != 0 && j != 0 && k == 0) {
            fprintf(stderr, "Structure file not in helix format\n");
            return -3;
        }

        if (i > 0 && i <= n && j > 0 && j <= n && k > 0 && k <= n && i != j && k <= (n - 3) / 2) {
            for (x = 0; x < k; x++) {
                fbp[i + x - 1] = j - x;
                fbp[j - x - 1] = i + x;
            }
        }

    } while (ch != NULL);

    fclose(fp);

    return count;
}

/*
 * Parse the GCG connect file and find the sequence length
 */
int getctlen(char *fn) {
    int count = 0, i, j;
    FILE *fp;
    char buf[256];
    char *ch, *tok;
    char nuc;

    if ((fp = file_open(fn, "r")) == NULL) {
        fprintf(stderr, "Unable to read from file '%s'\n", fn);
        return -1;
    }

    /* sequence data start here */
    do {
        i = j = 0;
        nuc = ' ';

        ch = fgets(buf, 256, fp);
        /* 1st column: position i */
        tok = strtok(ch, " ");
        if (tok != NULL)
            i = atoi(tok);
        /* 2nd column: Sequence alphabet */
        tok = strtok(NULL, " ");
        if (tok != NULL)
            nuc = tok[0];
        /* Skip column 3 and 4 */
        tok = strtok(NULL, " ");
        tok = strtok(NULL, " ");
        /* 5th column: position j */
        tok = strtok(NULL, " ");
        if (tok != NULL)
            j = atoi(tok);

        if (i > 0 && j >= 0 && i != j &&
                (nuc == 'A' || nuc == 'T' || nuc == 'C' || nuc == 'G' || nuc == 'U' ||
                nuc == 'a' || nuc == 't' || nuc == 'c' || nuc == 'g' || nuc == 'u')) {
            count++;
        }
    } while (ch != NULL);

    fclose(fp);

    if (count == 0)
        fprintf(stderr, "Sequence is of zero length\n");

    return count;
}

/*
 * Parse the input sequence file (in FASTA format) to find the sequence length
 */
int getseqlen(char *fn) {
    int count = 0, ch;
    FILE *fp;

    if ((fp = fopen(fn, "r")) == NULL) {
        fprintf(stderr, "Unable to read from file '%s'\n", fn);
        return -1;
    }

    ch = fgetc(fp);
    if (ch != '>') {
        fprintf(stderr, "Invalid FASTA file format\n");
        return -2;
    }

    /* skip sequence header */
    do {
        ch = fgetc(fp);
    } while (ch != EOF && ch != '\n');

        /* sequence data start here */
        while (ch != EOF) {
            ch = fgetc(fp);
            if (ch == 'A' || ch == 'T' || ch == 'C' || ch == 'G' || ch == 'U' ||
                    ch == 'a' || ch == 't' || ch == 'c' || ch == 'g' || ch == 'u' ||
                    ch == 'N' || ch == 'n')
                count++;
        }

    fclose(fp);

    if (count == 0)
        fprintf(stderr, "Sequence is of zero length\n");

    return count;
}

/*
 * Read FASTA file and store sequence data into an array
 */
int readFASTA(char *fn, int n, char *s) {
    int count = 0, ch;
    FILE *fp;

    if ((fp = fopen(fn, "r")) == NULL) {
        fprintf(stderr, "Unable to read from file '%s'\n", fn);
        return -1;
    }

    ch = fgetc(fp);
    if (ch != '>') {
        fprintf(stderr, "Invalid FASTA file format\n");
        return -2;
    }

    /* skip sequence header */
    do {
        ch = fgetc(fp);
    } while (ch != EOF && ch != '\n');

        /* sequence data start here */
        while (ch != EOF && count < n) {
            ch = fgetc(fp);
            if (ch == 'A' || ch == 'T' || ch == 'C' || ch == 'G' || ch == 'U' ||
                    ch == 'a' || ch == 't' || ch == 'c' || ch == 'g' || ch == 'u' ||
                    ch == 'N' || ch == 'n') {
                ch = toupper(ch);
                if (ch == 'T')
                    ch = 'U';
                s[count++] = ch;
            }
        }

    fclose(fp);

    if (count == 0)
        fprintf(stderr, "Sequence is of zero length\n");

    /* this bug cost me 2 hours... */
    /*
    s[count] = '\0';
     */

    return count;
}

struct ProgArgs * setArgs(int argc, char *argv[]) {
    struct ProgArgs *args;
    int i, itmp;

    args = (struct ProgArgs *) malloc(sizeof (struct ProgArgs));
    if (args == NULL) {
        fprintf(stderr, "Failed to allocate memory for program arguments\n");
        return NULL;
    }

    /* Initialize structure elements */
    args->format = DEF_FORMAT;
    args->strfile = NULL;
    args->seqfile = NULL;

    /* Do the argument matching here... */
    if (argc < 2) {
        return NULL;
    }

    for (i = 1; i < argc; i++)
        if (strcmp(argv[i], "-h") == 0)
            return NULL;

    i = 1;
    while (i < argc) {
        if (argv[i][0] == '-') {
            /* this is an option field */
            if (i >= argc - 1) {
                fprintf(stderr, "Argument to '%s' is missing\n", argv[i]);
                return NULL;
            }

            if (strcmp(argv[i] + 1, "f") == 0) {
                itmp = atoi(argv[++i]);
                if (itmp < 1 || itmp > 3) {
                    fprintf(stderr, "Unknown structure file format %d\n", itmp);
                    return NULL;
                } else args->format = itmp;
            } else if (strcmp(argv[i] + 1, "s") == 0) {
                if (!args->seqfile)
                    args->seqfile = argv[++i];
                else
                    fprintf(stderr, "Ignoring duplicate sequence file '%s'\n", argv[++i]);
            } else {
                fprintf(stderr, "Unknown argument '%s'\n", argv[i]);
                return NULL;
            }

        } else {
            /* this specifies an input structure file */
            if (!args->strfile)
                args->strfile = argv[i];
            else
                fprintf(stderr, "Ignoring duplicate structure file '%s'\n", argv[i]);
        }

        i++;
    }

    if (!args->strfile) {
        fprintf(stderr, "No structure file given.\n");
        return NULL;
    }

    if (strcasecmp(args->strfile, "STDIN") == 0 && args->format == 1) {
        fprintf(stderr, "Unable to read from standard input for GCG connect format.\n");
        fprintf(stderr, "Please provide a structure file name as input.\n");
        return NULL;
    }

    if (args->format != 1)
        if (!args->seqfile) {
            fprintf(stderr, "Sequence file required for this structure format.\n");
            return NULL;
        }

    return args;
}

int totalStruct(char * fn) {
    FILE * fh;
    char buf[1024];
    int count = 0;

    if (!(fh = fopen(fn, "r"))) {
        fprintf(stderr, "couldn't open file:  %s\n", fn);
        exit(-1);
    }

    while (fgets(buf, 1024, fh)) {
        if (strstr(buf, "Structure")) count++;
    }
    fclose(fh);

    return (count);
}
