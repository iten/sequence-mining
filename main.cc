#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <utility>
#include "protein-scanner.h"
#include "dna-scanner.h"

#define MAX_FILE_SIZE 4*1024*1024 // 4 MB
#define ARG_PRIOR 1000
#define ARG_MIN_DNA_FRACTION 1002
#define ARG_PROTEIN_REGEX 1003
#define ARG_DNA_POSSIBLE_REGEX 1004
#define ARG_MIN_PURE_LEN 1005
#define ARG_MIN_IMPURE_LEN 1006
#define ARG_MIN_OUTPUT_LEN 1007
#define ARG_HTML 1008
#define ARG_THRESHOLD 1009
using std::string;

void usage(char *name) {
    printf("Usage: %s positive-matrix negative-matrix file [file2]...\n"
           "Options:\n"
           " -h, --html: Remove HTML tags from input"
           " -p <prior>, --prior <prior>: Prior probability of seeing a protein\n"
           " --protein_regex <regex>: Regex that matches words that may contain"
           " a protein (needs a capture group)\n"
           " --dna_possible_regex <regex>: Regex that match possible DNA "
           "sequences\n"
           " --min_pure_len <length>: minimum length that a pure DNA fragment "
           "must have to be added to the growing sequence\n"
           " --min_dna_fraction <fraction>: minimum ratio of (DNA content)/"
           "(total length) that a fragment must have to be added to the growing"
           " sequence\n"
           " --min_impure_len <length>: minimum length that an impure DNA "
           "fragment must have in order to be added to the growing sequence\n"
           " --min_output_len <length>: minimum length that a final DNA "
           "sequence must have to be output\n", name);
}

int main(int argc, char *argv[]) {
    // options
    const struct option longopts[] = {
        { "prior", 1, NULL, ARG_PRIOR },
        { "min_dna_fraction", 1, NULL, ARG_MIN_DNA_FRACTION },
        { "protein_regex", 1, NULL, ARG_PROTEIN_REGEX },
        { "dna_possible_regex", 1, NULL, ARG_DNA_POSSIBLE_REGEX },
        { "min_pure_len", 1, NULL, ARG_MIN_PURE_LEN },
        { "min_impure_len", 1, NULL, ARG_MIN_IMPURE_LEN },
        { "min_output_len", 1, NULL, ARG_MIN_OUTPUT_LEN },
        { "html", 0, NULL, ARG_HTML },
        { "threshold", 1, NULL, ARG_THRESHOLD },
        { 0, 0, 0, 0 }
    };
    double prior = 0.01f, min_dna_fraction = 0.7f, threshold = 0.7f;
    const char *prot_regex = "[^A-Za-z]([ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]{10,})[^A-Za-z]";
    const char *dna_possible_regex = "([^- \n]+)";
    size_t min_pure_len = 3, min_impure_len = 10, min_output_len = 15;
    bool is_html = false;
    while(int arg = getopt_long(argc, argv, "p:h", longopts, NULL)) {
        if(arg == -1) {
            break;
        }
        switch(arg) {
        case 'p':
        case ARG_PRIOR:
            // Prior
            sscanf(optarg, "%lf", &prior);
            break;
        case 'h':
        case ARG_HTML:
            // HTML document instead of plaintext
            is_html = true;
            break;
        case ARG_PROTEIN_REGEX:
            prot_regex = optarg;
            break;
        case ARG_DNA_POSSIBLE_REGEX:
            dna_possible_regex = optarg;          
            break;
        case ARG_MIN_PURE_LEN:
            sscanf(optarg, "%lu", &min_pure_len);         
            break;
        case ARG_MIN_IMPURE_LEN:
            sscanf(optarg, "%lu", &min_impure_len);
            break;
        case ARG_MIN_OUTPUT_LEN:
            sscanf(optarg, "%lu", &min_output_len);
            break;
        case ARG_MIN_DNA_FRACTION:
            sscanf(optarg, "%lf", &min_dna_fraction);
            break;
        case ARG_THRESHOLD:
            sscanf(optarg, "%lf", &threshold);
            break;
        default:
            usage(argv[0]);
            exit(1);
        }
    }
    if(argv[optind] == NULL || argv[optind + 1] == NULL
       || argv[optind + 2] == NULL) {
        usage(argv[0]);
        exit(1);
    }
    ProteinScanner prot_scanner(prot_regex, argv[optind], argv[optind + 1],
                                prior, threshold);
    DnaScanner dna_scanner(dna_possible_regex, min_pure_len, min_impure_len,
                           min_output_len, min_dna_fraction);
    for(int i = optind + 2; argv[i] != NULL; ++i) {
        FILE *f = fopen(argv[i], "r");
        if(f == NULL) {
            perror("Couldn't open input file");
        }
        // Get file size
        struct stat stat;
        fstat(fileno(f), &stat);
        size_t size = stat.st_size;
        // Input buffer - NB: maps the entire file into memory -- we're working
        // with small files.
        char *buf = (char *) malloc(size + 1);
        if(buf==NULL) {
            perror("Buffer allocation failed");
        }
        size_t bytes_read = fread(buf, 1, size, f);
        buf[bytes_read] = '\0';
        if(is_html) {
            // GlobalReplace needs a string pointer rather than char **
            // Slow since this will copy all the data
            string replace = buf;
            RE2::GlobalReplace(&replace, "<.*?>", " ");
            // Copy it back into the buffer -- this is a little stupid
            strncpy(buf, replace.c_str(), size);
            buf[size] = '\0'; // Ensure null termination
        }
        printf("======== FILE %s ========\n", argv[i]);
        printf("======== PROTEINS ========\n");
        std::vector< std::pair<string, double> > proteins;
        prot_scanner.Scan(buf, &proteins);
        for(int j = 0; j < proteins.size(); ++j) {
            printf("%s %lf\n", proteins[j].first.data(), proteins[j].second);
        }
        printf("======== DNA ========\n");
        std::vector<string> dna;
        dna_scanner.Scan(buf, &dna);
        for(int j = 0; j < dna.size(); ++j) {
            printf("%s\n", dna[j].data());
        }
        free(buf);
        fclose(f);
    }
    return 0;
}
