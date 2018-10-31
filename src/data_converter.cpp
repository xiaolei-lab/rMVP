#ifdef _OPENMP
#include <omp.h>
#endif

#include <Rcpp.h>
#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/progress.hpp>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <fstream>

using namespace std;
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]

#include <bigmemory/isna.hpp>
#include <bigmemory/MatrixAccessor.hpp>


// ***** VCF *****

// [[Rcpp::export]]
List vcf_parser_map(std::string vcf_file, std::string out) {
    // Define
    const int MAP_INFO_N = 50;      // max length of "SNP, POS and CHROM"
    ifstream file(vcf_file);
    ofstream map(out + ".map");
    ofstream indfile(out + ".geno.ind");
    
    string line;
    vector<string> l;
    vector<string> ind;
    size_t n, m;
    
    // Skip Header
    string prefix("#CHROM");
    while (getline(file, line)) {
        if (!line.compare(0, prefix.size(), prefix)) { break; }
    }
    
    // Get inds
    boost::split(ind, line, boost::is_any_of("\t"));
    vector<string>(ind.begin() + 9, ind.end()).swap(ind);   // delete first 9 columns
    for (int i = 0; i < ind.size(); i++) {
        indfile << ind[i] << endl;
    }
    indfile.close();
    
    
    // Get num of ind / marker
    n = ind.size();
    
    // unsigned m = count(                 // 3-4s per 100Mb
    //     istream_iterator<char>(file),
    //     istream_iterator<char>(),
    //     '\n');
    // file.seekg(pos);
    
    // Map
    map << "SNP\tCHROM\tPOS"<< endl;
    m = 0;
    while (getline(file, line)) {
        string tmp = line.substr(0, MAP_INFO_N);
        boost::split(l, tmp, boost::is_any_of("\t"));
        
        // map
        if (l[2] == ".") {      // snp name missing
            l[2] = l[0] + '-' + l[1];
        }
        map << l[2] << '\t' << l[0] << '\t' << l[1] << endl;
        m++;
    }
    map.close();
    file.close();
    
    return List::create(_["n"] = n,
                        _["m"] = m);
}

template <typename T>
T vcf_marker_parser(string m, double NA_C) {
    if (m[0] != '.' && m[2] != '.') {
        return static_cast<T>(m[0] - '0' + m[2] - '0');
    } else {
        return static_cast<T>(NA_C);
    }
}

template <typename T>
void vcf_parser_genotype(std::string vcf_file, XPtr<BigMatrix> pMat, double NA_C, bool show_progress=true) {
    // define
    ifstream file(vcf_file);
    
    string line;
    vector<string> l;
    vector<char> markers;
    size_t m;
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    // progress bar
    boost::progress_display * progress = NULL;
    if(show_progress) {
        progress = new boost::progress_display(pMat->nrow());
    }
    
    // Skip Header
    string prefix("#CHROM");
    while (getline(file, line)) {
        if (!line.compare(0, prefix.size(), prefix)) { break; }
    }
    
    // parser genotype
    m = 0;
    while (getline(file, line)) {
        boost::split(l, line, boost::is_any_of("\t"));
        vector<string>(l.begin() + 9, l.end()).swap(l);
        markers.clear();
        transform(
            l.begin(), l.end(), 
            back_inserter(markers),
            boost::bind<T>(&vcf_marker_parser<T>, _1, NA_C)
        );
        for (size_t i = 0; i < markers.size(); i++) {
            // Rcout << m << ", " << i << "\t" << markers[i] << endl;
            // mat[i][m] = markers[i] == MISSING ? static_cast<T>(NA_C) : markers[i];
            mat[i][m] = markers[i];
        }
        m++;
        if(show_progress) {
            ++(*progress);
        }
    }
}

// [[Rcpp::export]]
void vcf_parser_genotype(std::string vcf_file, SEXP pBigMat, bool show_progress=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return vcf_parser_genotype<char>(vcf_file, xpMat, NA_CHAR, show_progress);
    case 2:
        return vcf_parser_genotype<short>(vcf_file, xpMat, NA_SHORT, show_progress);
    case 4:
        return vcf_parser_genotype<int>(vcf_file, xpMat, NA_INTEGER, show_progress);
    case 8:
        return vcf_parser_genotype<double>(vcf_file, xpMat, NA_REAL, show_progress);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}


// ***** HAPMAP *****

// [[Rcpp::export]]
List hapmap_parser_map(Rcpp::StringVector hmp_file, std::string out) {
    // Define
    const int MAP_INFO_N = 50;      // max length of "SNP, POS and CHROM"
    ofstream map(out + ".map");
    ofstream indfile(out + ".geno.ind");
    
    string line;
    vector<string> l;
    vector<string> ind;
    size_t n;
    size_t m;
    
    for (int i = 0; i < hmp_file.size(); i++) {
        ifstream file(hmp_file(i));
        
        // Skip Header
        string prefix("rs#");
        while (true) {
            getline(file, line);
            if (!line.compare(0, prefix.size(), prefix)) { break; }
        }
        
        // Get inds
        boost::split(ind, line, boost::is_any_of("\t"));
        vector<string>(ind.begin() + 11, ind.end()).swap(ind);   // delete first 11 columns
        for (int i = 0; i < ind.size(); i++) {
            indfile << ind[i] << endl;
        }
        indfile.close();
        
        // Get num of ind / marker
        n = ind.size();
        
        // Map
        map << "SNP\tCHROM\tPOS"<< endl;
        m = 0;
        while (getline(file, line)) {
            string tmp = line.substr(0, MAP_INFO_N);
            boost::split(l, tmp, boost::is_any_of("\t"));
            
            // map
            if (l[0] == ".") {      // snp name missing
                l[0] = l[2] + '-' + l[3];
            }
            map << l[0] << '\t' << l[2] << '\t' << l[3] << endl;
            m++;
        }
        map.close();
        file.close();
        
    }
    
    
    return List::create(_["n"] = n,
                        _["m"] = m);
}

template <typename T>
T hapmap_marker_parser(string m, char major, double NA_C) {
    int number = 0;
    if (m.length() != 2 || 
        // m[0] == 'N' || m[1] == 'N' ||
        (m[0] != 'A' && m[0] != 'T' && m[0] != 'G' && m[0] != 'C') ||
        (m[1] != 'A' && m[1] != 'T' && m[1] != 'G' && m[1] != 'C')
    ) {
        return static_cast<T>(NA_C);
    } else {
        number += (m[0] == major)?0:1;
        number += (m[1] == major)?0:1;
        return static_cast<T>(number);
    }
}

template <typename T>
void hapmap_parser_genotype(std::string hmp_file, XPtr<BigMatrix> pMat, double NA_C, bool show_progress=true) {
    // define
    ifstream file(hmp_file);
    
    string line;
    char major;
    vector<string> l;
    vector<char> markers;
    size_t m;
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    // progress bar
    boost::progress_display * progress = NULL;
    if(show_progress) {
        progress = new boost::progress_display(pMat->nrow());
    }
    
    // Skip Header
    string prefix("rs#");
    while (true) {
        getline(file, line);
        if (!line.compare(0, prefix.size(), prefix)) { break; }
    }
    
    // parser genotype
    m = 0;
    while (getline(file, line)) {
        boost::split(l, line, boost::is_any_of("\t"));
        major = l[1][0];
        vector<string>(l.begin() + 11, l.end()).swap(l);
        markers.clear();
        transform(
            l.begin(),
            l.end(),
            back_inserter(markers),
            boost::bind<T>(&hapmap_marker_parser<T>, _1, major, NA_C)
        );
        for (size_t i = 0; i < markers.size(); i++) {
            // Rcout << m << ", " << i << "\t" << markers[i] << endl;
            // mat[i][m] = markers[i] == MISSING ? static_cast<T>(NA_C) : markers[i];
            mat[i][m] = markers[i];
        }
        m++;
        if(show_progress) {
            ++(*progress);
        }
    }
}

// [[Rcpp::export]]
void hapmap_parser_genotype(std::string hmp_file, SEXP pBigMat, bool show_progress=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return hapmap_parser_genotype<char>(hmp_file, xpMat, NA_CHAR, show_progress);
    case 2:
        return hapmap_parser_genotype<short>(hmp_file, xpMat, NA_SHORT, show_progress);
    case 4:
        return hapmap_parser_genotype<int>(hmp_file, xpMat, NA_INTEGER, show_progress);
    case 8:
        return hapmap_parser_genotype<double>(hmp_file, xpMat, NA_REAL, show_progress);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}


// ***** NUMERIC *****

// [[Rcpp::export]]
List numeric_scan(std::string num_file) {
    // Define
    string line;
    vector<string> l;
    vector<string> ind;
    size_t n, m;
    
    ifstream file(num_file);
    
    getline(file, line);
    boost::split(l, line, boost::is_any_of("\t ,"));

    n = l.size();
    m = 1;
    while (getline(file, line))
        m++;
    
    return List::create(_["n"] = n,
                        _["m"] = m);
}


// ***** BFILE *****

template <typename T>
void write_bfile(XPtr<BigMatrix> pMat, std::string bed_file, double NA_C, bool show_progress=true) {
    // check input
    string ending = ".bed";
    if (0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending)) {
        bed_file += ending;
    }
    
    // define
    T c;
    long m = pMat->nrow();
    long n = pMat->ncol() / 4;  // 4 individual = 1 bit
    if (pMat->ncol() % 4 != 0) 
        n++;
    
    vector<uint8_t> geno(n);
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    ofstream fout(bed_file, ios::out | ios::binary);
    
    // progress bar
    boost::progress_display * progress = NULL;
    if(show_progress) {
        progress = new boost::progress_display(m);
    }
    
    // magic number of bfile
    const unsigned char magic_bytes[] = { 0x6c, 0x1b, 0x01 };
    fout.write((char*)magic_bytes , 3);
    
    // map
    std::map<T, int> code;
    code[0] = 3;
    code[1] = 2;
    code[2] = 0;
    code[static_cast<T>(NA_C)] = 1;
    
    // write bfile
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            uint8_t p = 0;
            for (size_t x = 0; x < 4 && (4 * j + x) < pMat->ncol(); x++) {
                c = mat[4 * j + x][i];
                if (isna(c)) {
                    p |= 1 << (x*2);
                } else {
                    p |= code[c] << (x*2);
                }
            }
            geno[j] = p;
        }
        fout.write((char*)geno.data(), geno.size());
        if(show_progress) {
            ++(*progress);
        }
    }
    fout.close();
    return;
}

// [[Rcpp::export]]
void write_bfile(SEXP pBigMat, std::string bed_file, bool show_progress=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return write_bfile<char>(xpMat, bed_file, NA_CHAR, show_progress);
    case 2:
        return write_bfile<short>(xpMat, bed_file, NA_SHORT, show_progress);
    case 4:
        return write_bfile<int>(xpMat, bed_file, NA_INTEGER, show_progress);
    case 8:
        return write_bfile<double>(xpMat, bed_file, NA_REAL, show_progress);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

template <typename T>
void read_bfile(std::string bed_file, XPtr<BigMatrix> pMat, long maxLine, double NA_C, int threads=2, bool show_progress=true) {
    // check input
    string ending = ".bed";
    if (0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending))
        bed_file += ending;
    
    // openmp
#ifdef _OPENMP
    if ( threads > 0 )
        omp_set_num_threads( threads );
    Rcerr << "Number of threads=" << omp_get_max_threads() << endl;
#endif
    
    // define
    long n = pMat->ncol() / 4;  // 4 individual = 1 bit
    if (pMat->ncol() % 4 != 0) 
        n++; 
    char * buffer;
    long buffer_size;
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    // map
    std::map<int, T> code;
    code[3] = 0;
    code[2] = 1;
    code[1] = static_cast<T>(NA_C);
    code[0] = 2;
    
    // open file
    ifstream fin(bed_file, ios::in | ios::binary);
    fin.seekg (0, fin.end);
    long length = fin.tellg();
    fin.seekg (0, fin.beg);
    
    // get buffer_size
    if (maxLine <= 0) {
        buffer_size = length - 3;
        show_progress = false;
    } else {    // memory
        buffer_size = maxLine * n;
    }
    
    // progress bar
    int n_block = (length - 3) / buffer_size;
    if ((length - 3) % buffer_size != 0) { n_block++; }
    Progress progress(n_block, show_progress);
    
    // magic number of bfile
    buffer = new char [3];
    fin.read(buffer, 3);
    
    // loop file
#pragma omp parallel for schedule(dynamic) 
    for (int i = 0; i < n_block; i++) {
        buffer = new char [buffer_size];
        fin.read(buffer, buffer_size);
        
        size_t r, c;
        // i: current block start, j: current bit.
        for (size_t j = 0; j < buffer_size && i * buffer_size + j < length - 3; j++) {
            // bit -> item in matrix
            r = (i * buffer_size + j) / n;
            c = (i * buffer_size + j) % n * 4;
            uint8_t p = buffer[j];
            
            for (size_t x = 0; x < 4 && (c + x) < pMat->ncol(); x++) {
                mat[c + x][r] = code[(p >> (2*x)) & 0x03];
            }
        }
        progress.increment();
    }
    fin.close();
    return;
}

// [[Rcpp::export]]
void read_bfile(std::string bed_file, SEXP pBigMat, long maxLine, int threads=2, bool show_progress=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return read_bfile<char>(bed_file, xpMat, maxLine, NA_CHAR, threads, show_progress);
    case 2:
        return read_bfile<short>(bed_file, xpMat, maxLine, NA_SHORT, threads, show_progress);
    case 4:
        return read_bfile<int>(bed_file, xpMat, maxLine, NA_INTEGER, threads, show_progress);
    case 8:
        return read_bfile<double>(bed_file, xpMat, maxLine, NA_REAL, threads, show_progress);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

/*** R
# setwd("~/code/MVP/src")
library(bigmemory)

# # Test case
# # MVP.Data.vcf_parser('/Volumes/Hyacz-WD/Biodata/demo/vcf/maize.vcf', 'maize1')
# MVP.Data.vcf_parser('/Users/hyacz/code/MVP/demo.data/VCF/myVCF.vcf', 'mvp1')
# system.time({MVP.Data.impute('mvp1.geno.desc')})
# x <- attach.big.matrix('mvp1.geno.desc')
# x[1:5,1:5]

c <- big.matrix(3093, 279)
read_bfile('demo.data/bed/plink.bed', c@address, -1)
c[1, 1:20]
library(plink2R)
a <- read_plink('demo.data/bed/plink')
b <- t(as.matrix(a$bed))
c <- as.matrix(c)
b[1, 1:10]
c[1, 1:20]
# numeric_scan('/Users/hyacz/code/MVP/demo.data/numeric/mvp.hmp.Numeric.txt')
*/
