// Data pre-processing module
// 
// Copyright (C) 2016-2018 by Xiaolei Lab
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "rMVP.h"
// #include <boost/bind/bind.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace std;
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(BH, bigmemory)]]
// [[Rcpp::depends(RcppProgress)]]

#include <bigmemory/isna.hpp>
#include <bigmemory/MatrixAccessor.hpp>

vector<string> split_line(const string& str, const string& whitespace = " ,\t\r\n")
{
    vector<string> result;
    size_t col_end = 0;
    while (true) {
        size_t col_begin = str.find_first_not_of(whitespace, col_end);
        // LOG_TRACE_FUNC << "col_begin: " << col_begin << endl;
        if (col_begin == std::string::npos)
            return result;
        col_end = str.find_first_of(whitespace, col_begin);
        // LOG_TRACE_FUNC << "col_end: " << col_end << endl;
        if (col_end == std::string::npos)
            col_end = str.size();
        string stri = str.substr(col_begin, col_end - col_begin);
        stri.erase(stri.find_last_not_of(" \n\r\t") + 1);
        result.push_back(stri);
    }
}

// ***** VCF *****

// [[Rcpp::export]]
List vcf_parser_map(std::string vcf_file, std::string out) {
    // Define
    const int MAP_INFO_N = 10000;      // max length of "SNP, POS and CHROM"
    ifstream file(vcf_file);
    ofstream map(out + ".geno.map");
    ofstream indfile(out + ".geno.ind");
    
    string line;
    vector<string> l;
    vector<string> ind;
    size_t n, m;
    
    // Skip Header
    string prefix("#CHROM");
    bool have_header = false;
    while (file) {
        getline(file, line);
        if (!line.compare(0, prefix.size(), prefix)) {
            have_header = true;
            break; 
        }
    }
    if (!have_header) {
        Rcpp::stop("wrong VCF file, no line begin with \"#CHROM\".");
    }
    
    // Write inds to file.
    ind = split_line(line);
    // boost::split(ind, line, boost::is_any_of("\t"));
    vector<string>(ind.begin() + 9, ind.end()).swap(ind);   // delete first 9 columns
    int indn = ind.size();
    for (int i = 0; i < indn; i++) {
        indfile << ind[i] << endl;
    }
    n = ind.size();
    indfile.close();

    // unsigned m = count(                 // 3-4s per 100Mb
    //     istream_iterator<char>(file),
    //     istream_iterator<char>(),
    //     '\n');
    // file.seekg(pos);
    
    // Write inds to file.
    map << "SNP\tCHROM\tPOS\tREF\tALT"<< endl;
    m = 0;
    while (getline(file, line)) {
        string tmp = line.substr(0, MAP_INFO_N);
        l = split_line(tmp);
        // boost::split(l, tmp, boost::is_any_of("\t"));

        if (l[2] == ".") {      // snp name missing
            l[2] = l[0] + '-' + l[1];
        }
        map << l[2] << '\t' << l[0] << '\t' << l[1] << '\t' << l[3] << '\t' << l[4] << endl;
        m++;
    }
    map.close();
    file.close();
    
    return List::create(_["n"] = n,
                        _["m"] = m);
}


double vcf_marker_parser(string m, double NA_C) {
    if (('0' == m[0] || '1' == m[0]) && ('0' == m[2] || '1' == m[2])) {
        return (m[0] - '0' + m[2] - '0');
    } else {
        return NA_C;
    }
}

template <typename T>
void vcf_parser_genotype(std::string vcf_file, XPtr<BigMatrix> pMat, long maxLine, double NA_C, int threads=0, bool verbose=true) {
    // define
    ifstream file(vcf_file);
    
    omp_setup(threads);
    string line;
    vector<string> l;
    vector<char> markers;
    size_t m;
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    // progress bar
    MinimalProgressBar_perc pb;
    Progress progress(pMat->ncol(), verbose, pb);
    
    // Skip Header
    string prefix("#CHROM");
    bool have_header = false;
    while (file) {
        getline(file, line);
        if (!line.compare(0, prefix.size(), prefix)) {
            have_header = true;
            break;
        }
    }
    if (!have_header) {
        Rcpp::stop("wrong VCF file, no line begin with \"#CHROM\".");
    }
    
    // parser genotype
    m = 0;
    vector<string> buffer;
    while (file) {
        buffer.clear();
        for (int i = 0; (maxLine <= 0 || i < maxLine) && getline(file, line); i++) {
            if (line.length() > 1) {    // Handling the blank line at the end of the file.
                buffer.push_back(line);
            }
        }
        #pragma omp parallel for private(l, markers)
        for (std::size_t i = 0; i < buffer.size(); i++) {
            l = split_line(buffer[i]);
            // boost::split(l, buffer[i], boost::is_any_of("\t"));
            // markers.clear();
            // vector<string>(l.begin() + 9, l.end()).swap(l);
            // transform(
            //     l.begin(), l.end(), 
            //     back_inserter(markers),
            //     boost::bind<T>(&vcf_marker_parser<T>, _1, NA_C)
            // );
            // for (std::size_t j = 0; j < markers.size(); j++) {
            //     mat[j][m + i] = markers[j];
            // }
            for(size_t j = 0; j < (l.size() - 9); j++) {
                mat[m + i][j] = static_cast<T>(vcf_marker_parser(l[j + 9], NA_C));
            }
        }
        progress.increment(buffer.size());
        m += buffer.size();
    }
}

// [[Rcpp::export]]
void vcf_parser_genotype(std::string vcf_file, SEXP pBigMat, long maxLine, int threads=0, bool verbose=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return vcf_parser_genotype<char>(vcf_file, xpMat, maxLine, NA_CHAR, threads, verbose);
    case 2:
        return vcf_parser_genotype<short>(vcf_file, xpMat, maxLine, NA_SHORT, threads, verbose);
    case 4:
        return vcf_parser_genotype<int>(vcf_file, xpMat, maxLine, NA_INTEGER, threads, verbose);
    case 8:
        return vcf_parser_genotype<double>(vcf_file, xpMat, maxLine, NA_REAL, threads, verbose);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}


// ***** HAPMAP *****

// [[Rcpp::export]]
List hapmap_parser_map(std::string hmp_file, std::string out) {
    // Define
    const int MAP_INFO_N = 1000;      // max length of "SNP, POS and CHROM"
    ofstream map(out + ".geno.map");
    ofstream indfile(out + ".geno.ind");
    
    string line;
    vector<string> l;
    vector<string> ind;
    vector<string> Major;

    size_t n;
    size_t m;
    
    // for (int i = 0; i < hmp_file.size(); i++) { // TODO(haohao): Support mulit hmp file.
        ifstream file(hmp_file);
        
        // Skip Header
        string prefix("rs#");
        bool have_header = false;
        while (file) {
            getline(file, line);
            if (!line.compare(0, prefix.size(), prefix)) {
                have_header = true;
                break; 
            }
        }
        if (!have_header) {
            Rcpp::stop("wrong HAPMAP file, no line begin with \"rs#\".");
        }
        
        // Write inds to file.
        ind = split_line(line);
        // boost::split(ind, line, boost::is_any_of(" \t"));
        vector<string>(ind.begin() + 11, ind.end()).swap(ind);   // delete first 11 columns
        for (size_t i = 0; i < ind.size(); i++) {
            indfile << ind[i] << endl;
        }
        n = ind.size();
        indfile.close();
        
        // Write SNPs to map file.
        map << "SNP\tCHROM\tPOS\tREF\tALT"<< endl;
        m = 0;
        while (getline(file, line)) {
            string tmp = line.substr(0, MAP_INFO_N);
            l = split_line(tmp);
            // boost::split(l, tmp, boost::is_any_of(" \t"));
            
            if (l[0] == ".") {      // snp name missing
                l[0] = l[2] + '-' + l[3];
            }
            
            vector<string> alleles1;
            vector<string> alleles;
            // vector<string> atcg{"A", "T", "C", "G"};
            std::unordered_set<string> atcg{"A", "T", "C", "G"};

            alleles1 = split_line(l[1], "/");
            // boost::split(alleles1, l[1], boost::is_any_of("/"));
            
            for(size_t ii = 0; ii < alleles1.size(); ii++){
                if (alleles1[ii].length() != 1) {
                    Rcpp::stop(("unknown variants [" + l[1] + "] at " + to_string(m+2) + "th row of second column in HAPMAP file.").c_str());
                }
                if (atcg.find(alleles1[ii]) != atcg.end()) {
                    alleles.push_back(alleles1[ii]);
                }
            }
            if (alleles.size() == 0) {
                Rcpp::stop(("unknown variants at " + to_string(m+2) + "th row of second column in HAPMAP file.").c_str());
            }
            if (alleles.size() == 1) {
                alleles.push_back(alleles[0]);
            }
            if (alleles.size() > 2) {
                Rcpp::stop(("variants with more than 2 alleles at " + to_string(m+2) + "th row in HAPMAP file.").c_str());
            }

            Major.push_back(alleles[0]);
            map << l[0] << '\t' << l[2] << '\t' << l[3] << 
                '\t' << alleles[0] << 
                '\t' << alleles[1] << endl;
            m++;
        }
        map.close();
        file.close();
    // }
    
    return List::create(_["n"] = n,
                        _["m"] = m,
                        _["Major"] = Major);
}

double hapmap_marker_parser(string m, char major, double NA_C) {
    if (m.length() == 1) {  // Hapmap
        // Rcout << "major: " << major << '\t' << "now: " << m[0] << endl;
        if (m[0] == '+' || m[0] == '0' || m[0] == '-' || m[0] == 'N') {
            return NA_C;
        } else if (m[0] == major) {
            return 0;
        } else if (m[0] == 'R' || m[0] == 'Y' || m[0] == 'S' || m[0] == 'W' || m[0] == 'K' || m[0] == 'M') {
            return 1;
        } else if (m[0] == 'A' || m[0] == 'T' || m[0] == 'G' || m[0] == 'C') {
            return 2;
        }
    } else if (m.length() == 2) {   // Hapmap Diploid
        if ((m[0] != 'A' && m[0] != 'T' && m[0] != 'G' && m[0] != 'C') ||
            (m[1] != 'A' && m[1] != 'T' && m[1] != 'G' && m[1] != 'C')) {
            return NA_C;
        } else {
            return ((m[0] == major) ? '0' : '1') + ((m[1] == major) ? '0' : '1');
        }
    }
    return NA_C;
}

template <typename T>
void hapmap_parser_genotype(std::string hmp_file, std::vector<std::string> Major, XPtr<BigMatrix> pMat, long maxLine, double NA_C, int threads=0, bool verbose=true) {
    // define
    ifstream file(hmp_file);

    omp_setup(threads);
    string line;
    char major;
    vector<string> l;
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    // progress bar
    MinimalProgressBar_perc pb;
    Progress progress(pMat->ncol(), verbose, pb);
    
    // Skip Header
    string prefix("rs#");
    bool have_header = false;
    while (file) {
        getline(file, line);
        if (!line.compare(0, prefix.size(), prefix)) {
            have_header = true;
            break; 
        }
    }

    if (!have_header) {
        Rcpp::stop("wrong HAPMAP file, no line begin with \"rs#\".");
    }

    l = split_line(line);
    size_t n_col = l.size();
    
    // parser genotype
    size_t m = 0;
    vector<string> buffer;
    while (file) {
        buffer.clear();
        for (int i = 0; (maxLine <= 0 || i < maxLine) && getline(file, line); i++) {
            // Rcout << i << endl << line << endl;
            if (line.length() > 1) {    // Handling the blank line at the end of the file.
                buffer.push_back(line);
            }
        }
        size_t n_marker = buffer.size();

        #pragma omp parallel for private(l, major)
        for (size_t i = 0; i < n_marker; i++) {
            l = split_line(buffer[i]);
            // boost::split(l, buffer[i], boost::is_any_of(" \t"));
            if(l.size() != n_col)
                Rcpp::stop(("line " + to_string(m+i+2) + " does not have " + to_string(n_col) + " elements in HAPMAP file.").c_str());
            major = Major[m + i][0];
            for(size_t j = 0; j < (l.size() - 11); j++) {
                mat[m + i][j] = static_cast<T>(hapmap_marker_parser(l[j + 11], major, NA_C));
            }
        }
        progress.increment(n_marker);
        m += n_marker;
    }
}

// [[Rcpp::export]]
void hapmap_parser_genotype(std::string hmp_file, std::vector<std::string> Major, SEXP pBigMat, long maxLine, int threads=0, bool verbose=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return hapmap_parser_genotype<char>(hmp_file, Major, xpMat, maxLine, NA_CHAR, threads, verbose);
    case 2:
        return hapmap_parser_genotype<short>(hmp_file, Major, xpMat, maxLine, NA_SHORT, threads, verbose);
    case 4:
        return hapmap_parser_genotype<int>(hmp_file, Major, xpMat, maxLine, NA_INTEGER, threads, verbose);
    case 8:
        return hapmap_parser_genotype<double>(hmp_file, Major, xpMat, maxLine, NA_REAL, threads, verbose);
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
    l = split_line(line);
    // boost::split(l, line, boost::is_any_of("\t ,"));

    n = l.size();
    m = 1;
    while (getline(file, line))
        m++;
    
    return List::create(_["n"] = n,
                        _["m"] = m);
}


// ***** BFILE *****

template <typename T>
void write_bfile(XPtr<BigMatrix> pMat, std::string bed_file, double NA_C, bool mrkbycol = true, int threads=0, bool verbose=true) {
    // check input
    string ending = ".bed";
    if (bed_file.length() <= ending.length() ||
        0 != bed_file.compare(bed_file.length() - ending.length(), ending.length(), ending)) {
        bed_file += ending;
    }
    
    // define
    T c;
    omp_setup(threads);
    int m = (mrkbycol ? pMat->ncol() : pMat->nrow());
    int nind = (mrkbycol ? pMat->nrow() : pMat->ncol());
    int n = nind / 4;  // 4 individual = 1 bit
    if (nind % 4 != 0) 
        n++;
    
    vector<uint8_t> geno(n);
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    FILE *fout;
    fout = fopen(bed_file.c_str(), "wb");
    
    // progress bar
    MinimalProgressBar_perc pb;
    Progress progress(m, verbose, pb);
    
    // magic number of bfile
    const unsigned char magic_bytes[] = { 0x6c, 0x1b, 0x01 };
    fwrite((char*)magic_bytes, 1, 3, fout);
    
    // map
    std::map<T, int> code;
    code[0] = 3;
    code[1] = 2;
    code[2] = 0;
    code[static_cast<T>(NA_C)] = 1;
    
    // write bfile
    if(mrkbycol){
        for (int i = 0; i < m; i++) {
            #pragma omp parallel for private(c)
            for (int j = 0; j < n; j++) {
                uint8_t p = 0;
                for (int x = 0; x < 4 && (4 * j + x) < nind; x++) {
                    c = mat[i][4 * j + x];
                    p |= code[c] << (x*2);
                }
                geno[j] = p;
            }
            fwrite((char*)geno.data(), 1, geno.size(), fout);
            progress.increment();
        }
    }else{
        for (int i = 0; i < m; i++) {
            #pragma omp parallel for private(c)
            for (int j = 0; j < n; j++) {
                uint8_t p = 0;
                for (int x = 0; x < 4 && (4 * j + x) < nind; x++) {
                    c = mat[4 * j + x][i];
                    p |= code[c] << (x*2);
                }
                geno[j] = p;
            }
            fwrite((char*)geno.data(), 1, geno.size(), fout);
            progress.increment();
        }
    }
    fclose(fout);
    return;
}

// [[Rcpp::export]]
void write_bfile(SEXP pBigMat, std::string bed_file, bool mrkbycol = true, int threads=0, bool verbose=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return write_bfile<char>(xpMat, bed_file, NA_CHAR, mrkbycol, threads, verbose);
    case 2:
        return write_bfile<short>(xpMat, bed_file, NA_SHORT, mrkbycol, threads, verbose);
    case 4:
        return write_bfile<int>(xpMat, bed_file, NA_INTEGER, mrkbycol, threads, verbose);
    case 8:
        return write_bfile<double>(xpMat, bed_file, NA_REAL, mrkbycol, threads, verbose);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}

template <typename T>
void read_bfile(std::string bed_file, XPtr<BigMatrix> pMat, long maxLine, double NA_C, int threads=0, bool verbose=true) {
    // check input
    if (!boost::ends_with(bed_file, ".bed")) {
        bed_file += ".bed";
    }
    
    // define
    omp_setup(threads);
    size_t ind = pMat->nrow();
    long n = ind / 4;  // 4 individual = 1 bit
    if (ind % 4 != 0) 
        n++; 
    char *buffer;
    long buffer_size;
    MatrixAccessor<T> mat = MatrixAccessor<T>(*pMat);
    
    // map
    std::map<int, T> code;
    code[3] = static_cast<T>(0);
    code[2] = static_cast<T>(1);
    code[1] = static_cast<T>(NA_C);
    code[0] = static_cast<T>(2);
    
    // open file
    FILE *fin;
    fin = fopen(bed_file.c_str(), "rb");
    fseek(fin, 0, SEEK_END);
    long length = ftell(fin);
    rewind(fin);
    
    // get buffer_size
    buffer_size = maxLine > 0 ? (maxLine * n) : (length - 3);
    
    // progress bar
    int n_block = (length - 3) / buffer_size;
    if ((length - 3) % buffer_size != 0) { n_block++; }
    MinimalProgressBar_perc pb;
    Progress progress(n_block, verbose, pb);
    
    // magic number of bfile
    buffer = new char [3];
    size_t n_bytes_read = static_cast<size_t>(fread(buffer, 1, 3, fin));
    
    // loop file
    size_t cond;
    long block_start;
    for (int i = 0; i < n_block; i++) {
        buffer = new char [buffer_size];
        n_bytes_read = static_cast<size_t>(fread(buffer, 1, buffer_size, fin));
        
        // i: current block, j: current bit.
        block_start = i * buffer_size;
        cond = min(buffer_size, length - 3 - block_start);

        #pragma omp parallel for
        for (size_t j = 0; j < cond; j++) {
            // bit -> item in matrix
            size_t r = j / n + i * maxLine;
            size_t c = j % n * 4;
            uint8_t p = buffer[j];
            
            for (size_t x = 0; x < 4 && (c + x) < ind; x++) {
                mat[r][c + x] = code[(p >> (2*x)) & 0x03];
            }
        }
        progress.increment();
    }
    fclose(fin);
    return;
}

// [[Rcpp::export]]
void read_bfile(std::string bed_file, SEXP pBigMat, long maxLine, int threads=0, bool verbose=true) {
    XPtr<BigMatrix> xpMat(pBigMat);
    
    switch(xpMat->matrix_type()) {
    case 1:
        return read_bfile<char>(bed_file, xpMat, maxLine, NA_CHAR, threads, verbose);
    case 2:
        return read_bfile<short>(bed_file, xpMat, maxLine, NA_SHORT, threads, verbose);
    case 4:
        return read_bfile<int>(bed_file, xpMat, maxLine, NA_INTEGER, threads, verbose);
    case 8:
        return read_bfile<double>(bed_file, xpMat, maxLine, NA_REAL, threads, verbose);
    default:
        throw Rcpp::exception("unknown type detected for big.matrix object!");
    }
}
