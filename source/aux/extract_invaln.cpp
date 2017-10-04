#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <map>
#include <set>
#include <cstdlib>
#include <array>

#include <bioinfo/util.h>
#include <bioinfo/multifastx.h>
#include <loonutil/util.h>
#include <loonutil/exception.h>

using namespace std;


/*********************************** class Interval **********************************************/
typedef array<int, 2> Interval;
istream& operator>>(istream& in, Interval& obj)
{
    return in >> obj[0] >> obj[1];
}

ostream& operator<<(ostream& out, const Interval& obj)
{
    return out << obj[0] << ' ' << obj[1];
}

/******************** class Interval variables ********************/
//vector<vector<Interval> > validated_segments;
//vector<map<Interval, size_t> > validated_segment_to_idx;
vector<set<Interval> > inversions;
Interval tmp_interval;
/*********************************** class Interval end **********************************************/


/*************** variables ************************/
map<string, string> fastx;
bool raw_file_is_known = false;
bool is_fasta = true;
loon::MultiFasta raw_reads_fasta;
loon::MultiFastq raw_reads_fastq;
loon::MultiFasta reference_genome;
set<size_t> used_ref_ids;
string common_outputprefix;
/*********************************** functions **********************************************/

void check_file_open(ifstream& fin)
{
    if(! fin.is_open())
    {
        cerr << "[ERROR]: Cannot read the file" << endl;
        exit(1);
    }
}

void check_file_open(ofstream& fout)
{
    if(! fout.is_open())
    {
        cerr << "[ERROR]: Cannot write the file" << endl;
        exit(1);
    }
}

// void read_graph_file_seg(const string& fname)
// {
//     ifstream fin(fname.c_str());
//     check_file_open( fin );
//     size_t n;
//     if(!(fin >> n))
//     {
//         cout << "No files to analyze" << endl;
//         exit(0);
//     }
//     validated_segments.resize( n );
//     validated_segment_to_idx.resize( n );
// 
//     size_t n_seg, ref_id;
//     for(size_t i = 0; i < n; ++i)
//     {
//         fin >> ref_id >> n_seg;
//         validated_segments[ ref_id ].resize( n_seg );
//         for(size_t j = 0; j < n; ++j)
//         {
//             fin >> validated_segments[ ref_id ][ j ];
//             validated_segment_to_idx[ validated_segments[ref_id][j] ] = j;
//         }
//     }
//     fin.close();
// }

size_t parse_qid(const string& qname)
{
    size_t ret = 0;
    size_t i = 0;
    while(qname[i] < '0' || qname[i] > '9')
        ++i;
    while(qname[i] >= '0' && qname[i] <= '9')
        ret = ret * 10 + qname[i++] - '0';
    return ret;
}

void read_inversion_report(const string& fname)
{
    ifstream fin( fname.c_str() );
    check_file_open( fin );
    size_t n, ref_id;
    while(fin >> n)
    {
        fin >> ref_id;
        if(ref_id + 1 >= inversions.size())
            inversions.resize( ref_id + 1 );
        for(size_t i = 0; i < n; ++i)
        {
            fin >> tmp_interval;
            inversions[ ref_id ].insert( tmp_interval );
        }
    }
    fin.close();
}

string build_filename(size_t ref_id, const string& suffix)
{
    return to_string(ref_id) + ".inv." + suffix;
}

void write_raw_fastx(const string& qname, const Interval& seg, size_t ref_id)
{
    if(is_fasta)
    {
        string fname = common_outputprefix + to_string(ref_id) + ".inv.fasta";
        ofstream fout;
        if(used_ref_ids.count( ref_id ) == 0)
        {
            used_ref_ids.insert( ref_id );
            fout.open(fname.c_str());
        }
        else
            fout.open(fname.c_str(), ios_base::app);

        size_t qid = parse_qid( qname );
        fout << ">" << raw_reads_fasta[ qid ].name << " " << raw_reads_fasta[qid].comment << " ; " << qname << ' ' << seg << endl;
        fout << raw_reads_fasta[qid].sequence << endl;
        fout.close();
    }
    else
    {
        string fname = common_outputprefix + to_string(ref_id) + ".inv.fastq";
        ofstream fout;
        if(used_ref_ids.count( ref_id ) == 0)
        {
            used_ref_ids.insert( ref_id );
            fout.open( fname.c_str() );
        }
        else
            fout.open( fname.c_str(), ios_base::app );

        size_t qid = parse_qid(qname);
        fout << "@" << raw_reads_fastq[ qid ].name << ' ' << raw_reads_fastq[qid].comment << " ; " << qname << ' ' << seg << endl;
        fout << raw_reads_fastq[qid].sequence << endl;
        fout << "+" << endl;
        fout << raw_reads_fastq[qid].quality << endl;
        fout.close();
    }
}

void generate_inversion_alignment(const string& graph_bridge_file, const string& outfile)
{
    ifstream fin( graph_bridge_file.c_str() );
    check_file_open( fin );
    ofstream fout( outfile.c_str() );
    check_file_open( fout );

    
    size_t ref_id;
    while(fin >> ref_id)
    {
        string qname[2];
        long long qlen[2];
        Interval r_aln[2], q_aln[2];
        short mapping_quality[2];
        char direction[2];
        size_t seg_id[2];
        Interval seg[2];

        char tmp_star;
        for(int ii = 0; ii < 2; ++ii)
        {
            if( ii == 1 )   fin >> tmp_star;
            fin >> qname[ii] >> qlen[ii] >> r_aln[ii] >> q_aln[ii]
                >> mapping_quality[ii] >> direction[ii]
                >> seg_id[ii] >> seg[ii];
        }
        bool tried_else = false;
        try
        {
            // if(inversions.at( ref_id ).count( r_aln[0] ) == 1 || inversions.at( ref_id ).count( r_aln[1] ) == 1)
            // {
            // }
            if(inversions.at( ref_id ).count( seg[0] ) == 1)
            {
                if(raw_file_is_known)
                    write_raw_fastx(qname[0], seg[0], ref_id);
                else
                    fout << qname[0] << ' ' << fastx[ qname[0] ] << ' ' << seg[0] << ' ' << ref_id << endl;
            }
            else
            {
                tried_else = true;
                if(inversions.at( ref_id ).count( seg[1] ) == 1)
                {
                    if(raw_file_is_known)
                        write_raw_fastx(qname[1], seg[1], ref_id);
                    else
                        fout << qname[1] << ' ' << fastx[ qname[1] ] << ' ' << seg[1] << ' ' << ref_id << endl;
                }
            }
        }catch(const out_of_range& oor)
        {
            try
            {
                if(!tried_else && inversions.at( ref_id ).count( seg[1] ) == 1)
                {
                    if(raw_file_is_known)
                        write_raw_fastx(qname[1], seg[1], ref_id);
                    else
                        fout << qname[1] << ' ' << fastx[ qname[1] ] << ' ' << seg[1] << ' ' << ref_id << endl;
                }
            }catch(const out_of_range& oor)
            {
            }
        }
    }

    if(raw_file_is_known)
    {
        ofstream fout2( common_outputprefix + "ref_maps.txt" );
        for(set<size_t>::const_iterator it = used_ref_ids.begin(); it != used_ref_ids.end(); ++it)
            fout2 << (*it) << ' ' << reference_genome[ *it ].name << endl;
        fout2.close();
    }

    fin.close();
    fout.close();
}

string parse_header(const string& header, char delimiter=' ')
{
    size_t space_pos = header.find(delimiter);
    if(space_pos != string::npos)
        return header.substr(1, space_pos - 1);
    else
        return header.substr(1);
}

void read_fastx( const string& fname )
{
    ifstream fin( fname.c_str() );
    string header, seq = "", line;
    bool dont_read_next = false;
    while(dont_read_next || getline(fin, line))
    {
        dont_read_next = false;
        if(line[0] == '>')
        {
            header = parse_header( line, '/' );
            while( getline(fin, line) )
            {
                if(line[0] == '>' || line[0] == '@')
                {
                    dont_read_next = true;
                    break;
                }
                seq += line; 
            }
        }
        else if(line[0] == '@')
        {
            header = parse_header( line, '/' );
            getline(fin, seq);
            getline(fin, line);
            getline(fin, line);
        }

        if(!seq.empty())
        {
            fastx[ header ] = seq;
            seq = "";
        }
    }
    if(!seq.empty())
        fastx[ header ] = seq;
    fin.close();
}

void read_raw_reads(const string& raw_reads_file)
{
    std::string ftype = loon::fileformat_by_suffix( raw_reads_file );
    if(ftype == "fasta")
    {
        is_fasta = true;
        raw_reads_fasta.read( raw_reads_file );
    }
    else if(ftype == "fastq")
    {
        is_fasta = false;
        raw_reads_fastq.read( raw_reads_file );
    }
    else
        throw loon::Exception(5, __LINE__, __FILE__, "File type error: [%s]", ftype.c_str());
}


int main(int argc, char* argv[])
{
    if(argc == 8)
    {
        raw_file_is_known = true;
        cerr << "[INFO]: read raw reads" << endl;
        read_raw_reads( argv[5] ); // load raw reads file
        cerr << "[INFO]: read raw reference" << endl;
        reference_genome.read( argv[6] ); // load raw reference file
        common_outputprefix = string(argv[7]) + "/";
    }
    cerr << "[INFO]: read fastx file" << endl;
    read_fastx( argv[1] );
    cerr << "[INFO]: read inversion report" << endl;
    read_inversion_report( argv[2] );
    cerr << "[INFO]: generate inversion alignment" << endl;
    generate_inversion_alignment( argv[3], argv[4]);
    cerr << "[INFO]: done!" << endl;
    return 0;
}

