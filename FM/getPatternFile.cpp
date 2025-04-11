#include<stdlib.h>
#include<fstream>
#include<vector>
#include<string>
#include<unistd.h>
#include<cstdlib>
#include<getopt.h>
#define i64 long long

using namespace std;
static const char *optString = "i:o:l:n:";
static const struct option longOpts[]={
    {"input",required_argument,NULL,'i'},
    {"output",required_argument,NULL,'o'},
    {"patternlength",required_argument,NULL,'l'},
    {"mount",required_argument,NULL,'n'},
};
void get_patterns(const char *file,vector<string> &patterns, vector<i64> &pos,int seed,int runtime,int patternlen)
{
    srand(unsigned(seed));
    FILE *fr = fopen(file, "r+");
    /* assert(fr); */
    fseek(fr, 0L, SEEK_END);
    i64 size = ftell(fr);
    char *text = new char[size + 1]; text[size] = '\0';
    rewind(fr);
    auto readnum = fread(text, sizeof(char), size, fr);
    string stext = text;
    pos.clear();
    unsigned long long  position;
    for(int i=0;i <runtime;i++)
    {
        position = rand() % (size-patternlen);
        /* cout<<position<<endl; */
        if(position + patternlen< size)
        {
            string tmp = stext.substr(position, patternlen);
            if(tmp.find("\n")!=tmp.npos||tmp.find("N")!=tmp.npos)
            {
                i--;
                continue;
            }
            // cout<<tmp<<endl;
            patterns.push_back(tmp);
            pos.push_back(position);
        }
    }
    if(text) delete [] text;
    if(fr) fclose(fr);
}

int main(int argc,char** argv)
{
    string inputfile = "";
    string outputfile = ".pattern";
    int patternLen = 4;
    int patternMount = 10;
    int ch;
    int longIndex =0;
    while((ch=getopt_long(argc,argv,optString,longOpts,&longIndex))!=-1){
        switch(ch){
            case'i':
                inputfile = string(optarg);
                break;
            case'o':
                outputfile = string(optarg);
                break;
            case'l':
                patternLen= atoi(optarg);
                break;
            case'n':
                patternMount= atoi(optarg);
                break;
        }
    }
    vector<string> pattern;
    vector<i64> pos;
    get_patterns(inputfile.c_str(),pattern,pos,123,patternMount,patternLen);
    fstream s(outputfile.c_str(),s.out);
    for(int i = 0 ;i<pattern.size();i++){
        s<<pattern[i]<<endl;
    }
return 0;
}
