#include<iostream>
#include<fstream>
#include<string>
using namespace std;
int main()
{
	string name1("AAAsources");
	string name2("../testfile/sources");
//	string name1("AAAsoureres.txt");
//	string name2("../FM-CIndex/AAAsoureres.txt");
	ifstream fin1,fin2;
	fin1.open(name1.c_str(),ios_base::in);
	fin2.open(name2.c_str(),ios_base::in);
	string tmp1(""), tmp2("");
	int i = 0;
	while(getline(fin1,tmp1) && getline(fin2,tmp2))
	{
		i++;
		if(tmp1 != tmp2)
		{
			cout<<"NO!!!!!!!!!!!!"<<endl;
			cout<<tmp1<<"	"<<tmp2<<"	"<<i<<endl;
//			return 0;
		}
		tmp1 = "";
		tmp2 = "";
		if(i % 10000 == 0)cout<<i<<endl;
	}
	cout<<"YES!!!!!!!!!!!!"<<endl;
	return 0;
}