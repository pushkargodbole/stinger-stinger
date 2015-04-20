#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>

template <class T>
void print(T a)
{
    cout << a;
}

template <class T>
void print(vector<T> & a)
{
    cout << "vector<";
    for(int i=0;i<a.size();++i)
    {
        print(a[i]);
        if(i<a.size()-1) cout << ", ";
    }
    cout << ">" << endl;
}
/*
template <class T, class A>
void print(map<T, A> & a)
{
    cout << "map<";
    for(map<T, A>::iterator it=a.begin(); it!=a.end(); ++it)
    {
        cout << a->first << ":";
        print(a->second);
        if(i<a.size()-1) cout << ", ";
    }
    cout << ">" << endl;
}
*/
