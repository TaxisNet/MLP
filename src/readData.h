#ifndef READDATA_H_INCLUDED
#define READDATA_H_INCLUDED

#include <vector>

extern void readData( int , char** , int* , double *** );
extern void readDataFromJson( const char* , int* , double *** , std::vector<double>& );
#endif // READDATA_H_INCLUDED