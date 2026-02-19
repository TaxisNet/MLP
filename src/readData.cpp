#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

#
#include <nlohmann/json.hpp>

using namespace std;

using json = nlohmann::json;

double CalcDistEuc ( double *X, double *Y, int I, int J );
double CalcDistAtt ( double *X, double *Y, int I, int J );
void CalcLatLong ( double *X, double *Y, int n, double *latit, double* longit );
double CalcDistGeo ( double *latit, double *longit, int I, int J );

void readDataFromJson(const char* jsonPath, int* Dimension, double*** Mdist, vector<double>& Mweights) {
    // Read JSON file
    std::ifstream file(jsonPath);
    if (!file.is_open()) {
        std::cout << "Cannot open JSON file: " << jsonPath << std::endl;
        exit(1);
    }
    

    // cout << "Reading data from JSON file: " << jsonPath << std::endl;
    json j;
    file >> j;
    
    // Get the distance matrix from JSON
    auto distMatrix = j["distance_matrix"];
    int N = distMatrix.size();
    
    // Check if node weights are provided
    if (j.contains("weights")) {
        auto nodeWeights = j["weights"];
        int M = nodeWeights.size();
        Mweights.resize(M+1); // 1-indexed
        for (int i = 0; i < M; i++) {
           Mweights[i+1] = nodeWeights[i].get<double>();
        }
    }



    // Allocate 2D array
    double** dist = new double*[N+1];
    for (int i = 0; i < N+1; i++) {
        dist[i] = new double[N+1];
    }
    
    // Copy data from JSON (1-indexed)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            dist[i+1][j+1] = distMatrix[i][j].get<double>();
        }
    }
    
    *Dimension = N;
    *Mdist = dist;
}


void readData( int argc, char** argv, int* Dimension, double ***Mdist )
{

     if (argc < 2) {
         cout << "\nFaltando parametros\n";
         cout << " ./exec [Instancia] "<< endl;
         exit(1);
     }

     if (argc > 2) {
          cout << "\nMuitos parametros\n";
          cout << " ./exec [Instancia] " << endl;
         exit(1);
     }

    int N;
    string arquivo, ewt;

     char *instancia;
     instancia = argv[1];

     //ifstream in(argv[1], ios::in);


    ifstream in( instancia, ios::in);

	if (!in) {
		cout << "arquivo nao pode ser aberto\n";
		exit(1);
    }

    while ( arquivo.compare("DIMENSION:") != 0 && arquivo.compare("DIMENSION" ) != 0 ) {
        in >> arquivo;
    }

    if ( arquivo.compare("DIMENSION" ) == 0 )  in >> arquivo;

    in >> N;

    while ( arquivo.compare("EDGE_WEIGHT_TYPE:") != 0 && arquivo.compare("EDGE_WEIGHT_TYPE" ) != 0 ) {
        in >> arquivo;
    }
    if ( arquivo.compare("EDGE_WEIGHT_TYPE" ) == 0 )  in >> arquivo;

    in >> ewt;

    double *x = new double [N+1];
    double *y = new double [N+1];

    // Alocar matriz 2D
    double **dist = new double*[N+1];

    for ( int i = 0; i < N+1; i++ ) {
        dist [i] = new double [N+1];
    }

    if ( ewt == "EXPLICIT" ) {

        while ( arquivo.compare("EDGE_WEIGHT_FORMAT:") != 0 && arquivo.compare("EDGE_WEIGHT_FORMAT" ) != 0 ) {
            in >> arquivo;
        }

        string ewf;
        if ( arquivo.compare("EDGE_WEIGHT_FORMAT" ) == 0 )  in >> arquivo;
        in >> ewf;

        if ( ewf == "FUNCTION" ) {
            cout << "FUNCTION - Not supported!" << endl; }

        else if ( ewf == "FULL_MATRIX" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            for ( int i = 1; i < N+1; i++ ) {
                for ( int j = 1; j < N+1; j++ ) {
                    in >> dist[i][j];
                }
            }
        }

        else if ( ewf == "UPPER_ROW" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            // Preencher Matriz Distancia
            for ( int i = 1; i < N; i++ ) {
                for ( int j = i+1; j < N+1; j++ ) {
                    in >> dist[i][j];
                    dist[j][i] = dist[i][j];
                }
            }

            for ( int i = 1; i < N+1; i++ ) {
                dist[i][i] = 0;
            }

        }

        else if ( ewf == "LOWER_ROW" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            // Preencher Matriz Distancia
            for ( int i = 2; i < N+1; i++ ) {
                for ( int j = 1; j < i; j++ ) {
                    in >> dist[i][j];
                    dist[j][i] = dist[i][j];
                }
            }

            for ( int i = 1; i < N+1; i++ ) {
                dist[i][i] = 0;
            }
        }

        else if ( ewf == "UPPER_DIAG_ROW" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            // Preencher Matriz Distancia
            for ( int i = 1; i < N+1; i++ ) {
                for ( int j = i; j < N+1; j++ ) {
                    in >> dist[i][j];
                    dist[j][i] = dist[i][j];
                }
            }
        }

        else if ( ewf == "LOWER_DIAG_ROW" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            // Preencher Matriz Distancia
            for ( int i = 1; i < N+1; i++ ) {
                for ( int j = 1; j <= i; j++ ) {
                    in >> dist[i][j];
                    dist[j][i] = dist[i][j];
                }
            }
        }

        else if ( ewf == "UPPER_COL" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            // Preencher Matriz Distancia
            for ( int j = 2; j < N+1; j++ ) {
                for ( int i = 1; i < j; i++ ) {
                    in >> dist[i][j];
                    dist[j][i] = dist[i][j];
                }
            }

            for ( int i = 1; i < N+1; i++ ) {
                dist[i][i] = 0;
            }

        }

        else if ( ewf == "LOWER_COL" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            // Preencher Matriz Distancia
            for ( int j = 1; j < N; j++ ) {
                for ( int i = j+1; i < N+1; i++ ) {
                    in >> dist[i][j];
                    dist[j][i] = dist[i][j];
                }
            }

            for ( int i = 1; i < N+1; i++ ) {
                dist[i][i] = 0;
            }

        }

        else if ( ewf == "UPPER_DIAG_COL" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            // Preencher Matriz Distancia
            for ( int j = 1; j < N+1; j++ ) {
                for ( int i = 1; i <= j; i++ ) {
                    in >> dist[i][j];
                    dist[j][i] = dist[i][j];
                }
            }
        }

        else if ( ewf == "LOWER_DIAG_COL" ) {

            while ( arquivo.compare("EDGE_WEIGHT_SECTION") != 0 ) {
                in >> arquivo;
            }

            // Preencher Matriz Distancia
            for ( int j = 1; j < N+1; j++ ) {
                for ( int i = j; i < N+1; i++ ) {
                    in >> dist[i][j];
                    dist[j][i] = dist[i][j];
                }
            }
        }

    }

    else if ( ewt == "EUC_2D" ) {

        while ( arquivo.compare("NODE_COORD_SECTION") != 0 ) {
            in >> arquivo;
        }
        // read coordinates
        int tempCity;
        for ( int i = 1; i < N+1; i++ ) {
            in >> tempCity >> x[i] >> y[i];
        }

        // Calculate Distance Matrix (Euclidian)
        for ( int i = 1; i < N+1; i++ ) {
            for ( int j = 1; j < N+1; j++ ) {
                dist[i][j] = floor ( CalcDistEuc ( x, y, i, j ) + 0.5 );
            }
        }
    }

    else if ( ewt == "EUD_3D" ) {
        cout << "EUC_3D - Not supported!" << endl; }

    else if ( ewt == "MAX_2D" ) {
        cout << "MAX_2D - Not supported!" << endl; }

    else if ( ewt == "MAX_3D" ) {
        cout << "MAX_3D - Not supported!" << endl; }

    else if ( ewt == "MAN_2D" ) {
        cout << "MAN_2D - Not supported!" << endl; }

    else if ( ewt == "MAN_3D" ) {
        cout << "MAN_3D - Not supported!" << endl; }

    else if ( ewt == "CEIL_2D" ) {

        while ( arquivo.compare("NODE_COORD_SECTION") != 0 ) {
            in >> arquivo;
        }
        // read coordinates
        int tempCity;
        for ( int i = 1; i < N+1; i++ ) {
            in >> tempCity >> x[i] >> y[i];
        }

        // Calculate Distance Matrix (Euclidiana)
        for ( int i = 1; i < N+1; i++ ) {
            for ( int j = 1; j < N+1; j++ ) {
                dist[i][j] = ceil ( CalcDistEuc ( x, y, i, j ) );
            }
        }
    }

    else if ( ewt == "GEO" ) {

        while ( arquivo.compare("NODE_COORD_SECTION") != 0 ) {
            in >> arquivo;
        }
        // read coordinates
        int tempCity;
        for ( int i = 1; i < N+1; i++ ) {
            in >> tempCity >> x[i] >> y[i];
        }

        double *latitude = new double [N+1];
        double *longitude = new double [N+1];

        CalcLatLong ( x, y, N, latitude, longitude );

        // Calculate Distance Matrix
        for ( int i = 1; i < N+1; i++ ) {
            for ( int j = 1; j < N+1; j++ ) {
                dist[i][j] = CalcDistGeo ( latitude, longitude, i, j );
            }
        }

    }

    else if ( ewt == "ATT" ) {

        while ( arquivo.compare("NODE_COORD_SECTION") != 0 ) {
            in >> arquivo;
        }

        // read coordinates
        int tempCity;
        int *tempX = new int [N+1];
        int *tempY = new int [N+1];

        for ( int i = 1; i < N+1; i++ ) {
            in >> tempCity >> tempX[i] >> tempY[i];
            x[i]=tempX[i];
            y[i]=tempY[i];
        }

        // Calculate Distance Matrix (Pesudo-Euclidian)
        for ( int i = 1; i < N+1; i++ ) {
            for ( int j = 1; j < N+1; j++ ) {
                dist[i][j] = CalcDistAtt ( x, y, i, j );
            }
        }

    }

    else if ( ewt == "XRAY1" ) {
        cout << "XRAY1 - Not supported!" << endl; }

    else if ( ewt == "XRAY2" ) {
        cout << "XRAY2 - Not supported!" << endl; }

    else if ( ewt == "SPECIAL" ) {
        cout << "SPECIAL - Not supported!" << endl; }

    *Dimension = N;
    *Mdist = dist;
}

double CalcDistEuc ( double *X, double *Y, int I, int J )
{
    return
    sqrt ( pow ( X[I] - X[J], 2 ) + pow ( Y[I] - Y[J], 2 ) );
}

double CalcDistAtt ( double *X, double *Y, int I, int J )
{
    // Calculate Pseudo Euclidian Distance 
    double rij, tij, dij;

    rij = sqrt ( ( pow ( X[I] - X[J], 2 ) + pow ( Y[I] - Y[J], 2 ) ) / 10 );
    tij = floor ( rij + 0.5 );

    if ( tij < rij )
        dij = tij + 1;
    else
        dij = tij;

    return dij;
}

void CalcLatLong ( double *X, double *Y, int n, double *latit, double* longit )
{
    double PI = 3.141592, min;
    int deg;

    for ( int i = 1; i < n+1; i++ ) {
        deg = (int) X[i];
        min = X[i] - deg;
        latit[i] = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
    }

    for ( int i = 1; i < n+1; i++ ) {
        deg = (int) Y[i];
        min = Y[i] - deg;
        longit[i] = PI * (deg + 5.0 * min / 3.0 ) / 180.0;
    }
}

double CalcDistGeo ( double *latit, double *longit, int I, int J )
{
    double q1, q2, q3, RRR = 6378.388;

    q1 = cos( longit[I] - longit[J] );
    q2 = cos( latit[I] - latit[J] );
    q3 = cos( latit[I] + latit[J] );

    return
    (int) ( RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3) ) + 1.0);
}
