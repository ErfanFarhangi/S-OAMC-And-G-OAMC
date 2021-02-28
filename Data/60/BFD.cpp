#include <iostream>
#include<fstream>
using namespace std;

#include <iostream> 
#include <cmath>
#include <math.h> 
#include<fstream>
#include<sstream>
#include <string> 
#include <stdlib.h>     /* srand, rand */
#include <vector>
#include <random>
#include <chrono>
#include <tuple>
#include <queue>
#include <time.h>
#include <stack> 
#include <list>
//#include<ilcplex/ilocplex.h>

//ILOSTLBEGIN

# define M_PI           3.14159265358979323846  /* pi */

using namespace std;
using namespace std::chrono;

struct Struct
{
	int p;
	int b;
	int value;
	bool temp[60];
};

struct ArrayDP
{
	double p;
	double b;
	bool temp[60];
};

struct Users
{
	double x;
	double y;
	double id;
	double md;
	double w;
	double p;
	double b;
	double total; /* summation of p and b */
	double max_p;
	double max_b;

	int index; /* to specify the number of users*/
};

struct Cloudlets
{
	double x;
	double y;
	double p;
	double b;
	int index;
};

struct values
{
	double delay;
	int index;
};
double ComputationTime = 0;
double MigrationTime = 0;
double OffloadingTime = 0;



double epsilon = 0.9; /* epsilon as an input ?*/

int user_total = 60; /* determines number of users */
int cloudlet_total = 12; /* determines number of cloudlets */

bool temp1[60];  /* user_total */
bool temp2[60];  /* user_total */

bool temp3[60]; /* inside Assignment function*/
double** UsageP = new double*[cloudlet_total]; /* user_total is the number of users */

double** UsageB = new double*[cloudlet_total]; /* user_total is the number of users */

ArrayDP* Array = new ArrayDP[user_total*10*user_total];


double total_size = 0;
int w = 5; /* window size*/
//double gamma = 0.3; /* gamma -- weight */
double theta = 0.005*0.67;  /*  propagation time (IN MICROSECOND) per each unit of distance (meter) */
int DPRAC_Limit = 10; /* for DPRAC */
int load_balance = 0; /* for load balancing */
Users* temporary = new Users[w];  /* window size -- used in function PREDICT */
Users** temporary_test = new Users*[user_total];  /* used in SpecPredict */

Users** first_prediction = new Users*[user_total];
Users** second_prediction = new Users*[user_total];
Users** third_prediction = new Users*[user_total];
Users** fourth_prediction = new Users*[user_total];
Users** fifth_prediction = new Users*[user_total];


int** V = new int*[user_total]; /* V is used in OAMC and changes in every t */
double** r = new double*[user_total]; /* r is used in OAMC and changes in every t -- Rounding */
int** Temp_V = new int*[user_total]; /* V is used in OAMC and changes in every t */
int** F = new int*[user_total]; /* F is used in Assignment */
int** N = new int*[user_total]; /* N is used in Assignment */

int* Migrations = new int[user_total]; /* records the number of migrations of each application */
double* Latency = new double[user_total]; /* records the latency of each application in its lifetime */

int migration = 0;
int migration1 = 0;

int MigrationNewYork = 0;
int MigrationOrlando = 0;
int MigrationDACT = 0;

double min_x = INT_MAX;
double min_y = INT_MAX;
double max_x = INT_MIN;
double max_y = INT_MIN;
int T = INT_MAX; /* is the intersection of number of time slots */
Users** users = new Users*[user_total]; /* user_total is the number of users */
Cloudlets* cloudlets = new Cloudlets[cloudlet_total];
Cloudlets* C = new Cloudlets[cloudlet_total]; /* used in OAMC */
int* count_time = new int[user_total]; /* shows how many time slots each user runs */
int** R = new int*[user_total]; /* holds the index of the assigned cloudlet */

bool** S = new bool*[cloudlet_total]; /* holds the index of the assigned cloudlet */
bool* I = new bool[user_total]; /* used in DPRAC */

Users** First_Predict = new Users*[user_total]; /* user_total is the number of users */
Users** Second_Predict = new Users*[user_total]; /* user_total is the number of users */
Users** Third_Predict = new Users*[user_total]; /* user_total is the number of users */
Users** Fourth_Predict = new Users*[user_total]; /* user_total is the number of users */
Users** Fifth_Predict = new Users*[user_total]; /* user_total is the number of users */

bool flag_test = false;

//typedef IloArray<IloBoolVarArray> BoolVarMatrix;
//typedef IloArray<BoolVarMatrix>   BoolVar3Matrix;
//typedef IloArray<BoolVar3Matrix>   BoolVar4Matrix;

double distance(double lat1, double lon1, double lat2, double lon2)
{
	return sqrt(pow(lat2 - lat1, 2) +
		pow(lon2 - lon1, 2) * 1.0);
}
static double haversine(double lat1, double lon1, double lat2, double lon2)
{
	/* distance between latitudes and longitudes */
	double dLat = (lat2 - lat1) * M_PI / 180.0;
	double dLon = (lon2 - lon1) *  M_PI / 180.0;

	/* convert to radians */
	lat1 = (lat1)* M_PI / 180.0;
	lat2 = (lat2)* M_PI / 180.0;

	/* apply formulae */
	double a = pow(sin(dLat / 2), 2) + pow(sin(dLon / 2), 2) * cos(lat1) * cos(lat2);
	double rad = 6371;
	double c = 2 * asin(sqrt(a));
	return rad * c;
}

void Initialize()
{
	for (int i = 0;i < user_total;i++)
	{
		Migrations[i] = 0;
		Latency[i] = 0;
	}
}

void Test()
{

	int d = 0;
	for (int i = 0;i<user_total;i++)
		d += Migrations[i];

	double d1 = 0;
	for (int i = 0;i<user_total;i++)
		d1 += Latency[i];
	cout << "migration  " << d << "  latency  " << d1 << endl;
}
void Print()
{
	for (int i = 0;i < user_total;i++)
	{
		cout<<"migration "<<Migrations[i]<<"  Latency "<<Latency[i]<<endl;
	}
}
bool compareTwoApplicationsBFD(Users user_1, Users user_2)
{
	return user_1.total >= user_2.total;
}
bool compareTwoCloudletsBFD(Cloudlets cloudlet_1, Cloudlets cloudlet_2)
{
	return cloudlet_1.p >= cloudlet_2.p && cloudlet_1.b >= cloudlet_2.b; 
}
bool compareTwoDelays(values value_1, values value_2)
{
	return value_1.delay <= value_2.delay;
}
void Seperator(string str, int User_Num, int timeslot)
{

	stringstream s(str); /* Used for breaking words */
	string word; /* to store individual words */

	int count = 0;
	while (s >> word)
	{

		stringstream geek(word);

		double read = 0;
		geek >> read;

		cout.precision(word.size());

		if (count == 0) /* latitude */
		{
			users[User_Num][timeslot - 1].x = read;
			
				if (read <= min_x)
					min_x = read;
				if (read>max_x)
					max_x = read;
					
		}
		if (count == 1) /* longitude */
		{
			users[User_Num][timeslot - 1].y = read;
			
			if (read <= min_y)
				min_y = read;
			if (read>max_y)
				max_y = read;

		}

		count++;
	}

}

void Results_MatrixCompletion()
{
	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{

		First_Predict[i] = new Users[200];
		Second_Predict[i] = new Users[200];
		Third_Predict[i] = new Users[200];
		Fourth_Predict[i] = new Users[200];
		Fifth_Predict[i] =  new Users[200];

	}

	cout.precision(15);

  int t=1;
    
  while(t<=19)
  { 
	stringstream out;
	out << t;
	ifstream data("results_Sample_Data_60_b_"+out.str()+".csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;


	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double b1 = 0;
		geek >> b1;

		stringstream geek2(two);
		double b2 = 0;
		geek2 >> b2;

		stringstream geek3(three);
		double b3 = 0;
		geek3 >> b3;

		stringstream geek4(four);
		double b4 = 0;
		geek4 >> b4;

		stringstream geek5(five);
		double b5 = 0;
		geek5 >> b5;

		stringstream geek6(six);
		double b6 = 0;
		geek6 >> b6;

		stringstream geek7(seven);
		double b7 = 0;
		geek7 >> b7;

		stringstream geek8(eight);
		double b8 = 0;
		geek8 >> b8;

		stringstream geek9(nine);
		double b9 = 0;
		geek9 >> b9;

		stringstream geek10(ten);
		double b10 = 0;
		geek10 >> b10;

		stringstream geek11(eleven);
		double b11 = 0;
		geek11 >> b11;

		stringstream geek12(twelve);
		double b12 = 0;
		geek12 >> b12;

		stringstream geek13(thirteen);
		double b13 = 0;
		geek13 >> b13;

		stringstream geek14(fourteen);
		double b14 = 0;
		geek14 >> b14;

		stringstream geek15(fifteen);
		double b15 = 0;
		geek15 >> b15;

		stringstream geek16(sixteen);
		double b16 = 0;
		geek16 >> b16;

		stringstream geek17(seventeen);
		double b17 = 0;
		geek17 >> b17;

		stringstream geek18(eighteen);
		double b18 = 0;
		geek18 >> b18;

		stringstream geek19(nineteen);
		double b19 = 0;
		geek19 >> b19;

		stringstream geek20(twenty);
		double b20 = 0;
		geek20 >> b20;

		if (t == 1)
		{
		  First_Predict[user_count][t].b = b2;
		  Second_Predict[user_count][t].b = b3;
		  Third_Predict[user_count][t].b = b4;
		  Fourth_Predict[user_count][t].b = b5;
		  Fifth_Predict[user_count][t].b = b6;
		}
		if (t == 2)
		{
			First_Predict[user_count][t].b = b3;
			Second_Predict[user_count][t].b = b4;
			Third_Predict[user_count][t].b = b5;
			Fourth_Predict[user_count][t].b = b6;
			Fifth_Predict[user_count][t].b = b7;
		}
		if (t == 3)
		{
			First_Predict[user_count][t].b = b4;
			Second_Predict[user_count][t].b = b5;
			Third_Predict[user_count][t].b = b6;
			Fourth_Predict[user_count][t].b = b7;
			Fifth_Predict[user_count][t].b = b8;
		}
		if (t == 4)
		{
			First_Predict[user_count][t].b = b5;
			Second_Predict[user_count][t].b = b6;
			Third_Predict[user_count][t].b = b7;
			Fourth_Predict[user_count][t].b = b8;
			Fifth_Predict[user_count][t].b = b9;
		}
		if (t == 5)
		{
			First_Predict[user_count][t].b = b6;
			Second_Predict[user_count][t].b = b7;
			Third_Predict[user_count][t].b = b8;
			Fourth_Predict[user_count][t].b = b9;
			Fifth_Predict[user_count][t].b = b10;
		}
		if (t == 6)
		{
			First_Predict[user_count][t].b = b7;
			Second_Predict[user_count][t].b = b8;
			Third_Predict[user_count][t].b = b9;
			Fourth_Predict[user_count][t].b = b10;
			Fifth_Predict[user_count][t].b = b11;
		}
		if (t == 7)
		{
			First_Predict[user_count][t].b = b8;
			Second_Predict[user_count][t].b = b9;
			Third_Predict[user_count][t].b = b10;
			Fourth_Predict[user_count][t].b = b11;
			Fifth_Predict[user_count][t].b = b12;
		}
		if (t == 8)
		{
			First_Predict[user_count][t].b = b9;
			Second_Predict[user_count][t].b = b10;
			Third_Predict[user_count][t].b = b11;
			Fourth_Predict[user_count][t].b = b12;
			Fifth_Predict[user_count][t].b = b13;
		}
		if (t == 9)
		{
			First_Predict[user_count][t].b = b10;
			Second_Predict[user_count][t].b = b11;
			Third_Predict[user_count][t].b = b12;
			Fourth_Predict[user_count][t].b = b13;
			Fifth_Predict[user_count][t].b = b14;
		}
		if (t == 10)
		{
			First_Predict[user_count][t].b = b11;
			Second_Predict[user_count][t].b = b12;
			Third_Predict[user_count][t].b = b13;
			Fourth_Predict[user_count][t].b = b14;
			Fifth_Predict[user_count][t].b = b15;
		}
		if (t == 11)
		{
			First_Predict[user_count][t].b = b12;
			Second_Predict[user_count][t].b = b13;
			Third_Predict[user_count][t].b = b14;
			Fourth_Predict[user_count][t].b = b15;
			Fifth_Predict[user_count][t].b = b16;
		}
		if (t == 12)
		{
			First_Predict[user_count][t].b = b13;
			Second_Predict[user_count][t].b = b14;
			Third_Predict[user_count][t].b = b15;
			Fourth_Predict[user_count][t].b = b16;
			Fifth_Predict[user_count][t].b = b17;
		}
		if (t == 13)
		{
			First_Predict[user_count][t].b = b14;
			Second_Predict[user_count][t].b = b15;
			Third_Predict[user_count][t].b = b16;
			Fourth_Predict[user_count][t].b = b17;
			Fifth_Predict[user_count][t].b = b18;
		}
		if (t == 14)
		{
			First_Predict[user_count][t].b = b15;
			Second_Predict[user_count][t].b = b16;
			Third_Predict[user_count][t].b = b17;
			Fourth_Predict[user_count][t].b = b18;
			Fifth_Predict[user_count][t].b = b19;
		}
		if (t == 15)
		{
			First_Predict[user_count][t].b = b16;
			Second_Predict[user_count][t].b = b17;
			Third_Predict[user_count][t].b = b18;
			Fourth_Predict[user_count][t].b = b19;
			Fifth_Predict[user_count][t].b = b20;
		}
		if (t == 16)
		{
			First_Predict[user_count][t].b = b17;
			Second_Predict[user_count][t].b = b18;
			Third_Predict[user_count][t].b = b19;
			Fourth_Predict[user_count][t].b = b20;
		}
		if (t == 17)
		{
			First_Predict[user_count][t].b = b18;
			Second_Predict[user_count][t].b = b19;
			Third_Predict[user_count][t].b = b20;

		}
		if (t == 18)
		{
			First_Predict[user_count][t].b = b19;
			Second_Predict[user_count][t].b = b20;
		}
		if (t == 19)
		{
			First_Predict[user_count][t].b = b20;
		}
		if (t == 20)
		{
			
		}

		user_count++;
	}
	t++;
   }

   t=1;
   /* id */
   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("results_Sample_Data_60_id_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;

	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double id1 = 0;
		   geek >> id1;

		   stringstream geek2(two);
		   double id2 = 0;
		   geek2 >> id2;

		   stringstream geek3(three);
		   double id3 = 0;
		   geek3 >> id3;

		   stringstream geek4(four);
		   double id4 = 0;
		   geek4 >> id4;

		   stringstream geek5(five);
		   double id5 = 0;
		   geek5 >> id5;

		   stringstream geek6(six);
		   double id6 = 0;
		   geek6 >> id6;

		   stringstream geek7(seven);
		   double id7 = 0;
		   geek7 >> id7;

		   stringstream geek8(eight);
		   double id8 = 0;
		   geek8 >> id8;

		   stringstream geek9(nine);
		   double id9 = 0;
		   geek9 >> id9;

		   stringstream geek10(ten);
		   double id10 = 0;
		   geek10 >> id10;

		   stringstream geek11(eleven);
		   double id11 = 0;
		   geek11 >> id11;

		   stringstream geek12(twelve);
		   double id12 = 0;
		   geek12 >> id12;

		   stringstream geek13(thirteen);
		   double id13 = 0;
		   geek13 >> id13;

		   stringstream geek14(fourteen);
		   double id14 = 0;
		   geek14 >> id14;

		   stringstream geek15(fifteen);
		   double id15 = 0;
		   geek15 >> id15;

		   stringstream geek16(sixteen);
		   double id16 = 0;
		   geek16 >> id16;

		   stringstream geek17(seventeen);
		   double id17 = 0;
		   geek17 >> id17;

		   stringstream geek18(eighteen);
		   double id18 = 0;
		   geek18 >> id18;

		   stringstream geek19(nineteen);
		   double id19 = 0;
		   geek19 >> id19;

		   stringstream geek20(twenty);
		   double id20 = 0;
		   geek20 >> id20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].id = id2*100;
			   Second_Predict[user_count][t].id = id3*100;
			   Third_Predict[user_count][t].id = id4*100;
			   Fourth_Predict[user_count][t].id = id5*100;
			   Fifth_Predict[user_count][t].id = id6*100;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].id = id3*100;
			   Second_Predict[user_count][t].id = id4*100;
			   Third_Predict[user_count][t].id = id5*100;
			   Fourth_Predict[user_count][t].id = id6*100;
			   Fifth_Predict[user_count][t].id = id7*100;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].id = id4*100;
			   Second_Predict[user_count][t].id = id5*100;
			   Third_Predict[user_count][t].id = id6*100;
			   Fourth_Predict[user_count][t].id = id7*100;
			   Fifth_Predict[user_count][t].id = id8*100;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].id = id5*100;
			   Second_Predict[user_count][t].id = id6*100;
			   Third_Predict[user_count][t].id = id7*100;
			   Fourth_Predict[user_count][t].id = id8*100;
			   Fifth_Predict[user_count][t].id = id9*100;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].id = id6*100;
			   Second_Predict[user_count][t].id = id7*100;
			   Third_Predict[user_count][t].id = id8*100;
			   Fourth_Predict[user_count][t].id = id9*100;
			   Fifth_Predict[user_count][t].id = id10*100;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].id = id7*100;
			   Second_Predict[user_count][t].id = id8*100;
			   Third_Predict[user_count][t].id = id9*100;
			   Fourth_Predict[user_count][t].id = id10*100;
			   Fifth_Predict[user_count][t].id = id11*100;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].id = id8*100;
			   Second_Predict[user_count][t].id = id9*100;
			   Third_Predict[user_count][t].id = id10*100;
			   Fourth_Predict[user_count][t].id = id11*100;
			   Fifth_Predict[user_count][t].id = id12*100;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].id = id9*100;
			   Second_Predict[user_count][t].id = id10*100;
			   Third_Predict[user_count][t].id = id11*100;
			   Fourth_Predict[user_count][t].id = id12*100;
			   Fifth_Predict[user_count][t].id = id13*100;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].id = id10*100;
			   Second_Predict[user_count][t].id = id11*100;
			   Third_Predict[user_count][t].id = id12*100;
			   Fourth_Predict[user_count][t].id = id13*100;
			   Fifth_Predict[user_count][t].id = id14*100;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].id = id11*100;
			   Second_Predict[user_count][t].id = id12*100;
			   Third_Predict[user_count][t].id = id13*100;
			   Fourth_Predict[user_count][t].id = id14*100;
			   Fifth_Predict[user_count][t].id = id15*100;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].id = id12*100;
			   Second_Predict[user_count][t].id = id13*100;
			   Third_Predict[user_count][t].id = id14*100;
			   Fourth_Predict[user_count][t].id = id15*100;
			   Fifth_Predict[user_count][t].id = id16*100;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].id = id13*100;
			   Second_Predict[user_count][t].id = id14*100;
			   Third_Predict[user_count][t].id = id15*100;
			   Fourth_Predict[user_count][t].id = id16*100;
			   Fifth_Predict[user_count][t].id = id17*100;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].id = id14*100;
			   Second_Predict[user_count][t].id = id15*100;
			   Third_Predict[user_count][t].id = id16*100;
			   Fourth_Predict[user_count][t].id = id17*100;
			   Fifth_Predict[user_count][t].id = id18*100;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].id = id15*100;
			   Second_Predict[user_count][t].id = id16*100;
			   Third_Predict[user_count][t].id = id17*100;
			   Fourth_Predict[user_count][t].id = id18*100;
			   Fifth_Predict[user_count][t].id = id19*100;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].id = id16*100;
			   Second_Predict[user_count][t].id = id17*100;
			   Third_Predict[user_count][t].id = id18*100;
			   Fourth_Predict[user_count][t].id = id19*100;
			   Fifth_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].id = id17*100;
			   Second_Predict[user_count][t].id = id18*100;
			   Third_Predict[user_count][t].id = id19*100;
			   Fourth_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].id = id18*100;
			   Second_Predict[user_count][t].id = id19*100;
			   Third_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].id = id19*100;
			   Second_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].id = id20*100;
		   }
		   if (t == 20)
		   {
			   //First_Predict[user_count][t].b = b2;
			   //Second_Predict[user_count][t].b = b3;
			   //Third_Predict[user_count][t].b = b4;
			   //Fourth_Predict[user_count][t].b = b5;
			   //Fifth_Predict[user_count][t].b = b6;
		   }

		   user_count++;
	   }
	   t++;
   }

   /* md */
   t = 1;

   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("results_Sample_Data_60_md_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;


	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double md1 = 0;
		   geek >> md1;

		   stringstream geek2(two);
		   double md2 = 0;
		   geek2 >> md2;


		   stringstream geek3(three);
		   double md3 = 0;
		   geek3 >> md3;

		   stringstream geek4(four);
		   double md4 = 0;
		   geek4 >> md4;

		   stringstream geek5(five);
		   double md5 = 0;
		   geek5 >> md5;

		   stringstream geek6(six);
		   double md6 = 0;
		   geek6 >> md6;

		   stringstream geek7(seven);
		   double md7 = 0;
		   geek7 >> md7;

		   stringstream geek8(eight);
		   double md8 = 0;
		   geek8 >> md8;

		   stringstream geek9(nine);
		   double md9 = 0;
		   geek9 >> md9;

		   stringstream geek10(ten);
		   double md10 = 0;
		   geek10 >> md10;

		   stringstream geek11(eleven);
		   double md11 = 0;
		   geek11 >> md11;

		   stringstream geek12(twelve);
		   double md12 = 0;
		   geek12 >> md12;

		   stringstream geek13(thirteen);
		   double md13 = 0;
		   geek13 >> md13;

		   stringstream geek14(fourteen);
		   double md14 = 0;
		   geek14 >> md14;

		   stringstream geek15(fifteen);
		   double md15 = 0;
		   geek15 >> md15;

		   stringstream geek16(sixteen);
		   double md16 = 0;
		   geek16 >> md16;

		   stringstream geek17(seventeen);
		   double md17 = 0;
		   geek17 >> md17;

		   stringstream geek18(eighteen);
		   double md18 = 0;
		   geek18 >> md18;

		   stringstream geek19(nineteen);
		   double md19 = 0;
		   geek19 >> md19;

		   stringstream geek20(twenty);
		   double md20 = 0;
		   geek20 >> md20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].md = md2*100;
			   Second_Predict[user_count][t].md = md3*100;
			   Third_Predict[user_count][t].md = md4*100;
			   Fourth_Predict[user_count][t].md = md5*100;
			   Fifth_Predict[user_count][t].md = md6*100;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].md = md3*100;
			   Second_Predict[user_count][t].md = md4*100;
			   Third_Predict[user_count][t].md = md5*100;
			   Fourth_Predict[user_count][t].md = md6*100;
			   Fifth_Predict[user_count][t].md = md7*100;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].md = md4*100;
			   Second_Predict[user_count][t].md = md5*100;
			   Third_Predict[user_count][t].md = md6*100;
			   Fourth_Predict[user_count][t].md = md7*100;
			   Fifth_Predict[user_count][t].md = md8*100;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].md = md5*100;
			   Second_Predict[user_count][t].md = md6*100;
			   Third_Predict[user_count][t].md = md7*100;
			   Fourth_Predict[user_count][t].md = md8*100;
			   Fifth_Predict[user_count][t].md = md9*100;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].md = md6*100;
			   Second_Predict[user_count][t].md = md7*100;
			   Third_Predict[user_count][t].md = md8*100;
			   Fourth_Predict[user_count][t].md = md9*100;
			   Fifth_Predict[user_count][t].md = md10*100;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].md = md7*100;
			   Second_Predict[user_count][t].md = md8*100;
			   Third_Predict[user_count][t].md = md9*100;
			   Fourth_Predict[user_count][t].md = md10*100;
			   Fifth_Predict[user_count][t].md = md11*100;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].md = md8*100;
			   Second_Predict[user_count][t].md = md9*100;
			   Third_Predict[user_count][t].md = md10*100;
			   Fourth_Predict[user_count][t].md = md11*100;
			   Fifth_Predict[user_count][t].md = md12*100;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].md = md9*100;
			   Second_Predict[user_count][t].md = md10*100;
			   Third_Predict[user_count][t].md = md11*100;
			   Fourth_Predict[user_count][t].md = md12*100;
			   Fifth_Predict[user_count][t].md = md13*100;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].md = md10*100;
			   Second_Predict[user_count][t].md = md11*100;
			   Third_Predict[user_count][t].md = md12*100;
			   Fourth_Predict[user_count][t].md = md13*100;
			   Fifth_Predict[user_count][t].md = md14*100;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].md = md11*100;
			   Second_Predict[user_count][t].md = md12*100;
			   Third_Predict[user_count][t].md = md13*100;
			   Fourth_Predict[user_count][t].md = md14*100;
			   Fifth_Predict[user_count][t].md = md15*100;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].md = md12*100;
			   Second_Predict[user_count][t].md = md13*100;
			   Third_Predict[user_count][t].md = md14*100;
			   Fourth_Predict[user_count][t].md = md15*100;
			   Fifth_Predict[user_count][t].md = md16*100;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].md = md13*100;
			   Second_Predict[user_count][t].md = md14*100;
			   Third_Predict[user_count][t].md = md15*100;
			   Fourth_Predict[user_count][t].md = md16*100;
			   Fifth_Predict[user_count][t].md = md17*100;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].md = md14*100;
			   Second_Predict[user_count][t].md = md15*100;
			   Third_Predict[user_count][t].md = md16*100;
			   Fourth_Predict[user_count][t].md = md17*100;
			   Fifth_Predict[user_count][t].md = md18*100;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].md = md15*100;
			   Second_Predict[user_count][t].md = md16*100;
			   Third_Predict[user_count][t].md = md17*100;
			   Fourth_Predict[user_count][t].md = md18*100;
			   Fifth_Predict[user_count][t].md = md19*100;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].md = md16*100;
			   Second_Predict[user_count][t].md = md17*100;
			   Third_Predict[user_count][t].md = md18*100;
			   Fourth_Predict[user_count][t].md = md19*100;
			   Fifth_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].md = md17*100;
			   Second_Predict[user_count][t].md = md18*100;
			   Third_Predict[user_count][t].md = md19*100;
			   Fourth_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].md = md18*100;
			   Second_Predict[user_count][t].md = md19*100;
			   Third_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].md = md19*100;
			   Second_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].md = md20*100;
		   }
		   if (t == 20)
		   {
			   //First_Predict[user_count][t].b = b2;
			   //Second_Predict[user_count][t].b = b3;
			   //Third_Predict[user_count][t].b = b4;
			   //Fourth_Predict[user_count][t].b = b5;
			   //Fifth_Predict[user_count][t].b = b6;
		   }

		   user_count++;
	   }
	   t++;
   }

   /* p */
   t = 1;
   
   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("results_Sample_Data_60_p_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;

	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double p1 = 0;
		   geek >> p1;

		   stringstream geek2(two);
		   double p2 = 0;
		   geek2 >> p2;


		   stringstream geek3(three);
		   double p3 = 0;
		   geek3 >> p3;

		   stringstream geek4(four);
		   double p4 = 0;
		   geek4 >> p4;

		   stringstream geek5(five);
		   double p5 = 0;
		   geek5 >> p5;

		   stringstream geek6(six);
		   double p6 = 0;
		   geek6 >> p6;

		   stringstream geek7(seven);
		   double p7 = 0;
		   geek7 >> p7;

		   stringstream geek8(eight);
		   double p8 = 0;
		   geek8 >> p8;

		   stringstream geek9(nine);
		   double p9 = 0;
		   geek9 >> p9;

		   stringstream geek10(ten);
		   double p10 = 0;
		   geek10 >> p10;

		   stringstream geek11(eleven);
		   double p11 = 0;
		   geek11 >> p11;

		   stringstream geek12(twelve);
		   double p12 = 0;
		   geek12 >> p12;

		   stringstream geek13(thirteen);
		   double p13 = 0;
		   geek13 >> p13;

		   stringstream geek14(fourteen);
		   double p14 = 0;
		   geek14 >> p14;

		   stringstream geek15(fifteen);
		   double p15 = 0;
		   geek15 >> p15;

		   stringstream geek16(sixteen);
		   double p16 = 0;
		   geek16 >> p16;

		   stringstream geek17(seventeen);
		   double p17 = 0;
		   geek17 >> p17;

		   stringstream geek18(eighteen);
		   double p18 = 0;
		   geek18 >> p18;

		   stringstream geek19(nineteen);
		   double p19 = 0;
		   geek19 >> p19;

		   stringstream geek20(twenty);
		   double p20 = 0;
		   geek20 >> p20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].p = p2;
			   Second_Predict[user_count][t].p = p3;
			   Third_Predict[user_count][t].p = p4;
			   Fourth_Predict[user_count][t].p = p5;
			   Fifth_Predict[user_count][t].p = p6;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].p = p3;
			   Second_Predict[user_count][t].p = p4;
			   Third_Predict[user_count][t].p = p5;
			   Fourth_Predict[user_count][t].p = p6;
			   Fifth_Predict[user_count][t].p = p7;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].p = p4;
			   Second_Predict[user_count][t].p = p5;
			   Third_Predict[user_count][t].p = p6;
			   Fourth_Predict[user_count][t].p = p7;
			   Fifth_Predict[user_count][t].p = p8;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].p = p5;
			   Second_Predict[user_count][t].p = p6;
			   Third_Predict[user_count][t].p = p7;
			   Fourth_Predict[user_count][t].p = p8;
			   Fifth_Predict[user_count][t].p = p9;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].p = p6;
			   Second_Predict[user_count][t].p = p7;
			   Third_Predict[user_count][t].p = p8;
			   Fourth_Predict[user_count][t].p = p9;
			   Fifth_Predict[user_count][t].p = p10;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].p = p7;
			   Second_Predict[user_count][t].p = p8;
			   Third_Predict[user_count][t].p = p9;
			   Fourth_Predict[user_count][t].p = p10;
			   Fifth_Predict[user_count][t].p = p11;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].p = p8;
			   Second_Predict[user_count][t].p = p9;
			   Third_Predict[user_count][t].p = p10;
			   Fourth_Predict[user_count][t].p = p11;
			   Fifth_Predict[user_count][t].p = p12;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].p = p9;
			   Second_Predict[user_count][t].p = p10;
			   Third_Predict[user_count][t].p = p11;
			   Fourth_Predict[user_count][t].p = p12;
			   Fifth_Predict[user_count][t].p = p13;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].p = p10;
			   Second_Predict[user_count][t].p = p11;
			   Third_Predict[user_count][t].p = p12;
			   Fourth_Predict[user_count][t].p = p13;
			   Fifth_Predict[user_count][t].p = p14;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].p = p11;
			   Second_Predict[user_count][t].p = p12;
			   Third_Predict[user_count][t].p = p13;
			   Fourth_Predict[user_count][t].p = p14;
			   Fifth_Predict[user_count][t].p = p15;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].p = p12;
			   Second_Predict[user_count][t].p = p13;
			   Third_Predict[user_count][t].p = p14;
			   Fourth_Predict[user_count][t].p = p15;
			   Fifth_Predict[user_count][t].p = p16;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].p = p13;
			   Second_Predict[user_count][t].p = p14;
			   Third_Predict[user_count][t].p = p15;
			   Fourth_Predict[user_count][t].p = p16;
			   Fifth_Predict[user_count][t].p = p17;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].p = p14;
			   Second_Predict[user_count][t].p = p15;
			   Third_Predict[user_count][t].p = p16;
			   Fourth_Predict[user_count][t].p = p17;
			   Fifth_Predict[user_count][t].p = p18;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].p = p15;
			   Second_Predict[user_count][t].p = p16;
			   Third_Predict[user_count][t].p = p17;
			   Fourth_Predict[user_count][t].p = p18;
			   Fifth_Predict[user_count][t].p = p19;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].p = p16;
			   Second_Predict[user_count][t].p = p17;
			   Third_Predict[user_count][t].p = p18;
			   Fourth_Predict[user_count][t].p = p19;
			   Fifth_Predict[user_count][t].p = p20;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].p = p17;
			   Second_Predict[user_count][t].p = p18;
			   Third_Predict[user_count][t].p = p19;
			   Fourth_Predict[user_count][t].p = p20;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].p = p18;
			   Second_Predict[user_count][t].p = p19;
			   Third_Predict[user_count][t].p = p20;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].p = p19;
			   Second_Predict[user_count][t].p = p20;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].p = p20;
		   }
		   if (t == 20)
		   {
	
		   }

		   user_count++;
	   }
	   t++;
   }

   /* w */
   t = 1;
   
   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("results_Sample_Data_60_w_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;


	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double w1 = 0;
		   geek >> w1;

		   stringstream geek2(two);
		   double w2 = 0;
		   geek2 >> w2;

		   stringstream geek3(three);
		   double w3 = 0;
		   geek3 >> w3;

		   stringstream geek4(four);
		   double w4 = 0;
		   geek4 >> w4;

		   stringstream geek5(five);
		   double w5 = 0;
		   geek5 >> w5;

		   stringstream geek6(six);
		   double w6 = 0;
		   geek6 >> w6;

		   stringstream geek7(seven);
		   double w7 = 0;
		   geek7 >> w7;

		   stringstream geek8(eight);
		   double w8 = 0;
		   geek8 >> w8;

		   stringstream geek9(nine);
		   double w9 = 0;
		   geek9 >> w9;

		   stringstream geek10(ten);
		   double w10 = 0;
		   geek10 >> w10;

		   stringstream geek11(eleven);
		   double w11 = 0;
		   geek11 >> w11;

		   stringstream geek12(twelve);
		   double w12 = 0;
		   geek12 >> w12;

		   stringstream geek13(thirteen);
		   double w13 = 0;
		   geek13 >> w13;

		   stringstream geek14(fourteen);
		   double w14 = 0;
		   geek14 >> w14;

		   stringstream geek15(fifteen);
		   double w15 = 0;
		   geek15 >> w15;

		   stringstream geek16(sixteen);
		   double w16 = 0;
		   geek16 >> w16;

		   stringstream geek17(seventeen);
		   double w17 = 0;
		   geek17 >> w17;

		   stringstream geek18(eighteen);
		   double w18 = 0;
		   geek18 >> w18;

		   stringstream geek19(nineteen);
		   double w19 = 0;
		   geek19 >> w19;

		   stringstream geek20(twenty);
		   double w20 = 0;
		   geek20 >> w20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].w = w2*100;
			   Second_Predict[user_count][t].w = w3*100;
			   Third_Predict[user_count][t].w = w4*100;
			   Fourth_Predict[user_count][t].w = w5*100;
			   Fifth_Predict[user_count][t].w = w6*100;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].w = w3*100;
			   Second_Predict[user_count][t].w = w4*100;
			   Third_Predict[user_count][t].w = w5*100;
			   Fourth_Predict[user_count][t].w = w6*100;
			   Fifth_Predict[user_count][t].w = w7*100;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].w = w4*100;
			   Second_Predict[user_count][t].w = w5*100;
			   Third_Predict[user_count][t].w = w6*100;
			   Fourth_Predict[user_count][t].w = w7*100;
			   Fifth_Predict[user_count][t].w = w8*100;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].w = w5*100;
			   Second_Predict[user_count][t].w = w6*100;
			   Third_Predict[user_count][t].w = w7*100;
			   Fourth_Predict[user_count][t].w = w8*100;
			   Fifth_Predict[user_count][t].w = w9*100;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].w = w6*100;
			   Second_Predict[user_count][t].w = w7*100;
			   Third_Predict[user_count][t].w = w8*100;
			   Fourth_Predict[user_count][t].w = w9*100;
			   Fifth_Predict[user_count][t].w = w10*100;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].w = w7*100;
			   Second_Predict[user_count][t].w = w8*100;
			   Third_Predict[user_count][t].w = w9*100;
			   Fourth_Predict[user_count][t].w = w10*100;
			   Fifth_Predict[user_count][t].w = w11*100;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].w = w8*100;
			   Second_Predict[user_count][t].w = w9*100;
			   Third_Predict[user_count][t].w = w10*100;
			   Fourth_Predict[user_count][t].w = w11*100;
			   Fifth_Predict[user_count][t].w = w12*100;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].w = w9*100;
			   Second_Predict[user_count][t].w = w10*100;
			   Third_Predict[user_count][t].w = w11*100;
			   Fourth_Predict[user_count][t].w = w12*100;
			   Fifth_Predict[user_count][t].w = w13*100;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].w = w10*100;
			   Second_Predict[user_count][t].w = w11*100;
			   Third_Predict[user_count][t].w = w12*100;
			   Fourth_Predict[user_count][t].w = w13*100;
			   Fifth_Predict[user_count][t].w = w14*100;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].w = w11*100;
			   Second_Predict[user_count][t].w = w12*100;
			   Third_Predict[user_count][t].w = w13*100;
			   Fourth_Predict[user_count][t].w = w14*100;
			   Fifth_Predict[user_count][t].w = w15*100;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].w = w12*100;
			   Second_Predict[user_count][t].w = w13*100;
			   Third_Predict[user_count][t].w = w14*100;
			   Fourth_Predict[user_count][t].w = w15*100;
			   Fifth_Predict[user_count][t].w = w16*100;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].w = w13*100;
			   Second_Predict[user_count][t].w = w14*100;
			   Third_Predict[user_count][t].w = w15*100;
			   Fourth_Predict[user_count][t].w = w16*100;
			   Fifth_Predict[user_count][t].w = w17*100;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].w = w14*100;
			   Second_Predict[user_count][t].w = w15*100;
			   Third_Predict[user_count][t].w = w16*100;
			   Fourth_Predict[user_count][t].w = w17*100;
			   Fifth_Predict[user_count][t].w = w18*100;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].w = w15*100;
			   Second_Predict[user_count][t].w = w16*100;
			   Third_Predict[user_count][t].w = w17*100;
			   Fourth_Predict[user_count][t].w = w18*100;
			   Fifth_Predict[user_count][t].w = w19*100;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].w = w16*100;
			   Second_Predict[user_count][t].w = w17*100;
			   Third_Predict[user_count][t].w = w18*100;
			   Fourth_Predict[user_count][t].w = w19*100;
			   Fifth_Predict[user_count][t].w = w20*100;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].w = w17*100;
			   Second_Predict[user_count][t].w = w18*100;
			   Third_Predict[user_count][t].w = w19*100;
			   Fourth_Predict[user_count][t].w = w20*100;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].w = w18*100;
			   Second_Predict[user_count][t].w = w19*100;
			   Third_Predict[user_count][t].w = w20*100;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].w = w19*100;
			   Second_Predict[user_count][t].w = w20*100;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].w = w20*100;

		   }
		   if (t == 20)
		   {
			
		   }

		   user_count++;
	   }
	   t++;
   }

   /* x */
   t = 1;

   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("results_Sample_Data_60_x_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;


	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double x1 = 0;
		   geek >> x1;

		   stringstream geek2(two);
		   double x2 = 0;
		   geek2 >> x2;


		   stringstream geek3(three);
		   double x3 = 0;
		   geek3 >> x3;

		   stringstream geek4(four);
		   double x4 = 0;
		   geek4 >> x4;

		   stringstream geek5(five);
		   double x5 = 0;
		   geek5 >> x5;

		   stringstream geek6(six);
		   double x6 = 0;
		   geek6 >> x6;

		   stringstream geek7(seven);
		   double x7 = 0;
		   geek7 >> x7;

		   stringstream geek8(eight);
		   double x8 = 0;
		   geek8 >> x8;

		   stringstream geek9(nine);
		   double x9 = 0;
		   geek9 >> x9;

		   stringstream geek10(ten);
		   double x10 = 0;
		   geek10 >> x10;

		   stringstream geek11(eleven);
		   double x11 = 0;
		   geek11 >> x11;

		   stringstream geek12(twelve);
		   double x12 = 0;
		   geek12 >> x12;

		   stringstream geek13(thirteen);
		   double x13 = 0;
		   geek13 >> x13;

		   stringstream geek14(fourteen);
		   double x14 = 0;
		   geek14 >> x14;

		   stringstream geek15(fifteen);
		   double x15 = 0;
		   geek15 >> x15;

		   stringstream geek16(sixteen);
		   double x16 = 0;
		   geek16 >> x16;

		   stringstream geek17(seventeen);
		   double x17 = 0;
		   geek17 >> x17;

		   stringstream geek18(eighteen);
		   double x18 = 0;
		   geek18 >> x18;

		   stringstream geek19(nineteen);
		   double x19 = 0;
		   geek19 >> x19;

		   stringstream geek20(twenty);
		   double x20 = 0;
		   geek20 >> x20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].x = x2;
			   Second_Predict[user_count][t].x = x3;
			   Third_Predict[user_count][t].x = x4;
			   Fourth_Predict[user_count][t].x = x5;
			   Fifth_Predict[user_count][t].x = x6;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].x = x3;
			   Second_Predict[user_count][t].x = x4;
			   Third_Predict[user_count][t].x = x5;
			   Fourth_Predict[user_count][t].x = x6;
			   Fifth_Predict[user_count][t].x = x7;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].x = x4;
			   Second_Predict[user_count][t].x = x5;
			   Third_Predict[user_count][t].x = x6;
			   Fourth_Predict[user_count][t].x = x7;
			   Fifth_Predict[user_count][t].x = x8;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].x = x5;
			   Second_Predict[user_count][t].x = x6;
			   Third_Predict[user_count][t].x = x7;
			   Fourth_Predict[user_count][t].x = x8;
			   Fifth_Predict[user_count][t].x = x9;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].x = x6;
			   Second_Predict[user_count][t].x = x7;
			   Third_Predict[user_count][t].x = x8;
			   Fourth_Predict[user_count][t].x = x9;
			   Fifth_Predict[user_count][t].x = x10;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].x = x7;
			   Second_Predict[user_count][t].x = x8;
			   Third_Predict[user_count][t].x = x9;
			   Fourth_Predict[user_count][t].x = x10;
			   Fifth_Predict[user_count][t].x = x11;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].x = x8;
			   Second_Predict[user_count][t].x = x9;
			   Third_Predict[user_count][t].x = x10;
			   Fourth_Predict[user_count][t].x = x11;
			   Fifth_Predict[user_count][t].x = x12;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].x = x9;
			   Second_Predict[user_count][t].x = x10;
			   Third_Predict[user_count][t].x = x11;
			   Fourth_Predict[user_count][t].x = x12;
			   Fifth_Predict[user_count][t].x = x13;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].x = x10;
			   Second_Predict[user_count][t].x = x11;
			   Third_Predict[user_count][t].x = x12;
			   Fourth_Predict[user_count][t].x = x13;
			   Fifth_Predict[user_count][t].x = x14;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].x = x11;
			   Second_Predict[user_count][t].x = x12;
			   Third_Predict[user_count][t].x = x13;
			   Fourth_Predict[user_count][t].x = x14;
			   Fifth_Predict[user_count][t].x = x15;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].x = x12;
			   Second_Predict[user_count][t].x = x13;
			   Third_Predict[user_count][t].x = x14;
			   Fourth_Predict[user_count][t].x = x15;
			   Fifth_Predict[user_count][t].x = x16;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].x = x13;
			   Second_Predict[user_count][t].x = x14;
			   Third_Predict[user_count][t].x = x15;
			   Fourth_Predict[user_count][t].x = x16;
			   Fifth_Predict[user_count][t].x = x17;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].x = x14;
			   Second_Predict[user_count][t].x = x15;
			   Third_Predict[user_count][t].x = x16;
			   Fourth_Predict[user_count][t].x = x17;
			   Fifth_Predict[user_count][t].x = x18;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].x = x15;
			   Second_Predict[user_count][t].x = x16;
			   Third_Predict[user_count][t].x = x17;
			   Fourth_Predict[user_count][t].x = x18;
			   Fifth_Predict[user_count][t].x = x19;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].x = x16;
			   Second_Predict[user_count][t].x = x17;
			   Third_Predict[user_count][t].x = x18;
			   Fourth_Predict[user_count][t].x = x19;
			   Fifth_Predict[user_count][t].x = x20;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].x = x17;
			   Second_Predict[user_count][t].x = x18;
			   Third_Predict[user_count][t].x = x19;
			   Fourth_Predict[user_count][t].x = x20;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].x = x18;
			   Second_Predict[user_count][t].x = x19;
			   Third_Predict[user_count][t].x = x20;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].x = x19;
			   Second_Predict[user_count][t].x = x20;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].x = x20;
		   }
		   if (t == 20)
		   {
			   
		   }

		   user_count++;
	   }
	   t++;
   }

   /* y */
   t = 1;

   while (t <= 19)
   {
	   stringstream out;
	   out << t;
	   ifstream data("results_Sample_Data_60_y_" + out.str() + ".csv");

	   if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	   string one;
	   string two;
	   string three;
	   string four;
	   string five;
	   string six;
	   string seven;
	   string eight;
	   string nine;
	   string ten;
	   string eleven;
	   string twelve;
	   string thirteen;
	   string fourteen;
	   string fifteen;
	   string sixteen;
	   string seventeen;
	   string eighteen;
	   string nineteen;
	   string twenty;

	   int user_count = 0;


	   while (data.good() && user_count<user_total)
	   {

		   getline(data, one, ',');
		   getline(data, two, ',');
		   getline(data, three, ',');
		   getline(data, four, ',');
		   getline(data, five, ',');
		   getline(data, six, ',');
		   getline(data, seven, ',');
		   getline(data, eight, ',');
		   getline(data, nine, ',');
		   getline(data, ten, ',');
		   getline(data, eleven, ',');
		   getline(data, twelve, ',');
		   getline(data, thirteen, ',');
		   getline(data, fourteen, ',');
		   getline(data, fifteen, ',');
		   getline(data, sixteen, ',');
		   getline(data, seventeen, ',');
		   getline(data, eighteen, ',');
		   getline(data, nineteen, ',');
		   getline(data, twenty, '\n');


		   /*  ---------------------------   */
		   stringstream geek(one);
		   double y1 = 0;
		   geek >> y1;

		   stringstream geek2(two);
		   double y2 = 0;
		   geek2 >> y2;


		   stringstream geek3(three);
		   double y3 = 0;
		   geek3 >> y3;

		   stringstream geek4(four);
		   double y4 = 0;
		   geek4 >> y4;

		   stringstream geek5(five);
		   double y5 = 0;
		   geek5 >> y5;

		   stringstream geek6(six);
		   double y6 = 0;
		   geek6 >> y6;

		   stringstream geek7(seven);
		   double y7 = 0;
		   geek7 >> y7;

		   stringstream geek8(eight);
		   double y8 = 0;
		   geek8 >> y8;

		   stringstream geek9(nine);
		   double y9 = 0;
		   geek9 >> y9;

		   stringstream geek10(ten);
		   double y10 = 0;
		   geek10 >> y10;

		   stringstream geek11(eleven);
		   double y11 = 0;
		   geek11 >> y11;

		   stringstream geek12(twelve);
		   double y12 = 0;
		   geek12 >> y12;

		   stringstream geek13(thirteen);
		   double y13 = 0;
		   geek13 >> y13;

		   stringstream geek14(fourteen);
		   double y14 = 0;
		   geek14 >> y14;

		   stringstream geek15(fifteen);
		   double y15 = 0;
		   geek15 >> y15;

		   stringstream geek16(sixteen);
		   double y16 = 0;
		   geek16 >> y16;

		   stringstream geek17(seventeen);
		   double y17 = 0;
		   geek17 >> y17;

		   stringstream geek18(eighteen);
		   double y18 = 0;
		   geek18 >> y18;

		   stringstream geek19(nineteen);
		   double y19 = 0;
		   geek19 >> y19;

		   stringstream geek20(twenty);
		   double y20 = 0;
		   geek20 >> y20;

		   if (t == 1)
		   {
			   First_Predict[user_count][t].y = y2;
			   Second_Predict[user_count][t].y = y3;
			   Third_Predict[user_count][t].y = y4;
			   Fourth_Predict[user_count][t].y = y5;
			   Fifth_Predict[user_count][t].y = y6;
		   }
		   if (t == 2)
		   {
			   First_Predict[user_count][t].y = y3;
			   Second_Predict[user_count][t].y = y4;
			   Third_Predict[user_count][t].y = y5;
			   Fourth_Predict[user_count][t].y = y6;
			   Fifth_Predict[user_count][t].y = y7;
		   }
		   if (t == 3)
		   {
			   First_Predict[user_count][t].y = y4;
			   Second_Predict[user_count][t].y = y5;
			   Third_Predict[user_count][t].y = y6;
			   Fourth_Predict[user_count][t].y = y7;
			   Fifth_Predict[user_count][t].y = y8;
		   }
		   if (t == 4)
		   {
			   First_Predict[user_count][t].y = y5;
			   Second_Predict[user_count][t].y = y6;
			   Third_Predict[user_count][t].y = y7;
			   Fourth_Predict[user_count][t].y = y8;
			   Fifth_Predict[user_count][t].y = y9;
		   }
		   if (t == 5)
		   {
			   First_Predict[user_count][t].y = y6;
			   Second_Predict[user_count][t].y = y7;
			   Third_Predict[user_count][t].y = y8;
			   Fourth_Predict[user_count][t].y = y9;
			   Fifth_Predict[user_count][t].y = y10;
		   }
		   if (t == 6)
		   {
			   First_Predict[user_count][t].y = y7;
			   Second_Predict[user_count][t].y = y8;
			   Third_Predict[user_count][t].y = y9;
			   Fourth_Predict[user_count][t].y = y10;
			   Fifth_Predict[user_count][t].y = y11;
		   }
		   if (t == 7)
		   {
			   First_Predict[user_count][t].y = y8;
			   Second_Predict[user_count][t].y = y9;
			   Third_Predict[user_count][t].y = y10;
			   Fourth_Predict[user_count][t].y = y11;
			   Fifth_Predict[user_count][t].y = y12;
		   }
		   if (t == 8)
		   {
			   First_Predict[user_count][t].y = y9;
			   Second_Predict[user_count][t].y = y10;
			   Third_Predict[user_count][t].y = y11;
			   Fourth_Predict[user_count][t].y = y12;
			   Fifth_Predict[user_count][t].y = y13;
		   }
		   if (t == 9)
		   {
			   First_Predict[user_count][t].y = y10;
			   Second_Predict[user_count][t].y = y11;
			   Third_Predict[user_count][t].y = y12;
			   Fourth_Predict[user_count][t].y = y13;
			   Fifth_Predict[user_count][t].y = y14;
		   }
		   if (t == 10)
		   {
			   First_Predict[user_count][t].y = y11;
			   Second_Predict[user_count][t].y = y12;
			   Third_Predict[user_count][t].y = y13;
			   Fourth_Predict[user_count][t].y = y14;
			   Fifth_Predict[user_count][t].y = y15;
		   }
		   if (t == 11)
		   {
			   First_Predict[user_count][t].y = y12;
			   Second_Predict[user_count][t].y = y13;
			   Third_Predict[user_count][t].y = y14;
			   Fourth_Predict[user_count][t].y = y15;
			   Fifth_Predict[user_count][t].y = y16;
		   }
		   if (t == 12)
		   {
			   First_Predict[user_count][t].y = y13;
			   Second_Predict[user_count][t].y = y14;
			   Third_Predict[user_count][t].y = y15;
			   Fourth_Predict[user_count][t].y = y16;
			   Fifth_Predict[user_count][t].y = y17;
		   }
		   if (t == 13)
		   {
			   First_Predict[user_count][t].y = y14;
			   Second_Predict[user_count][t].y = y15;
			   Third_Predict[user_count][t].y = y16;
			   Fourth_Predict[user_count][t].y = y17;
			   Fifth_Predict[user_count][t].y = y18;
		   }
		   if (t == 14)
		   {
			   First_Predict[user_count][t].y = y15;
			   Second_Predict[user_count][t].y = y16;
			   Third_Predict[user_count][t].y = y17;
			   Fourth_Predict[user_count][t].y = y18;
			   Fifth_Predict[user_count][t].y = y19;
		   }
		   if (t == 15)
		   {
			   First_Predict[user_count][t].y = y16;
			   Second_Predict[user_count][t].y = y17;
			   Third_Predict[user_count][t].y = y18;
			   Fourth_Predict[user_count][t].y = y19;
			   Fifth_Predict[user_count][t].y = y20;
		   }
		   if (t == 16)
		   {
			   First_Predict[user_count][t].y = y17;
			   Second_Predict[user_count][t].y = y18;
			   Third_Predict[user_count][t].y = y19;
			   Fourth_Predict[user_count][t].y = y20;
		   }
		   if (t == 17)
		   {
			   First_Predict[user_count][t].y = y18;
			   Second_Predict[user_count][t].y = y19;
			   Third_Predict[user_count][t].y = y20;
		   }
		   if (t == 18)
		   {
			   First_Predict[user_count][t].y = y19;
			   Second_Predict[user_count][t].y = y20;
		   }
		   if (t == 19)
		   {
			   First_Predict[user_count][t].y = y20;
		   }
		   if (t == 20)
		   {
			  
		   }

		   user_count++;
	   }
	   t++;
   }
   
}

void ReadReal()
{
	/* new dataset --- beginning-----  */
	for (int i = 0;i<user_total;i++)
	{
		users[i] = new Users[200];
		R[i] = new int[200];
		count_time[i] = 200;

	}
  {
	ifstream data("Sample_Data_60+21_b.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double b1 = 0;
		geek >> b1;
		users[user_count][1].b = b1;

		stringstream geek2(two);
		double b2 = 0;
		geek2 >> b2;
		users[user_count][2].b = b2;

		stringstream geek3(three);
		double b3 = 0;
		geek3 >> b3;
		users[user_count][3].b = b3;

		stringstream geek4(four);
		double b4 = 0;
		geek4 >> b4;
		users[user_count][4].b = b4;

		stringstream geek5(five);
		double b5 = 0;
		geek5 >> b5;
		users[user_count][5].b = b5;

		stringstream geek6(six);
		double b6 = 0;
		geek6 >> b6;
		users[user_count][6].b = b6;

		stringstream geek7(seven);
		double b7 = 0;
		geek7 >> b7;
		users[user_count][7].b = b7;

		stringstream geek8(eight);
		double b8 = 0;
		geek8 >> b8;
		users[user_count][8].b = b8;

		stringstream geek9(nine);
		double b9 = 0;
		geek9 >> b9;
		users[user_count][9].b = b9;

		stringstream geek10(ten);
		double b10 = 0;
		geek10 >> b10;
		users[user_count][10].b = b10;

		stringstream geek11(eleven);
		double b11 = 0;
		geek11 >> b11;
		users[user_count][11].b = b11;

		stringstream geek12(twelve);
		double b12 = 0;
		geek12 >> b12;
		users[user_count][12].b = b12;

		stringstream geek13(thirteen);
		double b13 = 0;
		geek13 >> b13;
		users[user_count][13].b = b13;

		stringstream geek14(fourteen);
		double b14 = 0;
		geek14 >> b14;
		users[user_count][14].b = b14;

		stringstream geek15(fifteen);
		double b15 = 0;
		geek15 >> b15;
		users[user_count][15].b = b15;

		stringstream geek16(sixteen);
		double b16 = 0;
		geek16 >> b16;
		users[user_count][16].b = b16;

		stringstream geek17(seventeen);
		double b17 = 0;
		geek17 >> b17;
		users[user_count][17].b = b17;

		stringstream geek18(eighteen);
		double b18 = 0;
		geek18 >> b18;
		users[user_count][18].b = b18;

		stringstream geek19(nineteen);
		double b19 = 0;
		geek19 >> b19;
		users[user_count][19].b = b19;

		stringstream geek20(twenty);
		double b20 = 0;
		geek20 >> b20;
		users[user_count][20].b = b20;

		user_count++;
	}

   }
	
   {
	ifstream data("Sample_Data_60+21_id.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double id1 = 0;
		geek >> id1;
		users[user_count][1].id = id1*100;

		stringstream geek2(two);
		double id2 = 0;
		geek2 >> id2;
		users[user_count][2].id = id2*100;

		stringstream geek3(three);
		double id3 = 0;
		geek3 >> id3;
		users[user_count][3].id = id3*100;

		stringstream geek4(four);
		double id4 = 0;
		geek4 >> id4;
		users[user_count][4].id = id4*100;

		stringstream geek5(five);
		double id5 = 0;
		geek5 >> id5;
		users[user_count][5].id = id5*100;

		stringstream geek6(six);
		double id6 = 0;
		geek6 >> id6;
		users[user_count][6].id = id6*100;

		stringstream geek7(seven);
		double id7 = 0;
		geek7 >> id7;
		users[user_count][7].id = id7*100;

		stringstream geek8(eight);
		double id8 = 0;
		geek8 >> id8;
		users[user_count][8].id = id8*100;

		stringstream geek9(nine);
		double id9 = 0;
		geek9 >> id9;
		users[user_count][9].id = id9*100;

		stringstream geek10(ten);
		double id10 = 0;
		geek10 >> id10;
		users[user_count][10].id = id10*100;

		stringstream geek11(eleven);
		double id11 = 0;
		geek11 >> id11;
		users[user_count][11].id = id11*100;

		stringstream geek12(twelve);
		double id12 = 0;
		geek12 >> id12;
		users[user_count][12].id = id12*100;

		stringstream geek13(thirteen);
		double id13 = 0;
		geek13 >> id13;
		users[user_count][13].id = id13*100;

		stringstream geek14(fourteen);
		double id14 = 0;
		geek14 >> id14;
		users[user_count][14].id = id14*100;

		stringstream geek15(fifteen);
		double id15 = 0;
		geek15 >> id15;
		users[user_count][15].id = id15*100;

		stringstream geek16(sixteen);
		double id16 = 0;
		geek16 >> id16;
		users[user_count][16].id = id16*100;

		stringstream geek17(seventeen);
		double id17 = 0;
		geek17 >> id17;
		users[user_count][17].id = id17*100;

		stringstream geek18(eighteen);
		double id18 = 0;
		geek18 >> id18;
		users[user_count][18].id = id18*100;

		stringstream geek19(nineteen);
		double id19 = 0;
		geek19 >> id19;
		users[user_count][19].id = id19*100;

		stringstream geek20(twenty);
		double id20 = 0;
		geek20 >> id20;
		users[user_count][20].id = id20*100;

		user_count++;
	}
	}

	{
	ifstream data("Sample_Data_60+21_md.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double md1 = 0;
		geek >> md1;
		users[user_count][1].md = md1*100;

		stringstream geek2(two);
		double md2 = 0;
		geek2 >> md2;
		users[user_count][2].md = md2*100;

		stringstream geek3(three);
		double md3 = 0;
		geek3 >> md3;
		users[user_count][3].md = md3*100;

		stringstream geek4(four);
		double md4 = 0;
		geek4 >> md4;
		users[user_count][4].md = md4*100;

		stringstream geek5(five);
		double md5 = 0;
		geek5 >> md5;
		users[user_count][5].md = md5*100;

		stringstream geek6(six);
		double md6 = 0;
		geek6 >> md6;
		users[user_count][6].md = md6*100;

		stringstream geek7(seven);
		double md7 = 0;
		geek7 >> md7;
		users[user_count][7].md = md7*100;

		stringstream geek8(eight);
		double md8 = 0;
		geek8 >> md8;
		users[user_count][8].md = md8*100;

		stringstream geek9(nine);
		double md9 = 0;
		geek9 >> md9;
		users[user_count][9].md = md9*100;

		stringstream geek10(ten);
		double md10 = 0;
		geek10 >> md10;
		users[user_count][10].md = md10*100;

		stringstream geek11(eleven);
		double md11 = 0;
		geek11 >> md11;
		users[user_count][11].md = md11*100;

		stringstream geek12(twelve);
		double md12 = 0;
		geek12 >> md12;
		users[user_count][12].md = md12*100;

		stringstream geek13(thirteen);
		double md13 = 0;
		geek13 >> md13;
		users[user_count][13].md = md13*100;

		stringstream geek14(fourteen);
		double md14 = 0;
		geek14 >> md14;
		users[user_count][14].md = md14*100;

		stringstream geek15(fifteen);
		double md15 = 0;
		geek15 >> md15;
		users[user_count][15].md = md15*100;

		stringstream geek16(sixteen);
		double md16 = 0;
		geek16 >> md16;
		users[user_count][16].md = md16*100;

		stringstream geek17(seventeen);
		double md17 = 0;
		geek17 >> md17;
		users[user_count][17].md = md17*100;

		stringstream geek18(eighteen);
		double md18 = 0;
		geek18 >> md18;
		users[user_count][18].md = md18*100;

		stringstream geek19(nineteen);
		double md19 = 0;
		geek19 >> md19;
		users[user_count][19].md = md19*100;

		stringstream geek20(twenty);
		double md20 = 0;
		geek20 >> md20;
		users[user_count][20].md = md20*100;

		user_count++;
	}
	}

	{
	ifstream data("Sample_Data_60+21_w.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double w1 = 0;
		geek >> w1;
		users[user_count][1].w = w1*100;

		stringstream geek2(two);
		double w2 = 0;
		geek2 >> w2;
		users[user_count][2].w = w2*100;

		stringstream geek3(three);
		double w3 = 0;
		geek3 >> w3;
		users[user_count][3].w = w3*100;

		stringstream geek4(four);
		double w4 = 0;
		geek4 >> w4;
		users[user_count][4].w = w4*100;

		stringstream geek5(five);
		double w5 = 0;
		geek5 >> w5;
		users[user_count][5].w = w5*100;

		stringstream geek6(six);
		double w6 = 0;
		geek6 >> w6;
		users[user_count][6].w = w6*100;

		stringstream geek7(seven);
		double w7 = 0;
		geek7 >> w7;
		users[user_count][7].w = w7*100;

		stringstream geek8(eight);
		double w8 = 0;
		geek8 >> w8;
		users[user_count][8].w = w8*100;

		stringstream geek9(nine);
		double w9 = 0;
		geek9 >> w9;
		users[user_count][9].w = w9*100;

		stringstream geek10(ten);
		double w10 = 0;
		geek10 >> w10;
		users[user_count][10].w = w10*100;

		stringstream geek11(eleven);
		double w11 = 0;
		geek11 >> w11;
		users[user_count][11].w = w11*100;

		stringstream geek12(twelve);
		double w12 = 0;
		geek12 >> w12;
		users[user_count][12].w = w12*100;

		stringstream geek13(thirteen);
		double w13 = 0;
		geek13 >> w13;
		users[user_count][13].w = w13*100;

		stringstream geek14(fourteen);
		double w14 = 0;
		geek14 >> w14;
		users[user_count][14].w = w14*100;

		stringstream geek15(fifteen);
		double w15 = 0;
		geek15 >> w15;
		users[user_count][15].w = w15*100;

		stringstream geek16(sixteen);
		double w16 = 0;
		geek16 >> w16;
		users[user_count][16].w = w16*100;

		stringstream geek17(seventeen);
		double w17 = 0;
		geek17 >> w17;
		users[user_count][17].w = w17*100;

		stringstream geek18(eighteen);
		double w18 = 0;
		geek18 >> w18;
		users[user_count][18].w = w18*100;

		stringstream geek19(nineteen);
		double w19 = 0;
		geek19 >> w19;
		users[user_count][19].w = w19*100;

		stringstream geek20(twenty);
		double w20 = 0;
		geek20 >> w20;
		users[user_count][20].w = w20*100;

		user_count++;
	}
	}

	{
	ifstream data("Sample_Data_60+21_p.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double p1 = 0;
		geek >> p1;
		users[user_count][1].p = p1;

		stringstream geek2(two);
		double p2 = 0;
		geek2 >> p2;
		users[user_count][2].p = p2;

		stringstream geek3(three);
		double p3 = 0;
		geek3 >> p3;
		users[user_count][3].p = p3;

		stringstream geek4(four);
		double p4 = 0;
		geek4 >> p4;
		users[user_count][4].p = p4;

		stringstream geek5(five);
		double p5 = 0;
		geek5 >> p5;
		users[user_count][5].p = p5;

		stringstream geek6(six);
		double p6 = 0;
		geek6 >> p6;
		users[user_count][6].p = p6;

		stringstream geek7(seven);
		double p7 = 0;
		geek7 >> p7;
		users[user_count][7].p = p7;

		stringstream geek8(eight);
		double p8 = 0;
		geek8 >> p8;
		users[user_count][8].p = p8;

		stringstream geek9(nine);
		double p9 = 0;
		geek9 >> p9;
		users[user_count][9].p = p9;

		stringstream geek10(ten);
		double p10 = 0;
		geek10 >> p10;
		users[user_count][10].p = p10;

		stringstream geek11(eleven);
		double p11 = 0;
		geek11 >> p11;
		users[user_count][11].p = p11;

		stringstream geek12(twelve);
		double p12 = 0;
		geek12 >> p12;
		users[user_count][12].p = p12;

		stringstream geek13(thirteen);
		double p13 = 0;
		geek13 >> p13;
		users[user_count][13].p = p13;

		stringstream geek14(fourteen);
		double p14 = 0;
		geek14 >> p14;
		users[user_count][14].p = p14;

		stringstream geek15(fifteen);
		double p15 = 0;
		geek15 >> p15;
		users[user_count][15].p = p15;

		stringstream geek16(sixteen);
		double p16 = 0;
		geek16 >> p16;
		users[user_count][16].p = p16;

		stringstream geek17(seventeen);
		double p17 = 0;
		geek17 >> p17;
		users[user_count][17].p = p17;

		stringstream geek18(eighteen);
		double p18 = 0;
		geek18 >> p18;
		users[user_count][18].p = p18;

		stringstream geek19(nineteen);
		double p19 = 0;
		geek19 >> p19;
		users[user_count][19].p = p19;

		stringstream geek20(twenty);
		double p20 = 0;
		geek20 >> p20;
		users[user_count][20].p = p20;

		user_count++;
	}
	}

	{
	ifstream data("Sample_Data_60+21_x.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double x1 = 0;
		geek >> x1;
		users[user_count][1].x = x1;
		if (x1 <= min_x)
			min_x = x1;
		if (x1>max_x)
			max_x = x1;

		stringstream geek2(two);
		double x2 = 0;
		geek2 >> x2;
		users[user_count][2].x = x2;
		if (x2 <= min_x)
			min_x = x2;
		if (x2>max_x)
			max_x = x2;
		stringstream geek3(three);
		double x3 = 0;
		geek3 >> x3;
		users[user_count][3].x = x3;
		if (x3 <= min_x)
			min_x = x3;
		if (x3>max_x)
			max_x = x3;

		stringstream geek4(four);
		double x4 = 0;
		geek4 >> x4;
		users[user_count][4].x = x4;
		if (x4 <= min_x)
			min_x = x4;
		if (x4>max_x)
			max_x = x4;

		stringstream geek5(five);
		double x5 = 0;
		geek5 >> x5;
		users[user_count][5].x = x5;
		if (x5 <= min_x)
			min_x = x5;
		if (x5>max_x)
			max_x = x5;
		stringstream geek6(six);
		double x6 = 0;
		geek6 >> x6;
		users[user_count][6].x = x6;
		if (x6 <= min_x)
			min_x = x6;
		if (x6>max_x)
			max_x = x6;
		stringstream geek7(seven);
		double x7 = 0;
		geek7 >> x7;
		users[user_count][7].x = x7;
		if (x7 <= min_x)
			min_x = x7;
		if (x7>max_x)
			max_x = x7;
		stringstream geek8(eight);
		double x8 = 0;
		geek8 >> x8;
		users[user_count][8].x = x8;
		if (x8 <= min_x)
			min_x = x8;
		if (x8>max_x)
			max_x = x8;
		stringstream geek9(nine);
		double x9 = 0;
		geek9 >> x9;
		users[user_count][9].x = x9;
		if (x9 <= min_x)
			min_x = x9;
		if (x9>max_x)
			max_x = x9;
		stringstream geek10(ten);
		double x10 = 0;
		geek10 >> x10;
		users[user_count][10].x = x10;
		if (x10 <= min_x)
			min_x = x10;
		if (x10>max_x)
			max_x = x10;
		stringstream geek11(eleven);
		double x11 = 0;
		geek11 >> x11;
		users[user_count][11].x = x11;
		if (x11 <= min_x)
			min_x = x11;
		if (x11>max_x)
			max_x = x11;
		stringstream geek12(twelve);
		double x12 = 0;
		geek12 >> x12;
		users[user_count][12].x = x12;
		if (x12 <= min_x)
			min_x = x12;
		if (x12>max_x)
			max_x = x12;
		stringstream geek13(thirteen);
		double x13 = 0;
		geek13 >> x13;
		users[user_count][13].x = x13;
		if (x13 <= min_x)
			min_x = x13;
		if (x13>max_x)
			max_x = x13;
		stringstream geek14(fourteen);
		double x14 = 0;
		geek14 >> x14;
		users[user_count][14].x = x14;
		if (x14 <= min_x)
			min_x = x14;
		if (x14>max_x)
			max_x = x14;
		stringstream geek15(fifteen);
		double x15 = 0;
		geek15 >> x15;
		users[user_count][15].x = x15;
		if (x15 <= min_x)
			min_x = x15;
		if (x15>max_x)
			max_x = x15;
		stringstream geek16(sixteen);
		double x16 = 0;
		geek16 >> x16;
		users[user_count][16].x = x16;
		if (x16 <= min_x)
			min_x = x16;
		if (x16>max_x)
			max_x = x16;
		stringstream geek17(seventeen);
		double x17 = 0;
		geek17 >> x17;
		users[user_count][17].x = x17;
		if (x17 <= min_x)
			min_x = x17;
		if (x17>max_x)
			max_x = x17;
		stringstream geek18(eighteen);
		double x18 = 0;
		geek18 >> x18;
		users[user_count][18].x = x18;
		if (x18 <= min_x)
			min_x = x18;
		if (x18>max_x)
			max_x = x18;
		stringstream geek19(nineteen);
		double x19 = 0;
		geek19 >> x19;
		users[user_count][19].x = x19;
		if (x19 <= min_x)
			min_x = x19;
		if (x19>max_x)
			max_x = x19;
		stringstream geek20(twenty);
		double x20 = 0;
		geek20 >> x20;
		users[user_count][20].x = x20;
		if (x20 <= min_x)
			min_x = x20;
		if (x20>max_x)
			max_x = x20;
		user_count++;
	}
	}

	{
	ifstream data("Sample_Data_60+21_y.csv");

	if (!data.is_open()) cout << "ERROR: File Open" << '\n';


	string one;
	string two;
	string three;
	string four;
	string five;
	string six;
	string seven;
	string eight;
	string nine;
	string ten;
	string eleven;
	string twelve;
	string thirteen;
	string fourteen;
	string fifteen;
	string sixteen;
	string seventeen;
	string eighteen;
	string nineteen;
	string twenty;

	int user_count = 0;

	while (data.good() && user_count<user_total)
	{

		getline(data, one, ',');
		getline(data, two, ',');
		getline(data, three, ',');
		getline(data, four, ',');
		getline(data, five, ',');
		getline(data, six, ',');
		getline(data, seven, ',');
		getline(data, eight, ',');
		getline(data, nine, ',');
		getline(data, ten, ',');
		getline(data, eleven, ',');
		getline(data, twelve, ',');
		getline(data, thirteen, ',');
		getline(data, fourteen, ',');
		getline(data, fifteen, ',');
		getline(data, sixteen, ',');
		getline(data, seventeen, ',');
		getline(data, eighteen, ',');
		getline(data, nineteen, ',');
		getline(data, twenty, '\n');


		/*  ---------------------------   */
		stringstream geek(one);
		double y1 = 0;
		geek >> y1;
		users[user_count][1].y = y1;
		if (y1 <= min_y)
			min_y = y1;
		if (y1>max_y)
			max_y = y1;

		stringstream geek2(two);
		double y2 = 0;
		geek2 >> y2;
		users[user_count][2].y = y2;
		if (y2 <= min_y)
			min_y = y2;
		if (y2>max_y)
			max_y = y2;

		stringstream geek3(three);
		double y3 = 0;
		geek3 >> y3;
		users[user_count][3].y = y3;
		if (y3 <= min_y)
			min_y = y3;
		if (y3>max_y)
			max_y = y3;

		stringstream geek4(four);
		double y4 = 0;
		geek4 >> y4;
		users[user_count][4].y = y4;
		if (y4 <= min_y)
			min_y = y4;
		if (y4>max_y)
			max_y = y4;
		stringstream geek5(five);
		double y5 = 0;
		geek5 >> y5;
		users[user_count][5].y = y5;
		if (y5 <= min_y)
			min_y = y5;
		if (y5>max_y)
			max_y = y5;
		stringstream geek6(six);
		double y6 = 0;
		geek6 >> y6;
		users[user_count][6].y = y6;
		if (y6 <= min_y)
			min_y = y6;
		if (y6>max_y)
			max_y = y6;
		stringstream geek7(seven);
		double y7 = 0;
		geek7 >> y7;
		users[user_count][7].y = y7;
		if (y7 <= min_y)
			min_y = y7;
		if (y7>max_y)
			max_y = y7;
		stringstream geek8(eight);
		double y8 = 0;
		geek8 >> y8;
		users[user_count][8].y = y8;
		if (y8 <= min_y)
			min_y = y8;
		if (y8>max_y)
			max_y = y8;
		stringstream geek9(nine);
		double y9 = 0;
		geek9 >> y9;
		users[user_count][9].y = y9;
		if (y9 <= min_y)
			min_y = y9;
		if (y9>max_y)
			max_y = y9;
		stringstream geek10(ten);
		double y10 = 0;
		geek10 >> y10;
		users[user_count][10].y = y10;
		if (y10 <= min_y)
			min_y = y10;
		if (y10>max_y)
			max_y = y10;
		stringstream geek11(eleven);
		double y11 = 0;
		geek11 >> y11;
		users[user_count][11].y = y11;
		if (y11 <= min_y)
			min_y = y11;
		if (y11>max_y)
			max_y = y11;
		stringstream geek12(twelve);
		double y12 = 0;
		geek12 >> y12;
		users[user_count][12].y = y12;
		if (y12 <= min_y)
			min_y = y12;
		if (y12>max_y)
			max_y = y12;
		stringstream geek13(thirteen);
		double y13 = 0;
		geek13 >> y13;
		users[user_count][13].y = y13;
		if (y13 <= min_y)
			min_y = y13;
		if (y13>max_y)
			max_y = y13;
		stringstream geek14(fourteen);
		double y14 = 0;
		geek14 >> y14;
		users[user_count][14].y = y14;
		if (y14 <= min_y)
			min_y = y14;
		if (y14>max_y)
			max_y = y14;
		stringstream geek15(fifteen);
		double y15 = 0;
		geek15 >> y15;
		users[user_count][15].y = y15;
		if (y15 <= min_y)
			min_y = y15;
		if (y15>max_y)
			max_y = y15;
		stringstream geek16(sixteen);
		double y16 = 0;
		geek16 >> y16;
		users[user_count][16].y = y16;
		if (y16 <= min_y)
			min_y = y16;
		if (y16>max_y)
			max_y = y16;
		stringstream geek17(seventeen);
		double y17 = 0;
		geek17 >> y17;
		users[user_count][17].y = y17;
		if (y17 <= min_y)
			min_y = y17;
		if (y17>max_y)
			max_y = y17;
		stringstream geek18(eighteen);
		double y18 = 0;
		geek18 >> y18;
		users[user_count][18].y = y18;
		if (y18 <= min_y)
			min_y = y18;
		if (y18>max_y)
			max_y = y18;
		stringstream geek19(nineteen);
		double y19 = 0;
		geek19 >> y19;
		users[user_count][19].y = y19;
		if (y19 <= min_y)
			min_y = y19;
		if (y19>max_y)
			max_y = y19;
		stringstream geek20(twenty);
		double y20 = 0;
		geek20 >> y20;
		users[user_count][20].y = y20;
		if (y20 <= min_y)
			min_y = y20;
		if (y20>max_y)
			max_y = y20;
		user_count++;
	}
	}

	   
}

void ReadCloudlets()
{
	
		ifstream data("Cloudlets.csv");

		if (!data.is_open()) cout << "ERROR: File Open" << '\n';


		string one;
		string two;
		string three;
		string four;
		
		int user_count = 0;

		while (data.good() && user_count<cloudlet_total)
		{

			getline(data, one, ',');
			getline(data, two, ',');
			getline(data, three, ',');
			getline(data, four, '\n');


			/*  ---------------------------   */
			stringstream geek(one);
			double x = 0;
			geek >> x;

			stringstream geek2(two);
			double y = 0;
			geek2 >> y;


			stringstream geek3(three);
			double p = 0;
			geek3 >> p;

			stringstream geek4(four);
			double b = 0;
			geek4 >> b;

			cloudlets[user_count].x = x;
			cloudlets[user_count].y = y;
			cloudlets[user_count].p = p;
			cloudlets[user_count].b = b;

			user_count++;
		}

}

bool compareTwoApplicationsForBFD(Users user_1, Users user_2)
{
	return user_1.p>=user_2.p &&  user_1.b >= user_2.b;
}


double BFD()
{

	
	double D = 0;


	T = 20; 
	for (int t = 1;t <=T;t++)
	{
		Users* u = new Users[user_total];
		for (int i = 0;i<user_total;i++)
		{
		     users[i][t].index = i;
			 u[i] = users[i][t];
		}
		std::sort(u, u + user_total, compareTwoApplicationsForBFD);
	

		for (int num = 0;num<cloudlet_total;num++)
			C[num] = cloudlets[num];

		for (int i = 0;i < user_total;i++)
		{
			double MAX = INT_MIN;
			int index = -1;
			
			
			/* Find the NEAREST cloudlet whose resources are enough for the current application */
			
			double distance = INT_MAX;
			int index_cloudlet = -1;
			
			for (int j = 0;j < cloudlet_total;j++)
			{
			
			  double CurrentDistance;
			    
			 // CurrentDistance = distance(u[i].x, u[i].y, cloudlets[j].x, cloudlets[j].y);
			  
			  CurrentDistance = sqrt(pow(cloudlets[j].x - u[i].x, 2) + pow(cloudlets[j].y - u[i].y, 2) * 1.0);
			  
			  if(CurrentDistance < distance &&  C[j].p >= u[i].p && C[j].b >= u[i].b)
			  {
			    index = j;
			    distance = CurrentDistance;
			  }
			}
			
			/* end of the new code! */
			
			
			/*for (int j = 0;j < cloudlet_total;j++)
			{
			    
				double fp = (u[i].p ) / (C[j].p + 1);
				double fb = (u[i].b) / (C[j].b + 1);

				if (fp + fb >= MAX && C[j].p >= u[i].p && C[j].b >= u[i].b)
				{
					index = j;
					MAX = fp + fb;
				}
			}*/
			

			R[u[i].index][t] = index;
			C[index].b -= u[i].b;
			C[index].p -= u[i].p;
		}

		/* TEST */
		bool* testarray = new bool[cloudlet_total];
		for (int i = 0;i<cloudlet_total;i++)
			testarray[i] = 0;

		for (int i = 0;i<user_total;i++)
		{
			int index;
			index = R[i][t];
			testarray[index] = 1;

			UsageP[index][t] += users[i][t].p;
			UsageB[index][t] += users[i][t].b;

		}
		int count = 0;
		for (int i = 0;i<cloudlet_total;i++)
			if (testarray[i] == 1)
				count++;

		load_balance += count;

		for (int i = 0;i < user_total;i++)
		{
			int index = R[i][t];
			D += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y)  * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;
			Latency[i] += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b + users[i][t].w / users[i][t].p) * 1000000;

			/* calculate metrics */
			OffloadingTime += distance(users[i][t].x, users[i][t].y, cloudlets[index].x, cloudlets[index].y) * theta + (users[i][t].id / users[i][t].b) * 1000000;
			ComputationTime += (users[i][t].w / users[i][t].p) * 1000000;
		}
		if (t >= 2)
			for (int i = 0;i < user_total;i++)
			{

				int index1 = R[i][t - 1];
				int index2 = R[i][t];

				if (index1 != index2)
				{

					if (i<user_total / 3)
						MigrationNewYork++;
					if (i >= user_total / 3 && i<(user_total * 2) / 3)
						MigrationOrlando++;
					if (i >= (user_total * 2) / 3 && i<user_total)
						MigrationDACT++;

				    Migrations[i]++;
					migration++; D += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					Latency[i] += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;
					MigrationTime += distance(cloudlets[index1].x, cloudlets[index1].y, cloudlets[index2].x, cloudlets[index2].y) * theta + (users[i][t].md / users[i][t].b) * 1000000;

				}
			}
	}
	return D;
}


int main()
{
    ReadReal();
    Results_MatrixCompletion();
	
	ReadCloudlets();

	/* allocating memory for V which will be used in OAMC */
	for (int i = 0;i<user_total;i++)
	{
		V[i] = new int[cloudlet_total];
		r[i] = new double[cloudlet_total];
		Temp_V[i] = new int[cloudlet_total];
		F[i] = new int[cloudlet_total];
		N[i] = new int[cloudlet_total];
	}
	/* allocating memory for S which will be used in Assignment */
	for (int i = 0;i<cloudlet_total;i++)
		{
		  S[i] = new bool[user_total];
		  UsageP[i] = new double[21]; /* UsageP[i][t] shows P which has already consumed from cloudlet i at time slot t*/
		  UsageB[i] = new double[21];
		}

    for(int i=0;i<cloudlet_total;i++)
	  for(int t=1;t<=20;t++)
	    {
		 UsageP[i][t] = 0;
		 UsageB[i][t] = 0;
		}
	
	
	load_balance = 0;
	for (int i = 0;i<cloudlet_total;i++)
		for (int t = 1;t <= 20;t++)
		{
			UsageP[i][t] = 0;
			UsageB[i][t] = 0;
		}

	ComputationTime = 0;
	OffloadingTime = 0;
	MigrationTime = 0;

	Initialize();

	auto start77 = high_resolution_clock::now();

	migration = 0;
	cout << BFD() << " BFD Migration  " << migration << " load balancing " << load_balance << endl;

	cout << "Migrations for different datasets\n" << "New York " << MigrationNewYork << " Orlando " << MigrationOrlando << " DACT " << MigrationDACT << endl;
	MigrationNewYork = MigrationOrlando = MigrationDACT = 0;

	cout << "metrics\n" << endl;
	migration = 0;
	cout << "computation time " << ComputationTime << endl;
	cout << "offloading time " << OffloadingTime << endl;
	cout << "migration time " << MigrationTime << endl;
	ComputationTime = OffloadingTime = MigrationTime = 0;
	// Get ending time point
	auto stop77 = high_resolution_clock::now();

	auto duration77 = duration_cast<microseconds>(stop77 - start77);

	cout << "Time taken by function: "
		<< duration77.count() << " microseconds" << endl;

	double AveragePConsumption = 0;
	double AverageBConsumption = 0;

	for (int t = 1;t <= 20;t++)
	{
		for (int j = 0;j<cloudlet_total;j++)
		{
			AveragePConsumption += (UsageP[j][t]) / cloudlets[j].p;
			AverageBConsumption += (UsageB[j][t]) / cloudlets[j].b;
		}

	}

	cout << "Average p consumption for BFD: " << AveragePConsumption / (20 * cloudlet_total) << endl;
	cout << "Average B consumption for BFD: " << AverageBConsumption / (20 * cloudlet_total) << endl;
	system("pause");

	
    return 0;
}